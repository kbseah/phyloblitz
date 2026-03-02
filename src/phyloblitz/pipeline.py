import json
import logging
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime
from functools import wraps
from multiprocessing import Pool
from pathlib import Path
from random import sample, seed
from subprocess import PIPE, Popen
from tempfile import NamedTemporaryFile, TemporaryDirectory

import matplotlib.pyplot as plt
import numpy as np
import oxli
import pyfastx
import pymarkovclustering as pymcl
import pysam

from phyloblitz.report import (
    generate_histogram,
    generate_report_html,
    generate_report_md,
)
from phyloblitz.utils import (
    CIGAROPS,
    filter_paf_overhang,
    lists_common_prefix,
    parse_cigar_ops,
)

logger = logging.getLogger(__name__)


def merge_intervals(intervals: list) -> list:
    """Merge overlapping numerical intervals.

    Intervals which touch will be merged, e.g. (0, 10), (10, 20) --> (0, 20).

    :param intervals: List of tuples of [start, end) 0-based pythonic,
        non-negative coordinates.
    :returns: List of tuples of intervals with overlaps merged, sorted by start
    :rtype: list
    """
    intervals.sort(key=lambda i: i[0])  # sort in place
    merged = [intervals[0]]
    for curr in intervals[1:]:
        if curr[0] <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(curr[1], merged[-1][1]))
        else:
            merged.append(curr)
    return merged


def sam_seq_generator(
    sam_file: str | Path,
    minlen: int = 1200,
    no_supp: bool = False,
    flanking: int = 0,
) -> dict:
    """Filter SAM alignment to get read segments for all-vs-all mapping.

    This uses pysam.AlignedSegment.query_sequence; if the read is mapped to
    reverse strand, the sequence and coordinates are already
    reverse-complemented.

    Secondary alignments are always skipped, because they represent hits to
    other reference sequences and hence not required for extracting read
    segments. Supplementary alignments potentially represent other copies of
    target gene on the same read, and should be extracted. The -Y flag in
    minimap2 must be set so that supplementary alignments are softclipped, not
    hardclipped, to allow correct extraction of the aligned segment.

    :param sam_file: Path to SAM file
    :param minlen: Minimum query alignment length; adjust if targeting a
        different gene, e.g. LSU rRNA
    :param no_supp: Ignore supplementary alignments
    :param flanking: Additional flanking sequence to extract
    :returns: dict of read name, start, stop, sequence, sequence quality
        scores, leading flanking sequence, leading flanking sequence qualities,
        trailing flanking sequence, trailing flanking sequence qualities.
    :rtype: tuple
    """
    logger.info("Filtering alignment for primary mappings with length >= %d", minlen)
    if flanking > 0:
        logger.info("Including flanking %d bp sequence when extracting reads", flanking)
    sam = pysam.AlignmentFile(sam_file, "r")
    for i in sam.fetch():
        if i.query_alignment_length < minlen or i.query_alignment_sequence is None:
            continue
        if i.is_secondary or (i.is_supplementary and no_supp):
            continue
        if i.is_supplementary and not no_supp:
            logger.debug("Parsing supplementary alignment for read %s", i.query_name)
        mystart = max(0, i.query_alignment_start)
        myend = min(i.query_length, i.query_alignment_end)
        pre_start = max(0, i.query_alignment_start - flanking)
        post_end = min(i.query_length, i.query_alignment_end + flanking)
        yield dict(
            zip(
                [
                    "rname",
                    "start",
                    "end",
                    "seq",
                    "quals",
                    "pre",
                    "pre_quals",
                    "post",
                    "post_quals",
                ],
                [
                    i.query_name,
                    mystart,
                    myend,
                    i.query_sequence[mystart:myend],
                    i.query_qualities_str[mystart:myend],
                    i.query_sequence[pre_start:mystart],
                    i.query_qualities_str[pre_start:mystart],
                    i.query_sequence[myend:post_end],
                    i.query_qualities_str[myend:post_end],
                ],
                strict=True,
            )
        )


def spoa_assemble_fasta(label_fastq: tuple) -> tuple[str, str]:
    """Run spoa assembly on a Fastq input file.

    :param label_fastq: Tuple of input label (str) and path to Fastq file with
        reads to assemble
    :returns: Tuple of input label and alignment of consensus and input
        sequences in Fasta format (stdout from spoa -r 2)
    :rtype: tuple
    """
    label, fastq = label_fastq
    cmd = ["spoa", "-r", "2", fastq]
    logger.debug("spoa command: %s", " ".join([str(i) for i in cmd]))
    proc = Popen(cmd, stdout=PIPE, text=True)
    # TODO directing stderr to PIPE and logger prevents mp.Pool from closing?
    # ignoring stderr from spoa for now
    out = proc.communicate()[0]
    logger.info("Assembly complete for cluster %s", str(label))
    return label, out


def parse_spoa_r2(fasta: str) -> dict:
    """Parse spoa gapped alignment + consensus in Fasta format (-r 2 output).

    :param fasta: Assembly from spoa as list of strings
    :returns: dict of sequences (multiple lines concatenated) keyed by headers.
        The conensus sequence assembled by spoa must have key "Consensus".
    :rtype: dict
    """
    seqs = {}
    prev_hdr = ""
    prev_seq = ""
    for line in fasta.split("\n"):
        if line.startswith(">"):
            if len(prev_hdr) == 0 and len(prev_seq) == 0:
                prev_hdr = line.rstrip()[1:]
            elif len(prev_hdr) > 0 and len(prev_seq) > 0:
                seqs[prev_hdr] = prev_seq
                prev_seq = ""
                prev_hdr = line.rstrip()[1:]
        else:
            prev_seq += line.rstrip()
    # catch last
    if len(prev_hdr) > 0 and len(prev_seq) > 0:
        seqs[prev_hdr] = prev_seq
    return seqs


def count_spoa_aln_vars(seqs: dict) -> dict:
    """Count mismatches and gaps vs consensus for each sequence in a cluster.

    :param seqs: Dict of sequences parsed by parse_spoa_r2; the consensus
        sequence must have key "Consensus"
    :returns: Count of base matches and mismatches, gaps relative to query and
        consensus, leading and trailing query gaps, for each sequence in the
        dict relative to consensus
    :rtype: dict
    """
    var = defaultdict(lambda: defaultdict(int))
    cons = seqs["Consensus"]
    for hdr, seq in seqs.items():
        if hdr == "Consensus":
            continue
        # get coordinates without trailing and leading gaps
        span = re.search(r"^-*([^-].*[^-])-*$", seq).span(1)
        var[hdr]["query_lead_gap"] = span[0]
        var[hdr]["query_trail_gap"] = len(seq) - span[1]
        for i in range(span[0], span[1]):
            if seq[i] == cons[i] and seq[i] != "-":
                # Matching base but not common gap
                var[hdr]["match"] += 1
            elif seq[i] == "-" and cons[i] != "-":
                var[hdr]["query_gap"] += 1
            elif seq[i] != "-" and cons[i] == "-":
                var[hdr]["cons_gap"] += 1
            elif seq[i] != "-" and cons[i] != "-" and seq[i] != cons[i]:
                var[hdr]["mismatch"] += 1
    return var


def count_spoa_aln_persite_vars(seqs: dict) -> dict:
    """[WIP] Calculate entropy per alignment position for clustered sequences vs consensus.

    If mismatches/gaps between sequences and the consensus are solely due to
    technical sequencing error, rather than erroneous clustering of divergent
    underlying reads, we expect fraction of variants per site to be roughly the
    sequencing error. If for example two sequences are misclustered, then we
    should see a secondary peak of ~50% variant coverage.

    :param seqs: Dict of aligned sequences parsed by parse_spoa_r2; the
        consensus sequence must have key "Consensus"
    :returns: dict keyed by alignment position, value is the sequence entropy
        per site in bits.
    """
    var = {}
    for i in range(len(seqs["Consensus"])):
        column = [seqs[hdr][i] for hdr in seqs if hdr != "Consensus"]
        _values, counts = np.unique(column, return_counts=True)
        counts_norm = counts / counts.sum()
        var[i] = -(counts_norm * np.log(counts_norm) / np.log(2)).sum()
    return var


class Pipeline:
    """phyloblitz pipeline run object."""

    def __init__(self, args: dict) -> None:
        """Construct Pipeline object.

        Attributes:
        * _ref : Path to reference database Fasta file
        * _refindex : Optional path to reference database index
        * _reads : Path to input reads in Fastq or Fastq.gz format
        * _outdir : Output folder path
        * _prefix : Filename prefix for output files
        * _platform : Sequencing mode, passed as mapping preset to minimap2
        * _resume : If True, resume partially executed run
        * _stats : Dict with run metadata and processed results, for
            troubleshooting and comparison of runs

        """
        self._ref = args["db"]
        self._refindex = args["dbindex"]
        self._reads = args["reads"]
        self._outdir = args["outdir"]
        self._prefix = args["prefix"]
        self._platform = args["platform"]
        self._resume = args["resume"]
        self._stats = {}
        self._stats["args"] = args
        self._stats["runstats"] = {}
        self._stats["starttime"] = str(datetime.now())

    OUTFILE_SUFFIX = {
        "initial_map": "_minimap_initial.sam",
        "intervals_fastq": "_intervals.fastq",
        "mapped_segments": "_mapped.fastq",
        "ava_map": "_ava.paf",
        "ava_filter": "_ava_filter.paf",
        "ava_seqtab": "_ava_seq.tab",
        "mcl_cluster": "_mcl.out",
        "isonclust3_cluster": "_isonclust3_out/clustering/final_clusters.tsv",
        "cluster_asm": "_final.fasta",
        "cluster_tophits": "_final_tophits.paf",
        "report_json": "_report.json",
        "report_md": "_report.md",
        "report_html": "_report.html",
        "report_dvs_hist": "_report_dvs_hist.png",
        "report_kmercount_plot": "_report_kmercount_plot.png",
    }

    def check_run_file(self, stage: str) -> bool:
        """Check if intermediate output file has been created.

        :param stage: Name of run stage, must be a key of self.OUTFILE_SUFFIX
        :returns: True if file already exists at expected path
        :rtype: bool
        """
        return Path.is_file(self.pathto(stage))

    def pathto(self, stage: str, basename_only: bool = False) -> Path:
        """Combine output directory prefix and filenames to intermediate file path.

        :param stage: Name of run stage, must be a key of self.OUTFILE_SUFFIX
        :param basename_only: Only report the base filename if True
        :returns: Expected path to intermediate output file
        :rtype: Path
        """
        try:
            if basename_only:
                return Path(Path(self._prefix + self.OUTFILE_SUFFIX[stage]).name)
            return Path(self._outdir) / Path(self._prefix + self.OUTFILE_SUFFIX[stage])
        except KeyError as e:
            e.add_not(f"Unknown intermediate file {stage}")
            raise

    def check_stage_file(stage: str, message: str):
        """Check whether outputs for each stage of Pipeline already exist.

        Output files are checked by their expected filenames. If the expected
        outputs already exist, the stage is skipped entirely. Use this function
        as a decorator for individual stage methods. Enables resuming the
        pipeline from partial output.

        :param stage: Name of run stage, must be a key of self.OUTFILE_SUFFIX
        :param message: Logging message to emit on starting each stage
        """

        # TODO specify more than one output file; check required input files
        def check_stage_decorator(func):
            @wraps(func)
            def wrapped_function(self, *args, **kwargs):
                if self.check_run_file(stage):
                    if self._resume:
                        logger.info(
                            "Stage %s file output already present, skipping",
                            stage,
                        )
                        return None
                    logger.error(
                        "Stage %s file output already present and option --resume not used, exiting",
                        stage,
                    )
                    sys.exit()
                else:
                    logger.info(message)
                    return func(self, *args, **kwargs)

            return wrapped_function

        return check_stage_decorator

    @check_stage_file(
        stage="initial_map",
        message="Initial mapping of reads to identify target intervals",
    )
    def run_minimap(self, threads=12, mode="map-ont", sample=None) -> int:
        """Map reads to reference database with minimap2.

        This initial mapping is used to extract the marker sequence segments
        from input long reads, and to summarize taxonomic composition of
        library based on taxonomy of the reference hits.

        Use mode presets for minimap2. Do not report unmapped reads or
        secondary mappings. Report supplementary mappings as soft-clipped
        (default hard-clipped) so the flanking sequence can be correctly parsed
        by pysam.

        :param threads: Number of threads for minimap2 to use
        :param mode: Mapping preset for minimap2, one of: `map-ont`, `map-pb`,
            "lr:hq", "map-hifi"
        :param sample: Number of reads to sample; if None then use all reads;
            also use this number as the random seed
        :returns: return code for samtools
        :rtype: tuple
        """
        # Don't use mappy because it doesn't support all-vs-all yet
        infile = self._reads
        if sample is not None:
            logger.info("Taking sample of %d reads", sample)
            # TODO use pyfastx Python API instead of CLI
            temp_infile = NamedTemporaryFile()  # TODO keeptmp
            cmd = [
                "pyfastx",
                "sample",
                "--sequential-read",
                "-n",
                str(sample),
                "-s",
                str(sample),
                "-o",
                temp_infile.name,
                self._reads,
            ]
            proc = Popen(cmd, stderr=PIPE)
            proc.wait()
            infile = temp_infile.name
            logger.debug("Sampled reads in temporary file %s", temp_infile.name)

        logger.info("Mapping reads to reference database with minimap2")
        with Path.open(self.pathto("initial_map"), "w") as sam_fh:
            cmd1 = [
                "minimap2",
                "-ax",
                mode,
                "--sam-hit-only",
                "--secondary=no",
                "-Y",
                f"-t {threads!s}",
            ]
            (
                cmd1.extend([self._refindex, infile])
                if self._refindex is not None
                else cmd1.extend([self._ref, infile])
            )
            logger.debug("minimap command: %s", " ".join([str(i) for i in cmd1]))
            proc1 = Popen(cmd1, stdout=sam_fh, stderr=PIPE)
            nreads = 0
            for l in proc1.stderr:
                nreads_s = re.search(r"mapped (\d+) sequences", l.decode())
                if nreads_s:
                    nreads += int(nreads_s.group(1))
                logger.debug("  minimap log: %s", l.decode().rstrip())
            self._stats["runstats"].update({"total input reads": nreads})
            return proc1.wait()

    @check_stage_file(
        stage="mapped_segments",
        message="Extracting read segments for all-vs-all mapping",
    )
    def extract_reads_for_ava(
        self,
        align_minlen: int = 1200,
        no_supp: bool = False,
        flanking: int = 0,
    ) -> None:
        """Extract read segments aligning to markers for all-vs-all mapping.

        Segments aligning to marker database are extracted in Fastq format to a
        new file, for the next clustering steps. Flanking sequences are also
        extracted, not to the Fastq file but Pipeline._stats["flanking"], so
        they do not influence the marker clustering and assembly.

        :param align_minlen: Minimum alignment length in bp
        :param no_supp: Ignore supplementary alignments if True
        :param flanking: Length of flanking sequence to extract, in bp
        """
        sam_file = self.pathto("initial_map")
        counter = 0
        self._stats["flanking"] = {}
        self._stats["reads2segment"] = {}
        with Path.open(self.pathto("mapped_segments"), "w") as fq_fh:
            for rec in sam_seq_generator(
                sam_file,
                minlen=align_minlen,
                no_supp=no_supp,
                flanking=flanking,
            ):
                name = ":".join([str(rec[i]) for i in ["rname", "start", "end"]])
                if rec["rname"] not in self._stats["reads2segment"]:
                    # Update reads2segment
                    self._stats["reads2segment"][rec["rname"]] = [
                        (rec["start"], rec["end"], name),
                    ]
                else:
                    # Check whether overlapping segment is already present
                    merged = [
                        (
                            other_name,
                            merge_intervals([(s, e), (rec["start"], rec["end"])]),
                        )
                        for (s, e, other_name) in self._stats["reads2segment"][
                            rec["rname"]
                        ]
                    ]
                    merged = [i for i in merged if len(i[1]) == 1]
                    if len(merged) > 0:
                        logger.debug(
                            "Segment %s overlaps with: %s ... skipping...",
                            name,
                            ", ".join([i[0] for i in merged]),
                        )
                        continue
                counter += 1
                fq_fh.write(
                    "@" + name + "\n" + rec["seq"] + "\n+\n" + rec["quals"] + "\n",
                )
                self._stats["flanking"][name] = {
                    i: rec[i] for i in ["pre", "pre_quals", "post", "post_quals"]
                }
        self._stats["runstats"].update({"mapped pass filter": counter})
        logger.info("Read segments extracted for all-vs-all mapping: %d", counter)

    @check_stage_file(
        stage="ava_map",
        message="All-vs-all mapping of mapped reads with minimap2",
    )
    def ava_map(self, mode: str = "map-ont", threads: int = 12) -> int:
        """All-vs-all mapping with minimap2 to generate clusters for assembly.

        Minimap2 presets for Nanopore and PacBio CLR are applied. For Nanopore
        Q20+ and PacBio HiFi/CCS reads, the "ava" presets were modified to
        reflect their lower sequence error. PAF alignment written to file.

        :param mode: Mapping preset mode used for initial minimap2 mapping step
        :param threads: Number of threads for minimap2 to use
        :returns: Return code of the minimap2 process
        :rtype: int
        """
        presets = {
            "map-ont": ["-x", "ava-ont"],
            "map-pb": ["-x", "ava-pb"],
            "lr:hq": ["-x", "lr:hq", "-Xw5", "-e0", "-m100", "-r2k"],
            "map-hifi": ["-x", "map-hifi", "-Xw5", "-e0", "-m100"],
        }
        with Path.open(self.pathto("ava_map"), "w") as paf_fh:
            cmd = (
                [
                    "minimap2",
                ]
                + presets[mode]
                + [
                    "-t",
                    str(threads),
                    self.pathto("mapped_segments"),
                    self.pathto("mapped_segments"),
                ]
            )
            logger.debug("minimap command: %s", " ".join([str(i) for i in cmd]))
            proc = Popen(cmd, stdout=paf_fh, stderr=PIPE)
            for l in proc.stderr:
                logger.debug("  minimap log: %s", l.decode().rstrip())
            return proc.wait()

    @check_stage_file(
        stage="ava_filter",
        message="Filtering overlapping incompatible overhangs from all-vs-all mappings",
    )
    def paf_file_filter_overhangs(self, max_overhang_frac=0.05) -> None:
        """Remove all-vs-all alignments with incompatible overhangs.

        All-vs-all alignments of extracted marker read segments will be
        clustered and used to assemble consensus marker sequences. However the
        alignments are not always end-to-end, but may have soft-clipped
        overhangs, when two sequences share a conserved middle portion but are
        divergent elsewhere. This leads to over-clustering, so such alignments
        with incompatible overhangs must be filtered out. Filtered alignments
        are written to a new file.

        :param max_overhang_frac: Max fraction of read length that same-side
            overhang is allowed to be, passed to filter_paf_overhang()
        """
        logger.info("Max overhang fraction: %f", max_overhang_frac)
        counter_all = 0
        counter = 0
        with (
            Path.open(self.pathto("ava_map")) as fh_in,
            Path.open(self.pathto("ava_filter"), "w") as fh_out,
        ):
            for line in fh_in:
                counter_all += 1
                out = filter_paf_overhang(line, max_overhang_frac=max_overhang_frac)
                if out is not None:
                    counter += 1
                    fh_out.write(out)
        logger.info("Retained %d out of %d all-vs-all alignments", counter, counter_all)

    def paf_get_dvs(self) -> None:
        """Read and summarize the per-base divergence scores from PAF alignments.

        Minimap2 reports per-base divergence for each PAF alignment in the dv:
        tag. We analyze the distribution of dv values from all-vs-all mappings
        as a reference-free measure of sequence error, and also use it to set
        the cutoffs for clustering with the mcl method. Updates Pipeline._stats
        with a dict of dv values per read ("dvs"), and a list of minimum dv
        value for each read ("min_dvs") to plot dv distribution.
        """
        logger.info("Reading dvs data from ava mapping")
        dvs = defaultdict(list)
        with Path.open(self.pathto("ava_filter")) as fh:
            for line in fh:
                spl = line.rstrip().split("\t")
                query = spl[0]
                dv = re.findall(r"dv:f:([\d\.]+)", line)
                if len(dv) != 1:
                    msg = f"Problem in PAF file {str(self.pathto('ava_filter'))}: more than one dv tag in entry"
                    raise ValueError(msg)
                dvs[query].append(float(dv[0]))
        self._stats.update({"dvs": dvs})
        # min dv for each read, to get dv distribution and median to
        # automatically set dv threshold for clustering
        min_dvs = [min(dvs[r]) for r in dvs]
        self._stats.update({"min_dvs": min_dvs})
        self._stats["runstats"].update(
            {
                "ava min dvs median": f"{np.median(min_dvs):.4f}",
            },
        )

    @check_stage_file(
        stage="mcl_cluster",
        message="Clustering with mcl",
    )
    def pymcl_cluster(
        self, dv_max: float = 0.03, dv_max_auto: bool = False, inflation: float = 2
    ) -> None:
        """Cluster marker read segments with MCL algorithm.

        Cluster read segments with the MCL algorithm, using one minus the
        sequence divergences from the (filtered) all-vs-all alignments as the
        edge score. The threshold divergence value can be hardcoded, or
        adjusted to the observed distribution, to account for different
        sequencing errors across libraries and platforms. The default inflation
        parameter was empirically chosen. Cluster output are files.

        :param dv_max: Maximum sequence divergence (dv:f tag in PAF file) to accept
        :param dv_max_auto: Set dv_max to max(0.001, 2 * median all-vs-all
            sequence divergence); parameter `dv_max` will be ignored). The value
            0.001 is to account for cases where median divergence is zero.
        :param inflation: Inflation parameter for MCL
        """
        if dv_max_auto:
            dv_max = max(
                0.001,
                2 * float(self._stats["runstats"]["ava min dvs median"]),
            )
        edges = []
        with Path.open(self.pathto("ava_filter")) as fh_paf:
            for line in fh_paf:
                spl = line.rstrip().split("\t")
                query = spl[0]
                target = spl[5]
                dv = re.findall(r"dv:f:([\d\.]+)", line)
                if len(dv) == 1 and float(dv[0]) < dv_max:
                    score = 1 - float(dv[0])
                    edges.append((query, target, score))
        self._stats["runstats"].update({"dv_max applied": dv_max})
        logger.info("Inflation parameter %f", inflation)
        matrix, labels = pymcl.edges_to_sparse_matrix(edges)
        mcl_matrix = pymcl.mcl(matrix, inflation=inflation)
        clusters = pymcl.extract_clusters(mcl_matrix, labels)
        with Path.open(self.pathto("mcl_cluster"), "w") as fh:
            fh.writelines("\t".join(c) + "\n" for c in clusters)

    @check_stage_file(
        stage="isonclust3_cluster",
        message="Clustering with isonclust3",
    )
    def isonclust3_cluster(self) -> int:
        """Cluster marker read segments with isONclust3.

        isONclust3 was originally designed for long-read transcriptome
        datasets, and uses minimizers with iterative cluster merging. The
        default preset for Nanopore, "ont", is applied here for both legacy and
        Q20+ Nanopore reads, while the "pacbio" preset is used for both PacBio
        CLR and HiFi. Cluster output is a directory with multiple files.

        :returns: Return code of the isONclust3 process.
        :rtype: int
        """
        if self._platform in ["lr:hq", "map-ont"]:
            isonclust3_mode = "ont"
        elif self._platform in ["map-hifi", "map-pb"]:
            isonclust3_mode = "pacbio"
        # isonclust3 takes output folder path as argument, automatically
        # creates `clustering` subfolder, so go two levels up
        outfolder = self.pathto("isonclust3_cluster").parent.parent
        cmd = [
            "isONclust3",
            "--no-fastq",
            "--fastq",
            self.pathto("mapped_segments"),
            "--mode",
            isonclust3_mode,
            "--outfolder",
            outfolder,
        ]
        if isonclust3_mode == "ont":
            cmd.append("--post-cluster")
        logger.debug("isonclust3 command: %s", " ".join([str(i) for i in cmd]))
        proc = Popen(cmd, stdout=PIPE)
        for l in proc.stdout:
            logger.debug("  isonclust3 log: %s", l.decode().rstrip())
        return proc.wait()

    def _cluster_seqs_from_mcl(
        self, mcl_out: Path, reads: Path, keeptmp: bool, min_clust_size: int = 5
    ) -> tuple:
        """Extract sequences from each MCL cluster to Fastq files.

        :param mcl_out: Path MCL output file
        :param reads: Path to reads from extract_reads_for_ava step
        :param keeptmp: Keep temporary files?
        :param min_clust_size: Minimum cluster size to return Fastq file
        :returns: Dict of file handles to each Fastq file keyed by cluster ID
        :returns: Dict of sequence IDs keyed by cluster ID
        """
        counter = 0
        fastq_handles = {}
        seq2cluster = {}
        with Path.open(mcl_out) as fh:
            for line in fh:
                seqs = line.rstrip().split("\t")
                logger.info("Cluster %d comprises %d sequences", counter, len(seqs))
                for seq in seqs:
                    seq2cluster[seq] = counter  # assume each seq in only one cluster
                if len(seqs) >= min_clust_size:
                    fastq_handles[counter] = NamedTemporaryFile(
                        suffix=".fastq",
                        mode="w",
                        delete=(not keeptmp),
                        delete_on_close=False,
                    )
                else:
                    logger.debug("Cluster %d below size cutoff", counter)
                counter += 1
        for name, seq, qual in pyfastx.Fastx(reads):
            if name in seq2cluster and seq2cluster[name] in fastq_handles:
                fastq_rec = "@" + name + "\n" + seq + "\n+\n" + qual + "\n"
                fastq_handles[seq2cluster[name]].write(fastq_rec)
        cluster2seq = defaultdict(list)
        for seq, cluster in seq2cluster.items():
            cluster2seq[cluster].append(seq)
        for handle in fastq_handles.values():
            handle.close()
        return fastq_handles, cluster2seq

    def _cluster_seqs_from_isonclust3(
        self, isonclust3_out: Path, reads: Path, keeptmp: bool, min_clust_size: int = 5
    ) -> tuple:
        """Extract sequences from each isONclust3 cluster to Fastq files.

        The numbering of clusters by isONclust3 itself is not consistent
        between the final_clusters.tsv file and the output fastq files; best to
        extract the sequences ourselves.

        :param isonclust3_out: Path to isONclust3 output file
        :param reads: Path to reads from extract_reads_for_ava step
        :param keeptmp: Keep temporary files?
        :param min_clust_size: Minimum cluster size to return Fastq file
        :returns: Dict of file handles to each Fastq file keyed by cluster ID
        :returns: Dict of sequence IDs keyed by cluster ID
        """
        fastq_handles = {}
        seq2cluster = {}
        with Path.open(isonclust3_out) as fh:
            for line in fh:
                (clust, seqname) = line.rstrip().split("\t")
                seq2cluster[seqname] = clust
        cluster2seq = defaultdict(list)
        for seqname in seq2cluster:
            cluster2seq[seq2cluster[seqname]].append(seqname)
        for cluster in cluster2seq:
            logger.info(
                "Cluster %s comprises %d sequences", cluster, len(cluster2seq[cluster])
            )
            if len(cluster2seq[cluster]) >= min_clust_size:
                fastq_handles[cluster] = NamedTemporaryFile(
                    suffix=".fastq",
                    mode="w",
                    delete=(not keeptmp),
                    delete_on_close=False,
                )
            else:
                logger.debug("Cluster %s below size cutoff", cluster)
        for seqname, seq, qual in pyfastx.Fastx(reads):
            if seqname in seq2cluster and seq2cluster[seqname] in fastq_handles:
                fastq_rec = "@" + seqname + "\n" + seq + "\n+\n" + qual + "\n"
                fastq_handles[seq2cluster[seqname]].write(fastq_rec)
        for handle in fastq_handles.values():
            handle.close()
        return fastq_handles, cluster2seq

    @check_stage_file(
        stage="cluster_asm",
        message="Extract cluster sequences and assemble with spoa",
    )
    def assemble_clusters(
        self,
        cluster_tool: str = "isonclust3",
        threads: int = 12,
        rseed: int = 12345,
        keeptmp: bool = False,
        min_clust_size: int = 5,
        max_clust_size: int = 500,
    ) -> None:
        """Extract cluster sequences and assemble with spoa

        Cluster reads with specified clustering tool, extract read segments for
        each cluster to separate Fastq files, and assemble consensus sequence
        for each with Spoa. Assembly step is embarrassingly parallelized.

        Updates Pipeline._stats["cluster2seq"] with lists of read names keyed
        by cluster id. Some summary stats on clusters are written to
        Pipeline._stats["runstats"].

        Other per-cluster summaries in Pipeline._stats: "cluster variant
        counts" -- number of variant sites by type for each read vs. its
        cluster consensus; "cluster persite variant counts" -- per-site entropy
        vs. consensus [WIP]; "cluster cons parsed" -- cluster alignment
        including consensus, produced by spoa.

        :param cluster_tool: Clustering method used, either "mcl" or "isonclust3".
        :param threads: Number of parallel assembly jobs to run.
        :param rseed: Random seed for downsampling reads in large clusters.
        :param keeptmp: If True, do not delete Fastq files with extracted reads.
        :param min_clust_size: Only assemble clusters containing at least this
            number of reads.
        :param max_clust_size: Downsample reads for clusters above this size.
        """
        if cluster_tool == "mcl":
            fastq_handles, cluster2seq = self._cluster_seqs_from_mcl(
                self.pathto("mcl_cluster"),
                self.pathto("mapped_segments"),
                keeptmp=keeptmp,
                min_clust_size=min_clust_size,
            )
        elif cluster_tool == "isonclust3":
            fastq_handles, cluster2seq = self._cluster_seqs_from_isonclust3(
                self.pathto("isonclust3_cluster"),
                self.pathto("mapped_segments"),
                keeptmp=keeptmp,
                min_clust_size=min_clust_size,
            )
        self._stats.update({"cluster2seq": cluster2seq})
        self._stats["runstats"].update(
            {
                "number of clusters": len(cluster2seq),
                "number of clusters > 5 reads": len(
                    [i for i in cluster2seq if len(cluster2seq[i]) > 5]
                ),
                "total reads in clusters": sum(
                    [len(cluster2seq[c]) for c in cluster2seq]
                ),
            },
        )
        # Downsample if >500 sequences in cluster
        logger.debug("Random seed %d", rseed)
        seed(rseed)
        for c, seqs in cluster2seq.items():
            if len(seqs) > max_clust_size:
                logger.debug("Cluster %s has %d reads, downsampling ...", c, len(seqs))
                logger.debug("Downsampling from file %s", fastq_handles[c].name)
                fq = pyfastx.Fastq(fastq_handles[c].name)
                logger.debug("Opened file handle %s", str(fq))
                selected_idx = sorted(sample(range(len(fq)), k=max_clust_size))
                # Do not use a context manager here because we need file later
                newhandle = NamedTemporaryFile(
                    suffix=".fastq",
                    mode="w",
                    delete=(not keeptmp),
                    delete_on_close=False,
                )
                logger.debug("Created new temporary file %s", str(newhandle.name))
                for idx in selected_idx:
                    newhandle.write(fq[idx].raw)
                fastq_handles[c] = newhandle
                newhandle.close()
        logger.info("Assemble consensus from clustered sequences with spoa")
        with Pool(threads) as pool:
            cluster_cons_tuples = pool.map(
                spoa_assemble_fasta,
                [(c, handle.name) for c, handle in fastq_handles.items()],
            )

        for handle in fastq_handles.values():  # close NamedTemporaryFile handles
            handle.close()

        cluster_cons = dict(cluster_cons_tuples)

        cluster_cons_parsed = {i: parse_spoa_r2(cluster_cons[i]) for i in cluster_cons}
        cluster_variant_counts = {
            i: count_spoa_aln_vars(cluster_cons_parsed[i]) for i in cluster_cons_parsed
        }
        cluster_persite_variant_counts = {  # TODO WIP
            i: count_spoa_aln_persite_vars(cluster_cons_parsed[i])
            for i in cluster_cons_parsed
        }
        self._stats.update(
            {
                "cluster variant counts": cluster_variant_counts,
                "cluster persite variant counts": cluster_persite_variant_counts,  # TODO WIP
                "cluster cons parsed": cluster_cons_parsed,
            },
        )
        with Path.open(self.pathto("cluster_asm"), "w") as fh:
            fh.writelines(
                f">cluster_{c!s} Consensus\n"
                + cluster_cons_parsed[c]["Consensus"].replace("-", "")
                + "\n"
                for c in cluster_cons_parsed
            )
            logger.info("Assembled sequences written to %s", self.pathto("cluster_asm"))

    def cluster_flanking_isonclust3(self) -> None:
        """Cluster flanking sequences with isonclust3 and report number of groups.

        rRNA marker genes are more conserved than most genes and intergenic
        sequences. Different strains of the same species may therefore have
        similar or even identical marker sequence, and end up represented in
        the same cluster. As a measure of the potential diversity in each
        cluster, sequence flanking the conserved markers is also extracted and
        clustered, and the resulting number of clusters is reported. Stored in
        Pipeline._stats["cluster flanking numclust"]
        """
        logger.info("Clustering flanking sequences with isonclust3")
        if self._platform in ["lr:hq", "map-ont"]:
            isonclust3_mode = "ont"
        elif self._platform in ["map-hifi", "map-pb"]:
            isonclust3_mode = "pacbio"
        cluster_flanking_numclust = {}
        # Iterate clusters and count number of clustered flanking seqs
        for cluster in self._stats["cluster2seq"]:
            with NamedTemporaryFile(suffix=".fastq", mode="w") as fq:
                for seqid in self._stats["cluster2seq"][cluster]:
                    seq = (
                        self._stats["flanking"][seqid]["pre"]
                        + "NNNNNNNNNN"
                        + self._stats["flanking"][seqid]["post"]
                    )
                    quals = (
                        self._stats["flanking"][seqid]["pre_quals"]
                        + "@@@@@@@@@@"
                        + self._stats["flanking"][seqid]["post_quals"]
                    )
                    fq.write("@" + seqid + "\n" + seq + "\n+\n" + quals + "\n")
                with TemporaryDirectory() as tmpdir:
                    cmd = [
                        "isONclust3",
                        "--no-fastq",
                        "--fastq",
                        fq.name,
                        "--mode",
                        isonclust3_mode,
                        "--outfolder",
                        tmpdir,
                    ]
                    if isonclust3_mode == "ont":
                        cmd.append("--post-cluster")
                    logger.debug(
                        "isonclust3 command: %s", " ".join([str(i) for i in cmd])
                    )
                    proc = Popen(cmd, stdout=PIPE, stderr=PIPE, text=True)
                    proc_stdout, _proc_stderr = proc.communicate()
                    num_clusters = re.search(
                        r"(\d+) different clusters identified",
                        proc_stdout,
                    )
                    if num_clusters:
                        logger.debug(
                            "Flanking sequences for target cluster %s fall into %d clusters",
                            cluster,
                            int(num_clusters.group(1)),
                        )
                        cluster_flanking_numclust[cluster] = int(num_clusters.group(1))
        self._stats["cluster flanking numclust"] = cluster_flanking_numclust

    def cluster_flanking_kmercount(self, k: int = 11, minlen: int = 500) -> None:
        """Count kmers from flanking sequences per cluster.

        As an alternative to clustering flanking sequences, report the k-mer
        multiplicity to be visualized with cluster_flanking_kmercount_plot().
        Stored in Pipeline._stats["flanking_kmer_histo"]

        :param k: K-mer length, in bp
        :param minlen: Only count k-mers for flanking sequence at least this length
        """
        logger.info("Counting k-mers from flanking sequences per cluster")
        kmer_tables = {}
        for cluster in self._stats["cluster2seq"]:
            kmer_tables[cluster] = oxli.KmerCountTable(k)
            for seqid in self._stats["cluster2seq"][cluster]:
                if seqid in self._stats["flanking"]:
                    if len(self._stats["flanking"][seqid]["pre"]) > minlen:
                        kmer_tables[cluster].consume(
                            self._stats["flanking"][seqid]["pre"],
                        )
                    if len(self._stats["flanking"][seqid]["post"]) > minlen:
                        kmer_tables[cluster].consume(
                            self._stats["flanking"][seqid]["post"],
                        )
        self._stats["flanking kmer histo"] = {
            cluster: kmer_table.histo() for cluster, kmer_table in kmer_tables.items()
        }

    @check_stage_file(
        stage="report_kmercount_plot",
        message="Generate k-mer multiplicity plot for report",
    )
    def cluster_flanking_kmercount_plot(self, min_clust_size: int = 5) -> None:
        """Plot k-mer multiplicity of flanking sequences for report.

        The k-mer multiplicity plot gives a visual impression of the sequence
        error and sequence depth for each cluster, and together with other
        metrics can help in judging whether there is adequate coverage for
        assembling the genome(s) represented by a cluster, in addition to the
        potential strain diversity. Image written to file.

        :param min_clust_size: Only report for clusters containing at least this
            number of reads.
        """
        cluster_ids = {
            i: len(self._stats["cluster2seq"][i])
            for i in self._stats["cluster2seq"]
            if len(self._stats["cluster2seq"][i]) > min_clust_size
        }
        histos = {
            i: self._stats["flanking kmer histo"][i]
            for i in self._stats["flanking kmer histo"]
            if i in cluster_ids
        }
        max_kmer_cov = max([max([j[0] for j in histos[i]]) for i in histos])
        fig_height = len(cluster_ids) * 0.6
        fig, axs = plt.subplots(
            len(histos),
            sharex=True,
            sharey=False,
            figsize=(4, fig_height),
            layout="constrained",
        )
        for idx, (cid, histo) in enumerate(histos.items()):
            axs[idx].bar(
                x=[i[0] for i in histo],
                height=[np.log1p(i[1]) for i in histo],
                width=1,
            )
            axs[idx].set_title(cid, y=1.0, pad=-15, loc="right")
            axs[idx].set_xlim(-2, max_kmer_cov + 2)
        fig.savefig(self.pathto("report_kmercount_plot"))

    def db_taxonomy(self) -> None:
        """Get taxonomy string from SILVA headers in database Fasta file.

        Parse SILVA-style headers to get dict of taxonomy strings keyed by
        accession, stored in the Pipeline._acc2tax attribute.
        """
        logger.info("Reading taxonomy from SILVA database file")
        self._acc2tax = {}
        with Path.open(self._ref) as fh:
            for line in fh:
                if line.startswith(">"):
                    spl = line.lstrip(">").rstrip().split(" ")
                    acc = spl[0]
                    taxstring = " ".join(spl[1:]).split(";")
                    self._acc2tax[acc] = taxstring
        logger.debug(" Accessions read: %d", len(self._acc2tax))

    def _per_read_consensus_taxonomy(
        self, sam_file: str | Path, minlen: int = 1200
    ) -> dict:
        """Consensus taxonomy of a single read from initial mapping with minimap2.

        :param sam_file: Path to SAM file of initial mapping
        :param minlen: Only consider alignments where query_alignment_length is
            this value and above.
        :returns: dict with consensus taxonomy (produced by
            lists_common_prefix()) keyed by accession.
        :rtype: dict
        """
        all_taxstrings = defaultdict(list)
        sam = pysam.AlignmentFile(sam_file, "r")
        for i in sam.fetch():
            if i.query_alignment_length >= minlen:
                all_taxstrings[i.query_name].append(self._acc2tax[i.reference_name])
        return {acc: lists_common_prefix(all_taxstrings[acc]) for acc in all_taxstrings}

    def summarize_initial_mapping_taxonomy(
        self, minlen: int, taxlevel: int = 4
    ) -> None:
        """Summarize taxonomy at a specified taxonomy level from initial mapping.

        Updates Pipeline._stats["initial_taxonomy"] with a Counter of the
        observed taxonomy strings at the requested taxon level.

        :param minlen: Only consider alignments where query_alignment_length is
            this value and above.
        :param taxlevel: Taxon level to generate summary (1-based, where 1
            represents the highest taxon level).
        """
        sam_file = self.pathto("initial_map")
        logger.info("Summarizing taxonomic composition of initial mapping")
        common_taxstrings = self._per_read_consensus_taxonomy(sam_file, minlen)
        self._stats.update(
            {
                "initial_taxonomy": Counter(
                    [";".join(i[0:taxlevel]) for i in common_taxstrings.values()],
                ),
            },
        )

    @check_stage_file(
        stage="cluster_tophits",
        message="Mapping assembled sequences to reference database with minimap2",
    )
    def cluster_asm_tophits(self, threads: int = 12) -> int:
        """Map assembled sequences to reference DB and get top hits.

        Find best reference match for cluster consensus sequences by aligning
        with minimap2. Alignment is written to file.

        :param threads: Number of threads for minimap2 to use
        :returns: Return code of the minimap2 process
        """
        cmd = [
            "minimap2",
            "-x",
            "asm5",
            "-c",
            "--eqx",
            "--secondary=no",
            "-t",
            str(threads),
            "-o",
            self.pathto("cluster_tophits"),
        ]
        (
            cmd.extend([self._refindex, self.pathto("cluster_asm")])
            if self._refindex is not None
            else cmd.extend([self._ref, self.pathto("cluster_asm")])
        )
        logger.debug("minimap command: %s", " ".join([str(i) for i in cmd]))
        proc = Popen(cmd, stderr=PIPE)
        for l in proc.stderr:
            logger.debug("  minimap log: %s", l.decode().rstrip())
        return proc.wait()

    def summarize_tophit_paf(self) -> None:
        """Summarize top hits of assembled seqs mapped to SILVA database by minimap2.

        Update Pipeline._stats with a dict of summary stats "cluster_tophits"
        for each hit, keyed by query sequence name.
        """
        out = {}

        with Path.open(self.pathto("cluster_tophits")) as fh:
            for line in fh:
                spl = line.rstrip().split("\t")
                hits = dict(
                    zip(
                        [
                            "qname",
                            "qlen",
                            "qstart",
                            "qend",
                            "strand",
                            "tname",
                            "tlen",
                            "tstart",
                            "tend",
                            "alnmatch",
                            "alnlen",
                        ],
                        spl[0:11],
                        strict=False,
                    )
                )
                cigar = next(i for i in spl if i.startswith("cg:Z:"))
                # Calculate derived metrics from PAF fields
                cigar_summary = parse_cigar_ops(cigar[5:])
                hits.update({CIGAROPS[c]: cigar_summary[c] for c in cigar_summary})
                # remove redundant % sign for display
                hits["align %id"] = "{:.2%}".format(
                    int(hits["alnmatch"]) / int(hits["alnlen"]),
                ).rstrip("%")
                hits["query %aln"] = "{:.2%}".format(
                    (int(hits["qend"]) - int(hits["qstart"])) / int(hits["qlen"]),
                ).rstrip("%")
                hits["target %aln"] = "{:.2%}".format(
                    (int(hits["tend"]) - int(hits["tstart"])) / int(hits["tlen"]),
                ).rstrip("%")
                out[spl[0]] = hits

        # Taxonomy of hit targets
        for rec in out.values():
            try:
                # hyperlink to ENA record
                # TODO: This only works for SILVA where accessions are derived
                # from ENA accessions. May not work with other databases e.g.
                # Greengenes
                rec["tophit"] = (
                    "["
                    + rec["tname"]
                    + "](https://www.ebi.ac.uk/ena/browser/view/"
                    + rec["tname"].split(".")[0]
                    + ")"
                )
                rec["tophit taxonomy"] = ";".join(self._acc2tax[rec["tname"]])
                rec["tophit species"] = self._acc2tax[rec["tname"]][-1]
                # Higher taxonomy to class level, except for chloroplast and mitochondria
                # Assumes SILVA taxonomy is in use
                if rec["tophit taxonomy"].startswith(
                    "Bacteria;Cyanobacteria;Cyanobacteriia;Chloroplast",
                ):
                    rec["higher taxonomy"] = "[Eukaryotic organelle Chloroplast]"
                elif rec["tophit taxonomy"].startswith(
                    "Bacteria;Proteobacteria;Alphaproteobacteria;Rickettsiales;Mitochondria",
                ):
                    rec["higher taxonomy"] = "[Eukaryotic organelle Mitochondria]"
                else:
                    rec["higher taxonomy"] = "; ".join(
                        self._acc2tax[rec["tname"]][0:-1],
                    )
            except KeyError as e:
                e.add_note(f"Accession {rec['tname']} not found in database?")
                raise
        self._stats.update({"cluster_tophits": out})

    @check_stage_file(
        stage="report_json",
        message="Writing report stats in JSON format",
    )
    def write_report_json(self) -> None:
        """Dump run stats file in JSON format.

        Dump the Pipeline._stats attribute to a JSON file for troubleshooting
        and comparison of different phyloblitz runs.
        """
        self._stats.update({"endtime": str(datetime.now())})
        with Path.open(self.pathto("report_json"), "w") as fh:
            json.dump(self._stats, fh, indent=4)

    def write_cluster_alns(self) -> None:
        """Dump cluster alignments for troubleshooting.

        Write alignments in Fasta format with suffix "_aln_cluster_{c}.fasta"
        where `c` is the cluster number.
        """
        try:
            for c in self._stats["cluster cons parsed"]:
                with Path.open(
                    Path(self._outdir)
                    / Path(self._prefix + "_aln_cluster_" + str(c) + ".fasta"),
                    "w",
                ) as fh:
                    for hdr, seq in self._stats["cluster cons parsed"][c].items():
                        fh.write(">" + hdr + "\n" + seq + "\n")
        except KeyError as e:
            e.add_note("Key `cluster cons parsed` not found, has spoa been run?")
            raise

    @check_stage_file(
        stage="report_dvs_hist",
        message="Generating plots for report",
    )
    def write_report_histogram(self) -> None:
        """Write histogram graphic required for report."""
        if "dv_max_applied" in self._stats["runstats"]:
            vline = (self._stats["runstats"]["dv_max applied"],)
        else:
            vline = None
        generate_histogram(
            vals=self._stats["min_dvs"],
            vline=vline,
            title="Histogram of ava min dvs",
            outfile=self.pathto("report_dvs_hist"),
            figsize=(3, 2),
        )

    @check_stage_file(
        stage="report_md",
        message="Writing report markdown",
    )
    def write_report_markdown(self) -> None:
        """Write final report in markdown format."""
        with Path.open(self.pathto("report_md"), "w") as fh:
            fh.write(
                generate_report_md(
                    self._stats,
                    self.pathto("report_dvs_hist", basename_only=True),
                    self.pathto("report_kmercount_plot", basename_only=True),
                ),
            )

    @check_stage_file(
        stage="report_html",
        message="Writing report HTML",
    )
    def write_report_html(self) -> None:
        """Write final report in HTML format."""
        with Path.open(self.pathto("report_html"), "w") as fh:
            fh.write(
                generate_report_html(
                    self._stats,
                    self.pathto("report_dvs_hist", basename_only=True),
                    self.pathto("report_kmercount_plot", basename_only=True),
                ),
            )
