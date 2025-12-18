#!/usr/bin/env python

import pysam
import re
import pyfastx
import logging
import os.path
import json

import numpy as np

from phyloblitz.report import (
    generate_histogram,
    generate_report_md,
    generate_report_html,
)
from phyloblitz.utils import (
    CIGAROPS,
    lists_common_prefix,
    parse_cigar_ops,
    filter_paf_overhang,
)
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile
from collections import defaultdict, Counter
from functools import wraps
from multiprocessing import Pool
from datetime import datetime

logger = logging.getLogger(__name__)


def get_firstpass_intervals(sam_file, minlen=1200):
    """Filter initial alignment to get primary mappings for all-vs-all mapping

    :param sam_file: Path to SAM file from initial mapping step
    :param minlen: Minimum query alignment length; adjust if targeting a different gene, e.g. LSU rRNA
    :returns: Lists of intervals with hits per read, keyed by read name
    :rtype: dict
    """
    logger.info(f"Filtering initial alignment to get aligned regions to extract")
    regions = defaultdict(list)
    sam = pysam.AlignmentFile(sam_file, "r")
    for i in sam.fetch():
        if i.query_alignment_length >= minlen:
            # include secondary and supplementary alignments
            # if sequence is left-hardclipped, correct the coordinates
            # because pysam's query_alignment_start does not count hardclips!
            qstart, qend = i.query_alignment_start, i.query_alignment_end
            if i.cigartuples[0][0] == 5:
                qstart += i.cigartuples[0][1]
                qend += i.cigartuples[0][1]
            regions[i.query_name].append((qstart, qend))
    logger.debug(f"Reads processed {str(len(regions))}")
    # merge intervals for each read
    merged_intervals = {i: merge_intervals(regions[i]) for i in regions}
    return merged_intervals


def merge_intervals(intervals: list):
    """Merge overlapping numerical intervals

    :param intervals: List of tuples of (start, end) coordinates
    :returns: List of tuples of intervals with overlaps merged
    :rtype: list
    """
    intervals.sort(key=lambda i: i[0])  # sort in place
    merged = [intervals[0]]
    for curr in intervals[1:]:
        if curr[0] < merged[-1][1]:
            merged[-1] = (merged[-1][0], max(curr[1], merged[-1][1]))
        else:
            merged.append(curr)
    return merged


def sam_seq_generator(sam_file, minlen=1200, flanking=0):
    """Filter SAM alignment to get primary mappings for all-vs-all mapping

    :param sam_file: Path to SAM file
    :param minlen: Minimum query alignment length; adjust if targeting a
        different gene, e.g. LSU rRNA
    :param flanking: [EXPERIMENTAL] Additional flanking sequence to extract
    """
    logger.info(
        f"Filtering alignment for primary mappings with length >= {str(minlen)}"
    )
    if flanking > 0:
        logger.info(
            f"Including flanking {str(flanking)} bp sequence when extracting reads"
        )
    sam = pysam.AlignmentFile(sam_file, "r")
    for i in sam.fetch():
        if i.query_alignment_length >= minlen and i.query_alignment_sequence:
            # Primary mappings only to ensure only one mapping per input read
            if not i.is_secondary and not i.is_supplementary:
                mystart = max(0, i.query_alignment_start - flanking)
                myend = min(i.query_length, i.query_alignment_end + flanking)
                name = ":".join([str(i) for i in [i.query_name, mystart, myend]])
                seq = i.query_sequence[mystart:myend]
                quals = i.query_qualities_str[mystart:myend]
                yield (name, seq, quals)


def spoa_assemble_fasta(label_fastq):
    """Run SPOA assembly on a Fastq input file

    :param label_fastq: Tuple of input label (str) and path to Fastq file with
        reads to assemble
    :returns: Alignment of consensus and input sequences in Fasta format
        (stdout from spoa -r 2)
    :rtype: str
    """
    label, fastq = label_fastq
    cmd = [
        "spoa",
        "-r",
        "2",
        fastq,
    ]
    logger.debug("spoa command: " + " ".join(cmd))
    proc = Popen(cmd, stdout=PIPE, text=True)
    # TODO directing stderr to PIPE and logger prevents mp.Pool from closing?
    # ignoring stderr from spoa for now
    out = proc.communicate()[0]
    logger.info(f"Assembly complete for cluster {str(label)}")
    return label, out


def parse_spoa_r2(fasta):
    """Parse SPOA gapped alignment + consensus in Fasta format (-r 2 output)

    :param fasta: Assembly from SPOA as list of strings
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


def count_spoa_aln_vars(seqs):
    """Count mismatches and gaps vs consensus for each sequence in a cluster

    :param seqs: Dict of sequences parsed by parse_spoa_r2; the consensus
        sequence must have key "Consensus"
    :returns: Count of base matches and mismatches, gaps relative to query and
        consensus, leading and trailing query gaps, for each sequence in the
        dict relative to consensus
    :rtype: dict
    """
    var = defaultdict(lambda: defaultdict(int))
    for hdr in seqs:
        if hdr != "Consensus":
            # get coordinates without trailing and leading gaps
            span = re.search(r"^-*([^-].*[^-])-*$", seqs[hdr]).span(1)
            var[hdr]["query_lead_gap"] = span[0]
            var[hdr]["query_trail_gap"] = len(seqs[hdr]) - span[1]
            for i in range(span[0], span[1]):
                if seqs[hdr][i] == seqs["Consensus"][i] and seqs[hdr][i] != "-":
                    # Matching base but not common gap
                    var[hdr]["match"] += 1
                elif seqs[hdr][i] == "-" and seqs["Consensus"][i] != "-":
                    var[hdr]["query_gap"] += 1
                elif seqs[hdr][i] != "-" and seqs["Consensus"][i] == "-":
                    var[hdr]["cons_gap"] += 1
                elif (
                    seqs[hdr][i] != "-"
                    and seqs["Consensus"][i] != "-"
                    and seqs[hdr][i] != seqs["Consensus"][i]
                ):
                    var[hdr]["mismatch"] += 1
    return var


def count_spoa_aln_persite_vars(seqs):  # TODO WIP
    """Count mismatches/gaps per alignment position for clustered sequences vs consensus

    If mismatches/gaps between sequences and the consensus are solely due to
    technical sequencing error, rather than erroneous clustering of divergent
    underlying reads, we expect fraction of variants per site to be roughly the
    sequencing error. If for example two sequences are misclustered, then we
    should see a secondary peak of ~50% variant coverage.
    """
    var = {}
    for i in range(len(seqs["Consensus"])):
        column = [seqs[hdr][i] for hdr in seqs if hdr != "Consensus"]
        values, counts = np.unique(column, return_counts=True)
        counts_norm = counts / counts.sum()
        var[i] = -(counts_norm * np.log(counts_norm) / np.log(2)).sum()
    return var


class Pipeline(object):
    """phyloblitz pipeline run object"""

    def __init__(self, args):
        """Constructor for Pipeline"""
        self._ref = args.db  # Path to reference database Fasta file
        self._refindex = (
            args.dbindex
        )  # Path to reference database Minimap2 index (optional)
        self._reads = (
            args.reads
        )  # Path to reads to be processed in Fastq or Fastq.gz format
        self._outdir = args.outdir  # Path to output folder
        self._prefix = args.prefix  # Filename prefix for output files
        self._platform = args.platform  # Sequencing mode, ont or pb
        self._resume = args.resume
        self._stats = {}
        self._stats["args"] = vars(args)
        self._stats["runstats"] = {}
        self._stats["starttime"] = str(datetime.now())

    OUTFILE_SUFFIX = {
        "initial_map": "_minimap_initial.sam",
        "intervals_fastq": "_intervals.fastq",
        "second_map": "_minimap_second.sam",
        "mapped_segments": "_mapped.fastq",
        "ava_map": "_ava.paf",
        "ava_filter": "_ava_filter.paf",
        "ava_abc": "_ava.abc",
        "ava_mci": "_ava.mci",
        "ava_seqtab": "_ava_seq.tab",
        "mcl_cluster": "_mcl.out",
        "isonclust3_cluster": "_isonclust3_out/clustering/final_clusters.tsv",
        "cluster_asm": "_final.fasta",
        "cluster_tophits": "_final_tophits.paf",
        "report_json": "_report.json",
        "report_md": "_report.md",
        "report_html": "_report.html",
        "report_dvs_hist": "_report_dvs_hist.png",
    }

    def check_run_file(self, stage):
        """Check if intermediate output file has been created

        :param stage: Name of run stage, must be a key of self.OUTFILE_SUFFIX
        :returns: True if file already exists at expected path
        """
        return os.path.isfile(self.pathto(stage))

    def pathto(self, stage, basename_only=False):
        """Combine output directory prefix and filenames to intermediate file path

        :param stage: Name of run stage, must be a key of self.OUTFILE_SUFFIX
        :param basename_only: Only report the base filename if True
        :returns: Expected path to intermediate output file
        """
        try:
            if basename_only:
                return os.path.basename(self._prefix + self.OUTFILE_SUFFIX[stage])
            else:
                return os.path.join(
                    self._outdir, self._prefix + self.OUTFILE_SUFFIX[stage]
                )
        except KeyError:
            raise Exception(f"Unknown intermediate file {stage}")

    def check_stage_file(stage, message):
        # TODO specify more than one output file; check required input files
        def check_stage_decorator(func):
            @wraps(func)
            def wrapped_function(self, *args, **kwargs):
                if self.check_run_file(stage):
                    if self._resume:
                        logger.info(
                            f"Stage {stage} file output already present, skipping"
                        )
                        return
                    else:
                        logger.error(
                            f"Stage {stage} file output already present and option --resume not used, exiting"
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
    def run_minimap(self, threads=12, mode="map-ont", sample=None):
        """Map reads to reference database with minimap2

        Only report reads that are mapped.

        :param threads: Number of threads for minimap2 to use
        :param mode: Mapping preset for minimap2, either `map-ont` or `map-pb`
        :param sample: Number of reads to sample; if None then use all reads;
            also use this number as the random seed
        :returns: return code for samtools
        :rtype: tuple
        """
        # Don't use mappy because it doesn't support all-vs-all yet
        infile = self._reads
        if sample is not None:
            logger.info(f"Taking sample of {str(sample)} reads")
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
            logger.debug(f"Sampled reads in temporary file {temp_infile.name}")

        logger.info("Mapping reads to reference database with minimap2")
        with open(self.pathto("initial_map"), "w") as sam_fh:
            cmd1 = ["minimap2", "-ax", mode, "--sam-hit-only", f"-t {str(threads)}"]
            (
                cmd1.extend([self._refindex, infile])
                if self._refindex is not None
                else cmd1.extend([self._ref, infile])
            )
            logger.debug("minimap command: " + " ".join(cmd1))
            proc1 = Popen(cmd1, stdout=sam_fh, stderr=PIPE)
            nreads = 0
            for l in proc1.stderr:
                nreads_s = re.search(r"mapped (\d+) sequences", l.decode())
                if nreads_s:
                    nreads += int(nreads_s.group(1))
                logger.debug(
                    "  minimap log: " + l.decode().rstrip()
                )  # output is in bytes
            self._stats["runstats"].update({"total input reads": nreads})
            return proc1.wait()

    @check_stage_file(
        stage="second_map",
        message="[EXPERIMENTAL] Twopass mode: Second mapping of extracted intervals for taxonomic summary",
    )
    def run_minimap_secondmap(self, threads=12, mode="map-ont"):
        """Map reads to reference database with minimap2

        Only report reads that are mapped.

        :param threads: Number of threads for minimap2 to use
        :param mode: Mapping preset for minimap2, either `map-ont` or `map-pb`
        :returns: return code for samtools
        :rtype: tuple
        """
        # Don't use mappy because it doesn't support all-vs-all yet
        with open(self.pathto("second_map"), "w") as sam_fh:
            cmd1 = ["minimap2", "-ax", mode, "--sam-hit-only", f"-t {str(threads)}"]
            (
                cmd1.extend([self._refindex, self.pathto("intervals_fastq")])
                if self._refindex is not None
                else cmd1.extend([self._ref, self.pathto("intervals_fastq")])
            )
            logger.debug("minimap command: " + " ".join(cmd1))
            proc1 = Popen(cmd1, stdout=PIPE, stderr=PIPE)
            nreads = 0
            for l in proc1.stderr:
                nreads_s = re.search(r"mapped (\d+) sequences", l.decode())
                if nreads_s:
                    nreads += int(nreads_s.group(1))
                logger.debug(
                    "  minimap log: " + l.decode().rstrip()
                )  # output is in bytes
            self._stats["runstats"].update({"total input reads secondmap": nreads})
            return proc1.wait()

    @check_stage_file(
        stage="intervals_fastq",
        message="[EXPERIMENTAL] Twopass mode: Extract aligned intervals on reads",
    )
    def twopass_extract_read_intervals(self, minlen=1200):
        merged_intervals = get_firstpass_intervals(
            self.pathto("initial_map"), minlen=minlen
        )
        self._stats.update({"merged_intervals": merged_intervals})
        counter = 0
        with open(self.pathto("intervals_fastq"), "w") as fh:
            for name, seq, qual in pyfastx.Fastx(self._reads):
                if name in merged_intervals:
                    for start, end in merged_intervals[name]:
                        newname = ":".join([str(i) for i in [name, start, end]])
                        fh.write("@" + newname + "\n")
                        fh.write(seq[start:end] + "\n")
                        fh.write("+\n")
                        fh.write(qual[start:end] + "\n")
                        counter += 1
        self._stats["runstats"].update({"firstpass intervals extracted": counter})
        logger.info(
            f"[EXPERIMENTAL] Twopass mode: Read intervals extracted: {str(counter)}"
        )

    @check_stage_file(
        stage="mapped_segments",
        message="Extracting read segments for all-vs-all mapping",
    )
    def extract_reads_for_ava(self, twopass=False, align_minlen=1200, flanking=0):
        sam_file = self.pathto("second_map") if twopass else self.pathto("initial_map")
        counter = 0
        with open(self.pathto("mapped_segments"), "w") as fq_fh:
            for name, seq, quals in sam_seq_generator(
                sam_file, minlen=align_minlen, flanking=flanking
            ):
                counter += 1
                fq_fh.write("@" + name + "\n")
                fq_fh.write(seq + "\n")
                fq_fh.write("+" + "\n")
                fq_fh.write(quals + "\n")
        self._stats["runstats"].update({"mapped pass filter": counter})
        logger.info(f"Read segments extracted for all-vs-all mapping: {str(counter)}")
        return

    @check_stage_file(
        stage="ava_map",
        message="All-vs-all mapping of mapped reads with minimap2",
    )
    def ava_map(self, mode="map-ont", threads=12):
        """All-vs-all mapping with minimap2 to generate clusters for assembly

        :param mode: Mapping preset mode used for initial minimap2 mapping step
        :param threads: Number of threads for minimap2 to use
        """
        presets = {
            "map-ont": ["-x", "ava-ont"],
            "map-pb": ["-x", "ava-pb"],
            "lr:hq": ["-x", "lr:hq", "-Xw5", "-e0", "-m100", "-r2k"],
            "map-hifi": ["-x", "map-hifi", "-Xw5", "-e0", "-m100"],
        }
        with open(self.pathto("ava_map"), "w") as paf_fh:
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
            logger.debug("minimap command: " + " ".join(cmd))
            proc = Popen(cmd, stdout=paf_fh, stderr=PIPE)
            for l in proc.stderr:
                logger.debug("  minimap log: " + l.decode().rstrip())
            return proc.wait()

    @check_stage_file(
        stage="ava_filter",
        message="Filtering overlapping incompatible overhangs from all-vs-all mappings",
    )
    def paf_file_filter_overhangs(self, max_overhang_frac=0.05):
        logger.info(f"Max overhang fraction: {str(max_overhang_frac)}")
        counter_all = 0
        counter = 0
        with open(self.pathto("ava_map"), "r") as fh_in:
            with open(self.pathto("ava_filter"), "w") as fh_out:
                for line in fh_in:
                    counter_all += 1
                    out = filter_paf_overhang(line, max_overhang_frac=max_overhang_frac)
                    if out is not None:
                        counter += 1
                        fh_out.write(out)
        logger.info(
            f"Retained {str(counter)} out of {str(counter_all)} all-vs-all alignments"
        )
        return

    def paf_get_dvs(self):
        logger.info("Reading dvs data from ava mapping")
        dvs = defaultdict(list)
        with open(self.pathto("ava_filter"), "rt") as fh:
            for line in fh:
                spl = line.rstrip().split("\t")
                query = spl[0]
                dv = re.findall(r"dv:f:([\d\.]+)", line)
                assert (
                    len(dv) == 1
                ), f"Problem in PAF file {self.pathto('ava_filter')}: more than one dv tag in entry"
                dvs[query].append(float(dv[0]))
        self._stats.update({"dvs": dvs})
        # min dv for each read, to get dv distribution and median to
        # automatically set dv threshold for clustering
        min_dvs = [min(dvs[r]) for r in dvs]
        self._stats.update({"min_dvs": min_dvs})
        self._stats["runstats"].update(
            {
                "ava min dvs median": "{:.4f}".format(np.median(min_dvs)),
            }
        )
        return

    @check_stage_file(
        stage="ava_abc",
        message="Converting PAF ava alignment to ABC format for clustering",
    )
    def paf_abc(self, dv_max=0.03, dv_max_auto=False):
        """Convert PAF alignment to ABC format for MCL clustering

        :param dv_max: Maximum sequence divergence (dv:f tag in PAF file) to accept
        :param dv_max_auto: Set dv_max to max(0.001, 2 * median all-vs-all
            sequence divergence); parameter `dv_max` will be ignored). The value
            0.001 is to account for cases where median divergence is zero.
        """
        if dv_max_auto:
            dv_max = max(
                0.001, 2 * float(self._stats["runstats"]["ava min dvs median"])
            )
        fh_abc = open(self.pathto("ava_abc"), "w")
        with open(self.pathto("ava_filter"), "rt") as fh_paf:
            for line in fh_paf:
                spl = line.rstrip().split("\t")
                query = spl[0]
                target = spl[5]
                dv = re.findall(r"dv:f:([\d\.]+)", line)
                if len(dv) == 1 and float(dv[0]) < dv_max:
                    score = 1 - float(dv[0])
                    fh_abc.write(
                        "\t".join([str(i) for i in [query, target, score]]) + "\n"
                    )
        self._stats["runstats"].update({"dv_max applied": dv_max})
        return fh_abc.close()

    @check_stage_file(
        stage="ava_mci",  # and ava_seqtab
        message="Preparing clustering files with mcxload",
    )
    def mcxload(self):
        """Convert ABC file to cluster input files for MCL"""
        cmd = [
            "mcxload",
            "-abc",
            self.pathto("ava_abc"),
            "--stream-mirror",
            "-o",
            self.pathto("ava_mci"),
            "-write-tab",
            self.pathto("ava_seqtab"),
        ]
        logger.debug("mcxload command: " + " ".join(cmd))
        proc = Popen(cmd, stderr=PIPE)
        for l in proc.stderr:
            logger.debug("  mcxload log: " + l.decode().rstrip())
        return proc.wait()

    @check_stage_file(
        stage="mcl_cluster",
        message="Clustering with mcl",
    )
    def mcl_cluster(self, inflation=2):
        """Clustering with MCL to get sequence clusters for assembly

        :param inflation: Inflation parameter passed to mcl -I option
        """
        logger.info(f"Inflation parameter {str(inflation)}")
        cmd = [
            "mcl",
            self.pathto("ava_mci"),
            "-I",
            str(inflation),
            "-use-tab",
            self.pathto("ava_seqtab"),
            "-o",
            self.pathto("mcl_cluster"),
        ]
        logger.debug("mcl command: " + " ".join(cmd))
        proc = Popen(cmd, stderr=PIPE)
        for l in proc.stderr:
            logger.debug("  mcl log: " + l.decode().rstrip())
        return proc.wait()

    @check_stage_file(
        stage="isonclust3_cluster",
        message="Clustering with isonclust3",
    )
    def isonclust3_cluster(self):
        if self._platform in ["lr:hq", "map-ont"]:
            isonclust3_mode = "ont"
        elif self._platform in ["map-hifi", "map-pb"]:
            isonclust3_mode = "pacbio"
        # isonclust3 takes output folder path as argument, automatically
        # creates `clustering` subfolder
        outfolder = os.path.dirname(os.path.dirname(self.pathto("isonclust3_cluster")))
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
        logger.debug("isonclust3 command: " + " ".join(cmd))
        proc = Popen(cmd, stderr=PIPE)
        for l in proc.stderr:
            logger.debug("  isonclust3 log: " + l.decode().rstrip())
        return proc.wait()

    def _cluster_seqs_from_mcl(self, mcl_out, reads, keeptmp):
        """Extract sequences from each MCL cluster to separate Fastq files for assembly

        :param mcl_out: Path MCL output file
        :param reads: Path to reads from extract_reads_for_ava step
        :param keeptmp: Keep temporary files?
        :returns: Dict of file handles to each Fastq file keyed by cluster ID
        """
        counter = 0
        cluster_files = []
        fastq_handles = {}
        seq2cluster = {}
        with open(mcl_out, "r") as fh:
            for line in fh:
                seqs = line.rstrip().split("\t")
                logger.info(
                    f"Cluster {str(counter)} comprises {str(len(seqs))} sequences"
                )
                for seq in seqs:
                    seq2cluster[seq] = counter  # assume each seq in only one cluster
                fastq_handles[counter] = NamedTemporaryFile(
                    suffix=".fastq",
                    mode="w",
                    delete=(not keeptmp),
                    delete_on_close=False,
                )
                counter += 1
        for name, seq, qual in pyfastx.Fastx(reads):
            if name in seq2cluster:
                fastq_rec = "@" + name + "\n" + seq + "\n+\n" + qual + "\n"
                fastq_handles[seq2cluster[name]].write(fastq_rec)
        cluster2seq = defaultdict(list)
        for seq in seq2cluster:
            cluster2seq[seq2cluster[seq]].append(seq)
        for i in fastq_handles:
            fastq_handles[i].close()
        return fastq_handles, cluster2seq

    def _cluster_seqs_from_isonclust3(self, isonclust3_out, reads, keeptmp):
        """Extract sequences from each isONclust3 cluster to separate Fastq files for assembly

        The numbering of clusters by isONclust3 itself is not consistent
        between the final_clusters.tsv file and the output fastq files; best to
        extract the sequences ourselves.

        :param isonclust3_out: Path isONclust3 output file
        :param reads: Path to reads from extract_reads_for_ava step
        :param keeptmp: Keep temporary files?
        :returns: Dict of file handles to each Fastq file keyed by cluster ID
        """
        cluster_files = []
        fastq_handles = {}
        seq2cluster = {}
        with open(isonclust3_out, "r") as fh:
            for line in fh:
                (clust, seqname) = line.rstrip().split("\t")
                seq2cluster[seqname] = clust
        cluster2seq = defaultdict(list)
        for seqname in seq2cluster:
            cluster2seq[seq2cluster[seqname]].append(seqname)
        for cluster in cluster2seq:
            logger.info(
                f"Cluster {str(cluster)} comprises {str(len(cluster2seq[cluster]))} sequences"
            )
            fastq_handles[cluster] = NamedTemporaryFile(
                suffix=".fastq",
                mode="w",
                delete=(not keeptmp),
                delete_on_close=False,
            )
        for seqname, seq, qual in pyfastx.Fastx(reads):
            if seqname in seq2cluster:
                fastq_rec = "@" + seqname + "\n" + seq + "\n+\n" + qual + "\n"
                fastq_handles[seq2cluster[seqname]].write(fastq_rec)
        for i in fastq_handles:
            fastq_handles[i].close()
        return fastq_handles, cluster2seq

    @check_stage_file(
        stage="cluster_asm",
        message="Extract cluster sequences and assemble with SPOA",
    )
    def assemble_clusters(self, cluster_tool="mcl", threads=12, keeptmp=False):
        if cluster_tool == "mcl":
            fastq_handles, cluster2seq = self._cluster_seqs_from_mcl(
                self.pathto("mcl_cluster"),
                self.pathto("mapped_segments"),
                keeptmp,
            )
        elif cluster_tool == "isonclust3":
            fastq_handles, cluster2seq = self._cluster_seqs_from_isonclust3(
                self.pathto("isonclust3_cluster"),
                self.pathto("mapped_segments"),
                keeptmp,
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
            }
        )
        logger.info("Assemble consensus from clustered sequences with spoa")
        with Pool(threads) as pool:
            cluster_cons_tuples = pool.map(
                spoa_assemble_fasta,
                [(c, handle.name) for c, handle in fastq_handles.items()],
            )

        for i in fastq_handles.values():  # close NamedTemporaryFile handles
            i.close()

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
            }
        )
        with open(self.pathto("cluster_asm"), "w") as fh:
            for c in cluster_cons_parsed:
                fh.write(
                    f">cluster_{str(c)} Consensus\n"
                    + cluster_cons_parsed[c]["Consensus"].replace("-", "")
                    + "\n"
                )
            logger.info(f"Assembled sequences written to {self.pathto('cluster_asm')}")

    def db_taxonomy(self):
        """Get taxonomy string from SILVA headers in database Fasta file

        :returns: Dict of taxonomy strings keyed by SILVA accession
        :rtype: dict
        """
        logger.info("Reading taxonomy from SILVA database file")
        self._acc2tax = {}
        with open(self._ref, "r") as fh:
            for line in fh:
                if line.startswith(">"):
                    spl = line.lstrip(">").rstrip().split(" ")
                    acc = spl[0]
                    taxstring = " ".join(spl[1:]).split(";")
                    self._acc2tax[acc] = taxstring
        logger.debug(f" Accessions read: {str(len(self._acc2tax))}")

    def _per_read_consensus_taxonomy(self, sam_file, minlen=1200):
        """Summarize taxonomy from initial mapping with minimap2

        :param sam_file: Path to SAM file of initial mapping
        """
        all_taxstrings = defaultdict(list)
        sam = pysam.AlignmentFile(sam_file, "r")
        for i in sam.fetch():
            if i.query_alignment_length >= minlen:
                all_taxstrings[i.query_name].append(self._acc2tax[i.reference_name])
        common_taxstrings = {
            acc: lists_common_prefix(all_taxstrings[acc]) for acc in all_taxstrings
        }
        return common_taxstrings

    def summarize_initial_mapping_taxonomy(self, twopass, minlen, taxlevel=4):
        sam_file = self.pathto("second_map") if twopass else self.pathto("initial_map")
        logger.info("Summarizing taxonomic composition of initial mapping")
        common_taxstrings = self._per_read_consensus_taxonomy(sam_file, minlen)
        # TODO right-pad taxonomy strings if too short
        self._stats.update(
            {
                "initial_taxonomy": Counter(
                    [";".join(i[0:taxlevel]) for i in common_taxstrings.values()]
                )
            }
        )
        return

    @check_stage_file(
        stage="cluster_tophits",
        message="Mapping assembled sequences to reference database with minimap2",
    )
    def cluster_asm_tophits(self, threads=12):
        """Map assembled sequences to reference DB and get top hits

        :param threads: Number of threads for minimap2 to use
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
        logger.debug("minimap command: " + " ".join(cmd))
        proc = Popen(cmd, stderr=PIPE)
        for l in proc.stderr:
            logger.debug("  minimap log: " + l.decode().rstrip())
        return proc.wait()

    def summarize_tophit_paf(self):
        """Summarize top hits of assembled seqs mapped to SILVA database by minimap2

        :returns: Dict of summary stats for each hit, keyed by query sequence name
        :rtype: dict
        """
        out = {}

        with open(self.pathto("cluster_tophits"), "r") as fh:
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
                    )
                )
                cigar = [i for i in spl if i.startswith("cg:Z:")][0]
                # Calculate derived metrics from PAF fields
                cigar_summary = parse_cigar_ops(cigar[5:])
                hits.update({CIGAROPS[c]: cigar_summary[c] for c in cigar_summary})
                hits["align %id"] = "{:.2%}".format(
                    int(hits["alnmatch"]) / int(hits["alnlen"])
                ).rstrip(
                    "%"
                )  # remove redundant % sign for display
                hits["query %aln"] = "{:.2%}".format(
                    (int(hits["qend"]) - int(hits["qstart"])) / int(hits["qlen"])
                ).rstrip("%")
                hits["target %aln"] = "{:.2%}".format(
                    (int(hits["tend"]) - int(hits["tstart"])) / int(hits["tlen"])
                ).rstrip("%")
                out[spl[0]] = hits

        # Taxonomy of hit targets
        for c in out:
            try:
                # hyperlink to ENA record
                out[c]["tophit"] = (
                    "["
                    + out[c]["tname"]
                    + "](https://www.ebi.ac.uk/ena/browser/view/"
                    + out[c]["tname"].split(".")[0]
                    + ")"
                )
                out[c]["tophit taxonomy"] = ";".join(self._acc2tax[out[c]["tname"]])
                out[c]["tophit species"] = self._acc2tax[out[c]["tname"]][-1]
                # Higher taxonomy to class level, except for chloroplast and mitochondria
                # Assumes SILVA taxonomy is in use
                if out[c]["tophit taxonomy"].startswith(
                    "Bacteria;Cyanobacteria;Cyanobacteriia;Chloroplast"
                ):
                    out[c]["higher taxonomy"] = "[Chloroplast]"
                elif out[c]["tophit taxonomy"].startswith(
                    "Bacteria;Proteobacteria;Alphaproteobacteria;Rickettsiales;Mitochondria"
                ):
                    out[c]["higher taxonomy"] = "[Mitochondria]"
                else:
                    out[c]["higher taxonomy"] = ";".join(
                        self._acc2tax[out[c]["tname"]][0:3]
                    )
            except KeyError:
                KeyError(f"Accession {out[c]['tname']} not found in database?")
        self._stats.update({"cluster_tophits": out})

    @check_stage_file(
        stage="report_json",
        message="Writing report stats in JSON format",
    )
    def write_report_json(self):
        """Dump run stats file in JSON format for troubleshooting"""
        self._stats.update({"endtime": str(datetime.now())})
        with open(self.pathto("report_json"), "w") as fh:
            json.dump(self._stats, fh, indent=4)
        return

    def write_cluster_alns(self):
        """Dump cluster alignments for troubleshooting"""
        try:
            for c in self._stats["cluster cons parsed"]:
                with open(
                    os.path.join(
                        self._outdir, self._prefix + "_aln_cluster_" + str(c) + ".fasta"
                    ),
                    "w",
                ) as fh:
                    for hdr in self._stats["cluster cons parsed"][c]:
                        fh.write(">" + hdr + "\n")
                        fh.write(self._stats["cluster cons parsed"][c][hdr] + "\n")
        except KeyError:
            raise Exception("Key `cluster cons parsed` not found, has SPOA been run?")
        return

    @check_stage_file(
        stage="report_dvs_hist",
        message="Generating plots for report",
    )
    def write_report_histogram(self):
        """Write histogram graphic required for report"""
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
        return

    @check_stage_file(
        stage="report_md",
        message="Writing report markdown",
    )
    def write_report_markdown(self):
        """Write final report in markdown format"""
        with open(self.pathto("report_md"), "w") as fh:
            fh.write(
                generate_report_md(
                    self._stats, self.pathto("report_dvs_hist", basename_only=True)
                )
            )
        return

    @check_stage_file(
        stage="report_html",
        message="Writing report HTML",
    )
    def write_report_html(self):
        """Write final report in HTML format"""
        with open(self.pathto("report_html"), "w") as fh:
            fh.write(
                generate_report_html(
                    self._stats, self.pathto("report_dvs_hist", basename_only=True)
                )
            )
        return
