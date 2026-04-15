"""phyloblitz run pipeline."""

import logging
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from subprocess import PIPE, Popen
from tempfile import NamedTemporaryFile, TemporaryDirectory

import matplotlib.pyplot as plt
import numpy as np
import oxli
import pymarkovclustering as pymcl
import pysam

from phyloblitz.__about__ import __version__
from phyloblitz.report import (
    generate_histogram,
    generate_report_html,
    generate_report_md,
)
from phyloblitz.utils import (
    Pipeline,
    check_stage_file,
    filter_paf_overhang,
    lists_common_prefix,
    run_md5,
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
            ),
        )


class Run(Pipeline):
    """phyloblitz run pipeline object."""

    def __init__(self, args: dict) -> None:
        """Construct Run object.

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
        self._stats = {
            "version": __version__,
            "args": args,
            "runstats": {},
            "starttime": str(datetime.now()),
            "db_md5": run_md5(args["db"]),
        }
        logger.debug(
            "Database file %s has MD5 checksum %s",
            self._ref,
            self._stats["db_md5"],
        )

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

    @check_stage_file(
        stage="initial_map",
        message="Initial mapping of reads to identify target intervals",
    )
    def run_minimap(
        self,
        threads: int = 12,
        mode: str = "map-ont",
        sample: int | None = None,
        keeptmp: bool = False,
    ) -> int:
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
            `lr:hq`, `map-hifi`
        :param sample: Number of reads to sample; if None then use all reads;
            also use this number as the random seed
        :param keeptmp: Keep temporary file?
        :returns: return code for samtools
        :rtype: tuple
        """
        # Don't use mappy because it doesn't support all-vs-all yet
        infile = self._reads
        if sample is not None:
            logger.info("Taking sample of %d reads", sample)
            # TODO: use pyfastx Python API instead of CLI
            temp_infile = NamedTemporaryFile(
                suffix=".fastq",
                mode="w",
                delete=(not keeptmp),
                delete_on_close=False,
            )
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
        extracted, not to the Fastq file but Run._stats["segments"], so
        they do not influence the marker clustering and assembly.

        :param align_minlen: Minimum alignment length in bp
        :param no_supp: Ignore supplementary alignments if True
        :param flanking: Length of flanking sequence to extract, in bp
        """
        sam_file = self.pathto("initial_map")
        counter = 0
        self._stats["reads2segment"] = {}
        self._stats["segments"] = {}
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
                    f"@{name!s}\n{rec['seq']!s}\n+\n{rec['quals']!s}\n",
                )
                self._stats["segments"][name] = rec
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
            cmd = [
                "minimap2",
                *presets[mode],
                "-t",
                str(threads),
                self.pathto("mapped_segments"),
                self.pathto("mapped_segments"),
            ]
            logger.debug("minimap command: %s", " ".join([str(i) for i in cmd]))
            proc = Popen(cmd, stdout=paf_fh, stderr=PIPE)
            for l in proc.stderr:
                logger.debug("  minimap log: %s", l.decode().rstrip())
            return proc.wait()

    @check_stage_file(
        stage="ava_filter",
        message="Filtering overlapping incompatible overhangs from all-vs-all mappings",
    )
    def paf_file_filter_overhangs(self, max_overhang_frac: float = 0.05) -> None:
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
        the cutoffs for clustering with the mcl method. Updates Run._stats
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
                    msg = f"Problem in PAF file {self.pathto('ava_filter')!s}: more than one dv tag in entry"
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

    @check_stage_file(stage="mcl_cluster", message="Clustering with mcl")
    def pymcl_cluster(
        self,
        dv_max: float = 0.03,
        dv_max_auto: bool = False,
        inflation: float = 2,
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

    @check_stage_file(stage="isonclust3_cluster", message="Clustering with isonclust3")
    def isonclust3_cluster(self) -> int:
        """Cluster marker read segments with isonclust3."""
        # isonclust3 takes output folder path as argument, automatically
        # creates `clustering` subfolder, so go two levels up
        outfolder = self.pathto("isonclust3_cluster").parent.parent
        reads = self.pathto("mapped_segments")
        return super().isonclust3_cluster(outfolder, reads)

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
        """Extract cluster sequences and assemble with spoa."""
        if cluster_tool == "mcl":
            cluster_out = self.pathto("mcl_cluster")
        elif cluster_tool == "isonclust3":
            cluster_out = self.pathto("isonclust3_cluster")
        return super().assemble_clusters(
            cluster_out=cluster_out,
            reads=self.pathto("mapped_segments"),
            cluster_asm=self.pathto("cluster_asm"),
            cluster_tool=cluster_tool,
            threads=threads,
            rseed=rseed,
            keeptmp=keeptmp,
            min_clust_size=min_clust_size,
            max_clust_size=max_clust_size,
        )

    def cluster_flanking_isonclust3(self) -> None:
        """Cluster flanking sequences with isonclust3 and report number of groups.

        rRNA marker genes are more conserved than most genes and intergenic
        sequences. Different strains of the same species may therefore have
        similar or even identical marker sequence, and end up represented in
        the same cluster. As a measure of the potential diversity in each
        cluster, sequence flanking the conserved markers is also extracted and
        clustered, and the resulting number of clusters is reported. Stored in
        Run._stats["cluster flanking numclust"]
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
                        self._stats["segments"][seqid]["pre"]
                        + "NNNNNNNNNN"
                        + self._stats["segments"][seqid]["post"]
                    )
                    quals = (
                        self._stats["segments"][seqid]["pre_quals"]
                        + "@@@@@@@@@@"
                        + self._stats["segments"][seqid]["post_quals"]
                    )
                    fq.write(f"@{seqid!s}\n{seq!s}\n+\n{quals!s}\n")
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
                        "isonclust3 command: %s",
                        " ".join([str(i) for i in cmd]),
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
        Stored in Run._stats["flanking_kmer_histo"]

        :param k: K-mer length, in bp
        :param minlen: Only count k-mers for flanking sequence at least this length
        """
        logger.info("Counting k-mers from flanking sequences per cluster")
        kmer_tables = {}
        for cluster in self._stats["cluster2seq"]:
            kmer_tables[cluster] = oxli.KmerCountTable(k)
            for seqid in self._stats["cluster2seq"][cluster]:
                if seqid in self._stats["segments"]:
                    if len(self._stats["segments"][seqid]["pre"]) > minlen:
                        kmer_tables[cluster].consume(
                            self._stats["segments"][seqid]["pre"],
                        )
                    if len(self._stats["segments"][seqid]["post"]) > minlen:
                        kmer_tables[cluster].consume(
                            self._stats["segments"][seqid]["post"],
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
            squeeze=False,  # account for len(histos)==1
            figsize=(4, fig_height),
            layout="constrained",
        )
        for idx, (cid, histo) in enumerate(histos.items()):
            axs[idx, 0].bar(
                x=[i[0] for i in histo],
                height=[np.log1p(i[1]) for i in histo],
                width=1,
            )
            axs[idx, 0].set_title(cid, y=1.0, pad=-15, loc="right")
            axs[idx, 0].set_xlim(-2, max_kmer_cov + 2)
        fig.savefig(self.pathto("report_kmercount_plot"))

    def _per_read_consensus_taxonomy(
        self,
        sam_file: str | Path,
        minlen: int = 1200,
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
        self,
        minlen: int,
        taxlevel: int = 4,
    ) -> None:
        """Summarize taxonomy at a specified taxonomy level from initial mapping.

        Updates Run._stats["initial_taxonomy"] with a Counter of the
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
        """Map assembled cluster consensus sequences to reference database with minimap2."""
        return super().cluster_asm_tophits(
            self.pathto("cluster_tophits"),
            self.pathto("cluster_asm"),
            threads=threads,
        )

    def summarize_tophit_paf(self) -> None:
        """Summarize top hits of assembled seqs mapped to SILVA database by minimap2."""
        return super().summarize_tophit_paf(tophits=self.pathto("cluster_tophits"))

    @check_stage_file(stage="report_json", message="Writing report stats as JSON")
    def write_report_json(self) -> None:
        """Dump run stats file in JSON format."""
        return super().write_report_json(out=self.pathto("report_json"))

    def write_cluster_alns(self) -> None:
        """Dump cluster alignments for troubleshooting.

        Write alignments in Fasta format with suffix "_aln_cluster_{c}.fasta"
        where `c` is the cluster number.
        """
        try:
            for c in self._stats["cluster cons parsed"]:
                with Path.open(
                    Path(self._outdir)
                    / Path(f"{self._prefix!s}_aln_cluster_{c!s}.fasta"),
                    "w",
                ) as fh:
                    for hdr, seq in self._stats["cluster cons parsed"][c].items():
                        fh.write(f">{hdr!s}\n{seq!s}\n")
        except KeyError as e:
            e.add_note("Key `cluster cons parsed` not found, has spoa been run?")
            raise

    @check_stage_file(stage="report_dvs_hist", message="Generating plots for report")
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

    @check_stage_file(stage="report_md", message="Writing report as Markdown")
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

    @check_stage_file(stage="report_html", message="Writing report as HTML")
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
