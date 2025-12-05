#!/usr/bin/env python

import argparse
import pysam
import re
import pyfastx
import logging
import os.path

from phyloblitz.utils import CIGAROPS, lists_common_prefix, parse_cigar_ops
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile
from collections import defaultdict, Counter
from functools import wraps
from multiprocessing import Pool

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


def extract_fastq_read_intervals(intervals, in_file, out_file):
    logger.info(
        f"Extract read intervals that align with reference database to {out_file}"
    )
    counter = 0
    with open(out_file, "w") as fh:
        for name, seq, qual in pyfastx.Fastx(in_file):
            if name in intervals:
                for start, end in intervals[name]:
                    newname = ":".join([str(i) for i in [name, start, end]])
                    fh.write("@" + newname + "\n")
                    fh.write(seq[start:end] + "\n")
                    fh.write("+\n")
                    fh.write(qual[start:end] + "\n")
                    counter += 1
    logger.info(f"Read intervals extracted: {str(counter)}")


def sam_seq_generator(sam_file, minlen=1200):
    """Filter SAM alignment to get primary mappings for all-vs-all mapping

    :param sam_file: Path to SAM file
    :param minlen: Minimum query alignment length; adjust if targeting a different gene, e.g. LSU rRNA
    """
    logger.info(
        f"Filtering alignment for primary mappings with length >= {str(minlen)}"
    )
    sam = pysam.AlignmentFile(sam_file, "r")
    for i in sam.fetch():
        if i.query_alignment_length >= minlen and i.query_alignment_sequence:
            # Primary mappings only to ensure only one mapping per input read
            if not i.is_secondary and not i.is_supplementary:
                name = ":".join(
                    [
                        str(i)
                        for i in [
                            i.query_name,
                            i.query_alignment_start,
                            i.query_alignment_end,
                        ]
                    ]
                )
                seq = i.query_alignment_sequence
                quals = i.query_qualities_str[
                    i.query_alignment_start : i.query_alignment_end
                ]
                yield (name, seq, quals)


def spoa_assemble_fasta(fastq):
    """Run SPOA assembly on a Fastq input file

    :param fastq: Path to Fastq file with reads to assemble
    :returns: Assembled sequence in Fasta format (stdout from spoa)
    :rtype: str
    """
    logger.info("Assemble consensus sequence for cluster with spoa")
    cmd = [
        "spoa",
        fastq,
    ]
    logger.debug("spoa command: " + " ".join(cmd))
    proc = Popen(cmd, stdout=PIPE, text=True, stderr=PIPE)
    for l in proc.stderr:
        logger.debug("  mcl log: " + l.decode().rstrip())
    return proc.communicate()[0]


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

    OUTFILE_SUFFIX = {
        "initial_map": "_minimap_initial.sam",
        "intervals_fastq": "_intervals.fastq",
        "second_map": "_minimap_second.sam",
        "mapped_segments": "_mapped.fastq",
        "ava_map": "_ava.paf",
        "ava_abc": "_ava.abc",
        "ava_mci": "_ava.mci",
        "ava_seqtab": "_ava_seq.tab",
        "mcl_cluster": "_mcl.out",
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

    def check_stage_file(stage, message):  # TODO more than one input file
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
    def run_minimap(self, threads=12, mode="map-ont"):
        """Map reads to reference database with minimap2

        Filter with samtools -F 4 to discard reads that are not mapped.

        :param threads: Number of threads for minimap2 to use
        :param mode: Mapping preset for minimap2, either `map-ont` or `map-pb`
        :returns: return code for samtools
        :rtype: tuple
        """
        # Don't use mappy because it doesn't support all-vs-all yet
        logger.info("Mapping reads to reference database with minimap2")
        with open(self.pathto("initial_map"), "w") as sam_fh:
            cmd1 = ["minimap2", "-ax", mode, f"-t {str(threads)}"]
            (
                cmd1.extend([self._refindex, self._reads])
                if self._refindex is not None
                else cmd1.extend([self._ref, self._reads])
            )
            logger.debug("minimap command: " + " ".join(cmd1))
            proc1 = Popen(cmd1, stdout=PIPE, stderr=PIPE)
            proc2 = Popen(
                ["samtools", "view", "-h", "-F 4"], stdin=proc1.stdout, stdout=sam_fh
            )
            nreads = 0
            for l in proc1.stderr:
                nreads_s = re.search(r"mapped (\d+) sequences", l.decode())
                if nreads_s:
                    nreads += int(nreads_s.group(1))
                logger.debug(
                    "  minimap log: " + l.decode().rstrip()
                )  # output is in bytes
            proc1.stdout.close()
            self._stats["runstats"]["total input reads"] = nreads
            return proc2.wait()

    @check_stage_file(
        stage="mapped_segments",
        message="Extracting read segments for all-vs-all mapping",
    )
    def extract_reads_for_ava(self, twopass=False, align_minlen=1200):
        sam_file = self.pathto("second_map") if twopass else self.pathto("initial_map")
        counter = 0
        with open(self.pathto("mapped_segments"), "w") as fq_fh:
            for name, seq, quals in sam_seq_generator(sam_file, minlen=align_minlen):
                counter += 1
                fq_fh.write("@" + name + "\n")
                fq_fh.write(seq + "\n")
                fq_fh.write("+" + "\n")
                fq_fh.write(quals + "\n")
        self._stats["runstats"]["mapped pass filter"] = counter
        logger.info(f"Read segments extracted for all-vs-all mapping: {str(counter)}")

    @check_stage_file(
        stage="ava_map",
        message="All-vs-all mapping of mapped reads with minimap2",
    )
    def ava_map(self, mode="ava-ont", threads=12):
        """All-vs-all mapping with minimap2 to generate clusters for assembly

        :param mode: Mapping preset mode, either `ava-ont` or `ava-pb`
        :param threads: Number of threads for minimap2 to use
        """
        with open(self.pathto("ava_map"), "w") as paf_fh:
            cmd = [
                "minimap2",
                "-x",
                mode,
                "-t",
                str(threads),
                self.pathto("mapped_segments"),
                self.pathto("mapped_segments"),
            ]
            logger.debug("minimap command: " + " ".join(cmd))
            proc = Popen(cmd, stdout=paf_fh, stderr=PIPE)
            for l in proc.stderr:
                logger.debug("  minimap log: " + l.decode().rstrip())
            return proc.wait()

    def paf_get_dvs(self):
        logger.info("Reading dvs data from ava mapping")
        dvs = defaultdict(list)
        with open(self.pathto("ava_map"), "rt") as fh:
            for line in fh:
                spl = line.rstrip().split("\t")
                query = spl[0]
                dv = re.findall(r"dv:f:([\d\.]+)", line)
                assert (
                    len(dv) == 1
                ), f"Problem in PAF file {self.pathto('ava_map')}: more than one dv tag in entry"
                dvs[query].append(float(dv[0]))
        self._stats["dvs"] = dvs

    @check_stage_file(
        stage="ava_abc",
        message="Converting PAF ava alignment to ABC format for clustering",
    )
    def paf_abc(self, dv_max=0.03):
        """Convert PAF alignment to ABC format for MCL clustering

        :param dv_max: Maximum sequence divergence (dv:f tag in PAF file) to accept
        """
        fh_abc = open(self.pathto("ava_abc"), "w")
        with open(self.pathto("ava_map"), "rt") as fh_paf:
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

    def _cluster_seqs(self, mcl_out, reads):
        """Extract sequences from each MCL cluster to separate Fastq files for assembly

        :param mcl_out: Path MCL output file
        :param reads: Path to reads from extract_reads_for_ava step
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
                    suffix=".fastq", mode="w", delete=True, delete_on_close=False
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

    @check_stage_file(
        stage="cluster_asm",
        message="Extract cluster sequences and assemble with SPOA",
    )
    def assemble_clusters(self, threads=12):
        fastq_handles, cluster2seq = self._cluster_seqs(
            self.pathto("mcl_cluster"),
            self.pathto("mapped_segments"),
        )
        self._stats["cluster2seq"] = cluster2seq
        self._stats["runstats"]["number of clusters"] = len(cluster2seq)
        self._stats["runstats"]["total reads in clusters"] = sum(
            [len(cluster2seq[c]) for c in cluster2seq]
        )

        with Pool(threads) as pool:
            cluster_cons = pool.map(
                spoa_assemble_fasta, [i.name for i in fastq_handles.values()]
            )

        for i in fastq_handles.values():  # close NamedTemporaryFile handles
            i.close()

        with open(self.pathto("cluster_asm"), "w") as fh:
            for cluster, seq in zip(fastq_handles.keys(), cluster_cons):
                fh.write(
                    re.sub(r"^>Consensus", f">cluster_{str(cluster)} Consensus", seq)
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
        self._stats["initial_taxonomy"] = Counter(
            [";".join(i[0:taxlevel]) for i in common_taxstrings.values()]
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
        self._stats["cluster_tophits"] = out


def init_args():
    parser = argparse.ArgumentParser(
        prog="phyloblitz",
        description="SSU rRNA profile from ONT or PacBio long reads",
    )
    parser.add_argument(
        "-d",
        "--db",
        help="Path to preprocessed SILVA database fasta file",
        required=True,
    )
    parser.add_argument(
        "--dbindex",
        help="Path to minimap2 index of database file (optional)",
        default=None,
    )
    parser.add_argument(
        "-r", "--reads", help="Fastq or Fastq.gz read file to screen", required=True
    )
    parser.add_argument(
        "--platform",
        help="Sequencing platform used, either `pb` or `ont`",
        choices=["ont", "pb"],
        default="ont",
    )
    parser.add_argument("-p", "--prefix", help="Output filename prefix", default="pbz")
    parser.add_argument("-o", "--outdir", help="Output folder path", default="pbz_test")
    parser.add_argument(
        "-t", "--threads", help="Number of parallel threads", default=12, type=int
    )
    parser.add_argument(
        "--twopass",
        help="[EXPERIMENTAL] Extract read segments and map again to reference",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--align_minlen",
        help="Minimum length of aligned segment",
        default=1200,
        type=int,
    )
    parser.add_argument(
        "--summary_taxlevel",
        help="Depth of taxonomy string for summary in report",
        default=4,
        type=int,
    )
    parser.add_argument(
        "--dv_max",
        help="Maximum pairwise sequence divergence in all-vs-all mapping to retain for clustering",
        default=0.03,
        type=float,
    )
    parser.add_argument(
        "--resume",
        help="Resume partially completed run based on expected filenames",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--noreport",
        help="Do not generate report file",
        default=False,
        action="store_true",
    )
    parser.add_argument(  # TODO not yet implemented
        "--keeptmp", help="Do not delete temp files", default=False, action="store_true"
    )
    parser.add_argument(
        "--log",
        help="Write logging messages to this file",
    )

    return parser.parse_args()
