#!/usr/bin/env python

import argparse
import pysam
import re
import pyfastx
import logging

from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile
from collections import defaultdict

logger = logging.getLogger(__name__)


def run_minimap(ref, refindex, reads, sam_file, threads=12, mode="map-ont"):
    """Map reads to reference database with minimap2

    Filter with samtools -F 4 to discard reads that are not mapped.

    :param ref: Path to reference sequence file
    :param refindex: Path to reference index file (optional)
    :param reads: Path to read file
    :type reads: list
    :param sam_file: Path to write SAM file output from mapping
    :param threads: Number of threads for minimap2 to use
    :param mode: Mapping preset for minimap2, either `map-ont` or `map-pb`
    :returns: return code for samtools and number of mapped sequences parsed
    from minimap2 log message
    :rtype: tuple
    """
    # Don't use mappy because it doesn't support all-vs-all yet
    logger.info("Mapping reads to reference database with minimap2")
    with open(sam_file, "w") as sam_fh:
        cmd1 = ["minimap2", "-ax", mode, f"-t {str(threads)}"]
        (
            cmd1.extend([refindex, reads])
            if refindex is not None
            else cmd1.extend([ref, reads])
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
            logger.debug("  minimap log: " + l.decode().rstrip())  # output is in bytes
        proc1.stdout.close()
        return proc2.wait(), nreads


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
    """Filter initial alignment to get primary mappings for all-vs-all mapping

    :param sam_file: Path to SAM file from initial mapping step
    :param minlen: Minimum query alignment length; adjust if targeting a different gene, e.g. LSU rRNA
    """
    logger.info(
        f"Filtering initial alignment to get primary mappings with length >= {str(minlen)}"
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


def ava_map(reads, paf_file, mode="ava-ont", threads=12):
    """All-vs-all mapping with minimap2 to generate clusters for assembly

    :param reads: Path to Fastq reads to be clustered
    :param paf_file: Path to write output PAF alignment
    :param mode: Mapping preset mode, either `ava-ont` or `ava-pb`
    :param threads: Number of threads for minimap2 to use
    """
    logger.info("All-vs-all mapping of mapped reads with minimap2")
    with open(paf_file, "w") as paf_fh:
        cmd = ["minimap2", "-x", mode, "-t", str(threads), reads, reads]
        logger.debug("minimap command: " + " ".join(cmd))
        proc = Popen(cmd, stdout=paf_fh, stderr=PIPE)
        for l in proc.stderr:
            logger.debug("  minimap log: " + l.decode().rstrip())
        return proc.wait()


def paf_get_dvs(paf_file):
    dvs = defaultdict(list)
    with open(paf_file, "rt") as fh:
        for line in fh:
            spl = line.rstrip().split("\t")
            query = spl[0]
            dv = re.findall(r"dv:f:([\d\.]+)", line)
            assert (
                len(dv) == 1
            ), f"Problem in PAF file {paf_file}: more than one dv tag in entry"
            dvs[query].append(float(dv[0]))
    return dvs


def paf_abc(paf_file, abc_file, dv_max=0.03):
    """Convert PAF alignment to ABC format for MCL clustering

    :param paf_file: Path to PAF alignment from previous all-vs-all mapping step
    :param abc_file: Path to write ABC file
    :param dv_max: Maximum sequence divergence (dv:f tag in PAF file) to accept
    """
    logger.info("Converting PAF ava alignment to ABC format for clustering")
    fh_abc = open(abc_file, "w")
    with open(paf_file, "rt") as fh_paf:
        for line in fh_paf:
            spl = line.rstrip().split("\t")
            query = spl[0]
            target = spl[5]
            dv = re.findall(r"dv:f:([\d\.]+)", line)
            if len(dv) == 1 and float(dv[0]) < dv_max:
                score = 1 - float(dv[0])
                fh_abc.write("\t".join([str(i) for i in [query, target, score]]) + "\n")
    return fh_abc.close()


def mcxload(abc_file, mci_file, tab_file):
    """Convert ABC file to cluster input files for MCL

    :param abc_file: Path to ABC file from previous conversion step
    :param mci_file: Path to write MCI file
    :param tab_file: Path to write Tab file
    """
    logger.info("Preparing clustering files with mcxload")
    cmd = [
        "mcxload",
        "-abc",
        abc_file,
        "--stream-mirror",
        "-o",
        mci_file,
        "-write-tab",
        tab_file,
    ]
    logger.debug("mcxload command: " + " ".join(cmd))
    proc = Popen(cmd, stderr=PIPE)
    for l in proc.stderr:
        logger.debug("  mcxload log: " + l.decode().rstrip())
    return proc.wait()


def mcl_cluster(mci_file, tab_file, mcl_out, inflation=2):
    """Clustering with MCL to get sequence clusters for assembly

    :param mci_file: Path to MCI file produced by mcxload from previous step
    :param tab_file: Path to tab file produced by mcxload from previous step
    :param mcl_out: Path to write MCL output
    :param inflation: Inflation parameter passed to mcl -I option
    """
    logger.info(f"Clustering with mcl and inflation parameter {str(inflation)}")
    cmd = [
        "mcl",
        mci_file,
        "-I",
        str(inflation),
        "-use-tab",
        tab_file,
        "-o",
        mcl_out,
    ]
    logger.debug("mcl command: " + " ".join(cmd))
    proc = Popen(cmd, stderr=PIPE)
    for l in proc.stderr:
        logger.debug("  mcl log: " + l.decode().rstrip())
    return proc.wait()


def cluster_seqs(mcl_out, reads):
    """Cluster sequences with MCL and write to temporary Fastq files

    :param mcl_out: Path to write MCL output
    :param reads: Path to reads from sam_seq_generator step
    """
    counter = 0
    cluster_files = []
    fastq_handles = {}
    seq2cluster = {}
    with open(mcl_out, "r") as fh:
        for line in fh:
            seqs = line.rstrip().split("\t")
            logger.info(f"Cluster {str(counter)} comprises {str(len(seqs))} sequences")
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


def spoa_assemble(fastq):
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


def cluster_asm_tophits(ref, refindex, fasta_file, paf_file, threads=12):
    """Map assembled sequences to reference DB and get top hits

    :param ref: Path to reference sequence file
    :param refindex: Path to reference index file (optional)
    :param paf_file: Path to write PAF file output from mapping
    :param threads: Number of threads for minimap2 to use
    """
    logger.info("Mapping assembled sequences to reference database with minimap2")
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
        paf_file,
    ]
    (
        cmd.extend([refindex, fasta_file])
        if refindex is not None
        else cmd.extend([ref, fasta_file])
    )
    logger.debug("minimap command: " + " ".join(cmd))
    proc = Popen(cmd, stderr=PIPE)
    for l in proc.stderr:
        logger.debug("  minimap log: " + l.decode().rstrip())
    return proc.wait()


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
