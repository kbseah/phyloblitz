#!/usr/bin/env python3

import argparse
import pysam
import os.path
import sys
import re
import pyfastx
import logging
import json

from subprocess import Popen, PIPE
from multiprocessing import Pool
from os import makedirs
from tempfile import NamedTemporaryFile
from collections import defaultdict


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
    "report_json": "_report.json",
}

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="%(levelname)s: %(asctime)s : %(message)s",
    level=logging.DEBUG,
    datefmt="%Y-%m-%d %H:%M:%S",
)


def check_run_file(args, stage):
    """Check if intermediate output file has been created

    :param args: Command line arguments parsed by ArgumentParser.parse_args
    :param stage: Name of run stage, must be a key of OUTFILE_SUFFIX
    :returns: True if file already exists at expected path
    """
    return os.path.isfile(pathto(args, stage))


def pathto(args, stage):
    """Combine output directory prefix and filenames to intermediate file path

    :param args: Command line arguments parsed by ArgumentParser.parse_args
    :param stage: Name of run stage, must be a key of OUTFILE_SUFFIX
    :returns: Expected path to intermediate output file
    """
    try:
        return os.path.join(args.outdir, args.prefix + OUTFILE_SUFFIX[stage])
    except KeyError:
        raise Exception(f"Unknown intermediate file {stage}")


def run_minimap(ref, reads, sam_file, threads=12, mode="map-ont"):
    """Map reads to reference database with minimap2

    Filter with samtools -F 4 to discard reads that are not mapped.

    :param ref: Path to reference sequence file
    :param reads: Path to read file
    :type reads: list
    :param sam_file: Path to write SAM file output from mapping
    :param threads: Number of threads for minimap2 to use
    :param mode: Mapping preset for minimap2, either `map-ont` or `map-pb`
    """
    # Don't use mappy because it doesn't support all-vs-all yet
    logger.info("Mapping reads to reference database with minimap2")
    with open(sam_file, "w") as sam_fh:
        cmd1 = ["minimap2", "-ax", mode, f"-t {str(threads)}", ref, reads]
        logger.debug("minimap command: " + " ".join(cmd1))
        proc1 = Popen(cmd1, stdout=PIPE)  # TODO pipe stderr to log file
        proc2 = Popen(
            ["samtools", "view", "-h", "-F 4"], stdin=proc1.stdout, stdout=sam_fh
        )
        proc1.stdout.close()
        return proc2.wait()


def alignment_first_pass(sam_file, minlen=1200):
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
        proc = Popen(cmd, stdout=paf_fh)
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
            ), "Problem in PAF file {paf_file}: more than one dv tag in entry"
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
    proc = Popen(cmd)
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
    proc = Popen(cmd)
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
    proc = Popen(cmd, stdout=PIPE, text=True)
    return proc.communicate()[0]


def db_taxonomy(silvadb):
    """Get taxonomy string from SILVA headers in database Fasta file"""
    acc2tax = {}
    with open(silvadb, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                spl = line.rstrip().split(" ")
                acc = spl[0]
                taxstring = " ".join(spl[1:]).split(";")
                acc2tax[acc] = taxstring
    return acc2tax


def main():
    parser = argparse.ArgumentParser(
        prog="phyloblitz",
        description="SSU rRNA profile from ONT or PacBio long reads",
    )
    parser.add_argument(
        "-d", "--db", help="Path to preprocessed SILVA database fasta file"
    )
    parser.add_argument("-r", "--reads", help="Fastq or Fastq.gz read file to screen")
    parser.add_argument(
        "--platform",
        help="Sequencing platform used, either `pb` or `ont`",
        default="ont",
    )
    parser.add_argument("-p", "--prefix", help="Output filename prefix", default="pbz")
    parser.add_argument("-o", "--outdir", help="Output folder path", default="pbz_test")
    parser.add_argument(
        "-t", "--threads", help="Number of parallel threads", default=12, type=int
    )
    parser.add_argument(
        "--align_minlen",
        help="Minimum length of aligned segment",
        default=1200,
        type=int,
    )
    parser.add_argument(
        "--resume",
        help="Resume partially completed run based on expected filenames",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--keeptmp", help="Do not delete temp files", default=False, action="store_true"
    )
    parser.add_argument(
        "--log",
        help="Write logging messages to this file",
    )
    args = parser.parse_args()

    stats = {}
    stats["args"] = vars(args)

    logger.debug("Arguments:")
    for i in vars(args):
        logger.debug(f" {i} : {str(vars(args)[i])}")

    if not args.db or not args.reads:
        parser.print_help(sys.stderr)
        sys.exit(1)

    assert args.platform in ["ont", "pb"], "--platform must be either ont or pb"

    if args.log:
        logfh = logging.FileHandler(args.log)
        formatter = logging.Formatter(
            "%(levelname)s: %(asctime)s : %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )
        logfh.setFormatter(formatter)
        logger.addHandler(logfh)

    logger.info("Starting phyloblitz run ... ")

    logger.info(f"Creating output folder {args.outdir}")
    try:
        makedirs(args.outdir, exist_ok=False)
    except FileExistsError:
        if not args.resume:
            logger.error(
                f"Output folder {args.outdir} already exists, and option --resume not used"
            )
            sys.exit(1)
        else:
            logger.error(f"Output folder {args.outdir} already exists, resuming run")

    logger.info("Reading taxonomy from SILVA database file")
    acc2tax = db_taxonomy(args.db)
    logger.debug(f" Accessions read: {str(len(acc2tax))}")

    if not check_run_file(args, "initial_map"):
        logger.info("Initial mapping of reads to identify target intervals")
        map_ret = run_minimap(
            args.db,
            args.reads,
            pathto(args, "initial_map"),
            mode="map-" + args.platform,
            threads=args.threads,
        )

    if not check_run_file(args, "intervals_fastq"):
        logger.info("Retrieve aligned intervals on reads")
        merged_intervals = alignment_first_pass(pathto(args, "initial_map"))
        stats["merged_intervals"] = merged_intervals
        extract_fastq_read_intervals(
            merged_intervals, args.reads, pathto(args, "intervals_fastq")
        )

    if not check_run_file(args, "second_map"):
        logger.info("Second mapping of extracted intervals for taxonomic summary")
        map_ret = run_minimap(
            args.db,
            pathto(args, "intervals_fastq"),
            pathto(args, "second_map"),
            mode="map-" + args.platform,
            threads=args.threads,
        )

    if not check_run_file(args, "mapped_segments"):
        counter = 0
        with open(pathto(args, "mapped_segments"), "w") as fq_fh:
            for name, seq, quals in sam_seq_generator(
                pathto(args, "second_map"), minlen=args.align_minlen
            ):
                counter += 1
                fq_fh.write("@" + name + "\n")
                fq_fh.write(seq + "\n")
                fq_fh.write("+" + "\n")
                fq_fh.write(quals + "\n")
        logger.info(f"Read segments extracted by second mapping: {str(counter)}")

    if not check_run_file(args, "ava_map"):
        ava_ret = ava_map(
            pathto(args, "mapped_segments"),
            pathto(args, "ava_map"),
            mode="ava-" + args.platform,
            threads=args.threads,
        )
        stats["dvs"] = paf_get_dvs(pathto(args, "ava_map"))

    if not check_run_file(args, "ava_abc"):
        abc_ret = paf_abc(pathto(args, "ava_map"), pathto(args, "ava_abc"), dv_max=0.03)

    if not check_run_file(args, "ava_mci") and not check_run_file(args, "ava_seqtab"):
        mcx_ret = mcxload(
            pathto(args, "ava_abc"), pathto(args, "ava_mci"), pathto(args, "ava_seqtab")
        )

    if not check_run_file(args, "mcl_cluster"):
        mcl_ret = mcl_cluster(
            pathto(args, "ava_mci"),
            pathto(args, "ava_seqtab"),
            pathto(args, "mcl_cluster"),
        )

    if not check_run_file(args, "cluster_asm"):
        fastq_handles, cluster2seq = cluster_seqs(
            pathto(args, "mcl_cluster"),
            pathto(args, "mapped_segments"),
        )
        stats["cluster2seq"] = cluster2seq

        with Pool(args.threads) as pool:
            cluster_cons = pool.map(
                spoa_assemble, [i.name for i in fastq_handles.values()]
            )

        for i in fastq_handles.values():  # close NamedTemporaryFile handles
            i.close()

        with open(pathto(args, "cluster_asm"), "w") as fh:
            for cluster, seq in zip(fastq_handles.keys(), cluster_cons):
                fh.write(
                    re.sub(r"^>Consensus", f">cluster_{str(cluster)} Consensus", seq)
                )
            logger.info(f"Assembled sequences written to {pathto(args, 'cluster_asm')}")

    if not check_run_file(args, "report_json"):
        with open(pathto(args, "report_json"), "w") as fh:
            logger.info("Writing report stats to " + pathto(args, "report_json"))
            json.dump(stats, fh, indent=4)

    logger.info("-------------- phyloblitz run complete --------------")

    if args.log:
        logfh.close()
