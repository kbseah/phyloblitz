#!/usr/bin/env python3

import argparse
import pysam
import os.path
import sys
import re

from subprocess import Popen, PIPE


OUTFILE_SUFFIX = {
    "initial_map" : "_minimap.sam",
    "mapped_segments" : "_mapped.fastq",
    "ava_map" : "_ava.paf",
    "ava_abc" : "_ava.abc",
    "ava_mci" : "_ava.mci",
    "ava_seqtab" : "_ava_seq.tab",
    "mcl_cluster" : "_mcl.out",
}


def check_run_file(args, stage):
    return os.path.isfile(
        os.path.join(args.outdir, args.prefix + OUTFILE_SUFFIX[stage])
    )


def pathto(args, stage):
    try:
        return os.path.join(args.outdir, args.prefix + OUTFILE_SUFFIX[stage])
    except KeyError:
        raise Exception(f"Unknown intermediate file {stage}")


def run_minimap(ref, reads, sam_file, threads=12, mode="map-ont"):
    # TODO use mappy
    with open(sam_file, "w") as sam_fh:
        cmd1 = ["minimap2", "-ax", mode, f"-t {str(threads)}", ref] + reads
        print(" ".join(cmd1))
        proc1 = Popen(cmd1, stdout=PIPE) # TODO pipe stderr to log file
        proc2 = Popen(["samtools", "view", "-h", "-F 4"], stdin=proc1.stdout, stdout=sam_fh)
        proc1.stdout.close()
        return proc2.wait()


def sam_seq_generator(sam_file, minlen=1200):
    sam = pysam.AlignmentFile(sam_file, "r")
    for i in sam.fetch():
        if i.query_alignment_length >= minlen and i.query_alignment_sequence:
            # Primary mappings only to ensure only one mapping per input read
            if not i.is_secondary and not i.is_supplementary:
                name = ":".join([str(i) for i in [i.query_name, i.query_alignment_start, i.query_alignment_end]])
                seq = i.query_alignment_sequence
                quals = i.query_qualities_str[i.query_alignment_start:i.query_alignment_end]
                yield (name, seq, quals)


def ava_map(reads, paf_file, mode="ava-ont", threads=12):
    with open(paf_file, "w") as paf_fh:
        cmd = ["minimap2", "-x", mode, "-t", str(threads), reads, reads]
        print(" ".join(cmd))
        proc = Popen(cmd, stdout=paf_fh)
        return proc.wait()


def paf_abc(paf_file, abc_file, dv_max=0.03):
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
    proc = Popen(cmd)
    return proc.wait()


def mcl_cluster(mci_file, tab_file, mcl_out, inflation=2):
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
    proc = Popen(cmd)
    return proc.wait()


def main():
    parser = argparse.ArgumentParser(
        prog="phyloblitz",
        description="SSU rRNA profile from ONT or PacBio long reads",
    )
    parser.add_argument("-d", "--db", help="Path to preprocessed SILVA database fasta file")
    parser.add_argument("-r", "--reads", help="Comma-separated list of read files to screen")
    parser.add_argument("-p", "--prefix", help="Output filename prefix", default="pbz")
    parser.add_argument("-o", "--outdir", help="Output folder path", default="pbz_test")
    parser.add_argument("-t", "--threads", help="Number of parallel threads", default=12)
    parser.add_argument("--align_minlen", help="Minimum length of aligned segment", default=1200)
    args = parser.parse_args()

    if not args.db or not args.reads:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if not check_run_file(args, "initial_map"):
        map_ret = run_minimap(args.db, args.reads.split(","), pathto(args, "initial_map"), threads=args.threads, )

    if not check_run_file(args, "mapped_segments"):
        counter = 0
        with open(pathto(args, "mapped_segments"), "w") as fq_fh:
            for (name, seq, quals) in sam_seq_generator(pathto(args, "initial_map")):
                counter += 1
                fq_fh.write("@" + name + "\n")
                fq_fh.write(seq + "\n")
                fq_fh.write("+" + "\n")
                fq_fh.write(quals + "\n")
        print(f"Reads extracted: {str(counter)}") # TODO proper logging

    if not check_run_file(args, "ava_map"):
        ava_ret = ava_map(pathto(args, "mapped_segments"), pathto(args, "ava_map"), mode="ava-ont", threads=args.threads)

    if not check_run_file(args, "ava_abc"):
        abc_ret = paf_abc(pathto(args, "ava_map"), pathto(args, "ava_abc"), dv_max=0.03)

    if not check_run_file(args, "ava_mci") and not check_run_file(args, "ava_seqtab"):
        mcx_ret = mcxload(pathto(args, "ava_abc"), pathto(args, "ava_mci"), pathto(args, "ava_seqtab"))

    if not check_run_file(args, "mcl_cluster"):
        mcl_ret = mcl_cluster(pathto(args, "ava_mci"), pathto(args, "ava_seqtab"), pathto(args, "mcl_cluster"))
