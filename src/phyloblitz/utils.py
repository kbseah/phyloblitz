#!/usr/bin/env python3

import logging
import re
import os.path

from collections import defaultdict

logger = logging.getLogger(__name__)

CIGAROPS = {
    "M": "match",
    "I": "insertion",  # to reference
    "D": "deletion",  # to reference
    "N": "skipped",  # region to reference
    "S": "soft clipping",
    "H": "hard clipping",
    "P": "padding",
    "=": "seq match",
    "X": "seq mismatch",
}


def parse_cigar_ops(cigar):
    """Summarize operations in a CIGAR string"""
    summary = defaultdict(int)
    cigar = re.findall(r"(\d+)(\D)", cigar)
    for num, op in cigar:
        summary[op] += int(num)
    return summary


def lists_common_prefix(lol):
    """Get common prefix in a list of lists

    :param lol: list of lists of strings
    :returns: list of the common prefix
    :rtype: list
    """
    out = []
    for j in range(min([len(l) for l in lol])):
        s = set([i[j] for i in lol])
        if len(s) == 1:
            out.append(s.pop())
        else:
            break
    return out


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


def check_run_file(args, stage):
    """Check if intermediate output file has been created

    :param args: Command line arguments parsed by ArgumentParser.parse_args
    :param stage: Name of run stage, must be a key of OUTFILE_SUFFIX
    :returns: True if file already exists at expected path
    """
    return os.path.isfile(pathto(args, stage))


def pathto(args, stage, basename_only=False):
    """Combine output directory prefix and filenames to intermediate file path

    :param args: Command line arguments parsed by ArgumentParser.parse_args
    :param stage: Name of run stage, must be a key of OUTFILE_SUFFIX
    :param basename_only: Only report the base filename if True
    :returns: Expected path to intermediate output file
    """
    try:
        if basename_only:
            return os.path.basename(args.prefix + OUTFILE_SUFFIX[stage])
        else:
            return os.path.join(args.outdir, args.prefix + OUTFILE_SUFFIX[stage])
    except KeyError:
        raise Exception(f"Unknown intermediate file {stage}")
