#!/usr/bin/env python3

import logging
import os.path

logger = logging.getLogger(__name__)

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
