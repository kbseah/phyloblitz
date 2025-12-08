#!/usr/bin/env python3

import logging
import argparse
import json
import sys
import re

import phyloblitz.pipeline as pipeline

from os import makedirs
from datetime import datetime


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
        "--dv_max_auto",
        help="Set dv_max parameter automatically at 2 * median of all-vs-all divergence value",
        default=False,
        action="store_true",
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
    parser.add_argument(
        "--keeptmp", help="Do not delete temp files", default=False, action="store_true"
    )
    parser.add_argument(
        "--log",
        help="Write logging messages to this file",
    )
    parser.add_argument(
        "--debug",
        help="Display logging DEBUG level messages to console",
        default=False,
        action="store_true",
    )

    return parser.parse_args()


def main():
    logging.basicConfig(level=logging.DEBUG)
    root_logger = logging.getLogger()
    root_logger.handlers.clear()  # avoid duplicate handlers
    logger = logging.getLogger(__name__)  # Logger for this module

    args = init_args()

    console_handler = logging.StreamHandler(sys.stderr)
    if args.debug:
        console_handler.setLevel(logging.DEBUG)
    else:
        console_handler.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "%(levelname)s : %(module)s : %(asctime)s : %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    if args.log:
        logfile_handler = logging.FileHandler(args.log)
        logfile_handler.setFormatter(formatter)
        logfile_handler.setLevel(logging.DEBUG)
        root_logger.addHandler(logfile_handler)

    logger.debug("Arguments:")
    for i in vars(args):
        logger.debug(f" {i} : {str(vars(args)[i])}")

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

    p = pipeline.Pipeline(args)
    p.run_minimap(threads=args.threads, mode="map-" + args.platform)

    if args.twopass:
        p.twopass_extract_read_intervals(minlen=args.align_minlen)
        p.run_minimap_secondmap(threads=args.threads, mode="map-" + args.platform)

    # All-vs-all mapping
    p.extract_reads_for_ava(twopass=args.twopass, align_minlen=args.align_minlen)
    p.ava_map(mode="ava-" + args.platform, threads=args.threads)
    p.paf_get_dvs()
    p.paf_abc(dv_max=args.dv_max, dv_max_auto=args.dv_max_auto)

    # Clustering
    p.mcxload()
    p.mcl_cluster()
    p.assemble_clusters(threads=args.threads, keeptmp=args.keeptmp)

    # Taxonomy summary
    p.db_taxonomy()
    p.summarize_initial_mapping_taxonomy(
        twopass=args.twopass,
        minlen=args.align_minlen,
        taxlevel=args.summary_taxlevel,
    )
    p.cluster_asm_tophits(threads=args.threads)
    p.summarize_tophit_paf()

    # Write reports
    p.write_report_json()

    if not args.noreport:
        p.write_report_histogram()
        p.write_report_markdown()
        p.write_report_html()

    logger.info("-------------- phyloblitz run complete --------------")

    if args.log:
        logfile_handler.close()
