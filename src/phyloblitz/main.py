#!/usr/bin/env python3

import logging
import json
import sys
import re

import phyloblitz.pipeline as pipeline

from multiprocessing import Pool
from os import makedirs
from datetime import datetime


def main():
    logging.basicConfig(level=logging.DEBUG)
    root_logger = logging.getLogger()
    root_logger.handlers.clear()  # avoid duplicate handlers
    logger = logging.getLogger(__name__)  # Logger for this module

    args = pipeline.init_args()

    console_handler = logging.StreamHandler(sys.stderr)
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
    p.paf_abc(dv_max=args.dv_max)

    # Clustering
    p.mcxload()
    p.mcl_cluster()
    p.assemble_clusters(threads=args.threads)

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
