#!/usr/bin/env python3

import logging
import json
import sys
import re

import phyloblitz.pipeline as pipeline
import phyloblitz.report as report

from phyloblitz.utils import pathto, check_run_file
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

    stats = {}  # Collect data and metadata for reporting stats
    stats["args"] = vars(args)
    stats["runstats"] = {}  # initialize dict to capture run statistics

    if args.log:
        logfile_handler = logging.FileHandler(args.log)
        formatter = logging.Formatter(
            "%(levelname)s : %(module)s : %(asctime)s : %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        logfile_handler.setFormatter(formatter)
        root_logger.addHandler(logfile_handler)

    logger.debug("Arguments:")
    for i in vars(args):
        logger.debug(f" {i} : {str(vars(args)[i])}")

    logger.info("Starting phyloblitz run ... ")

    # TODO put this in pipeline initialization
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
    p.run_minimap(threads=args.threads)

    if args.twopass:  # TODO refactor this
        if not check_run_file(args, "intervals_fastq"):
            logger.info("Retrieve aligned intervals on reads")
            merged_intervals = pipeline.get_firstpass_intervals(
                pathto(args, "initial_map")
            )
            stats["merged_intervals"] = merged_intervals
            pipeline.extract_fastq_read_intervals(
                merged_intervals, args.reads, pathto(args, "intervals_fastq")
            )

        if not check_run_file(args, "second_map"):
            logger.info("Second mapping of extracted intervals for taxonomic summary")
            map_ret, nreads = pipeline.run_minimap(
                args.db,
                args.dbindex,
                pathto(args, "intervals_fastq"),
                pathto(args, "second_map"),
                mode="map-" + args.platform,
                threads=args.threads,
            )
        sam_for_read_extraction = pathto(args, "second_map")

    else:
        sam_for_read_extraction = pathto(args, "initial_map")

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
        if not check_run_file(args, "report_dvs_hist"):
            logger.info("Generating plots for report")
            report.generate_report_plots(stats, args)
        if not check_run_file(args, "report_md"):
            with open(pathto(args, "report_md"), "w") as fh:
                logger.info("Writing report markdown to " + pathto(args, "report_md"))
                fh.write(report.generate_report_md(stats, args))

        if not check_run_file(args, "report_html"):
            with open(pathto(args, "report_html"), "w") as fh:
                logger.info("Writing report HTML to " + pathto(args, "report_html"))
                fh.write(report.generate_report_html(stats, args))

    logger.info("-------------- phyloblitz run complete --------------")

    if args.log:
        logfile_handler.close()
