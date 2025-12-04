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
    stats["starttime"] = str(datetime.now())

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

    if not check_run_file(args, "initial_map"):
        logger.info("Initial mapping of reads to identify target intervals")
        map_ret, nreads = pipeline.run_minimap(
            args.db,
            args.dbindex,
            args.reads,
            pathto(args, "initial_map"),
            mode="map-" + args.platform,
            threads=args.threads,
        )
        stats["runstats"]["total input reads"] = nreads

    if args.twopass:
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

    if not check_run_file(args, "mapped_segments"):
        counter = 0
        with open(pathto(args, "mapped_segments"), "w") as fq_fh:
            for name, seq, quals in pipeline.sam_seq_generator(
                sam_for_read_extraction, minlen=args.align_minlen
            ):
                counter += 1
                fq_fh.write("@" + name + "\n")
                fq_fh.write(seq + "\n")
                fq_fh.write("+" + "\n")
                fq_fh.write(quals + "\n")
        stats["runstats"]["mapped pass filter"] = counter
        logger.info(f"Read segments extracted for all-vs-all mapping: {str(counter)}")

    if not check_run_file(args, "ava_map"):
        ava_ret = pipeline.ava_map(
            pathto(args, "mapped_segments"),
            pathto(args, "ava_map"),
            mode="ava-" + args.platform,
            threads=args.threads,
        )
        stats["dvs"] = pipeline.paf_get_dvs(pathto(args, "ava_map"))

    if not check_run_file(args, "ava_abc"):
        abc_ret = pipeline.paf_abc(
            pathto(args, "ava_map"), pathto(args, "ava_abc"), dv_max=args.dv_max
        )

    if not check_run_file(args, "ava_mci") and not check_run_file(args, "ava_seqtab"):
        mcx_ret = pipeline.mcxload(
            pathto(args, "ava_abc"), pathto(args, "ava_mci"), pathto(args, "ava_seqtab")
        )

    if not check_run_file(args, "mcl_cluster"):
        mcl_ret = pipeline.mcl_cluster(
            pathto(args, "ava_mci"),
            pathto(args, "ava_seqtab"),
            pathto(args, "mcl_cluster"),
        )

    if not check_run_file(args, "cluster_asm"):
        fastq_handles, cluster2seq = pipeline.cluster_seqs(
            pathto(args, "mcl_cluster"),
            pathto(args, "mapped_segments"),
        )
        stats["cluster2seq"] = cluster2seq
        stats["runstats"]["number of clusters"] = len(cluster2seq)
        stats["runstats"]["total reads in clusters"] = sum(
            [len(cluster2seq[c]) for c in cluster2seq]
        )

        with Pool(args.threads) as pool:
            cluster_cons = pool.map(
                pipeline.spoa_assemble, [i.name for i in fastq_handles.values()]
            )

        for i in fastq_handles.values():  # close NamedTemporaryFile handles
            i.close()

        with open(pathto(args, "cluster_asm"), "w") as fh:
            for cluster, seq in zip(fastq_handles.keys(), cluster_cons):
                fh.write(
                    re.sub(r"^>Consensus", f">cluster_{str(cluster)} Consensus", seq)
                )
            logger.info(f"Assembled sequences written to {pathto(args, 'cluster_asm')}")

    logger.info("Reading taxonomy from SILVA database file")
    acc2tax = report.db_taxonomy(args.db)
    logger.debug(f" Accessions read: {str(len(acc2tax))}")

    stats["initial_taxonomy"] = report.summarize_initial_mapping_taxonomy(
        pathto(args, "initial_map"),
        acc2tax,
        minlen=args.align_minlen,
        taxlevel=args.summary_taxlevel,
    )

    if not check_run_file(args, "cluster_tophits"):
        pipeline.cluster_asm_tophits(
            args.db,
            args.dbindex,
            pathto(args, "cluster_asm"),
            pathto(args, "cluster_tophits"),
            threads=args.threads,
        )
        stats["cluster_tophits"] = report.summarize_tophit_paf(
            pathto(args, "cluster_tophits"), acc2tax
        )

    stats["endtime"] = str(datetime.now())
    # comes after stats["endtime"] because it writes this information to report
    if not check_run_file(args, "report_json"):
        with open(pathto(args, "report_json"), "w") as fh:
            logger.info("Writing report stats to " + pathto(args, "report_json"))
            json.dump(stats, fh, indent=4)

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
