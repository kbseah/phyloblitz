#!/usr/bin/env python3

import logging
import json
import sys
import re

import phyloblitz.pipeline as pipeline
import rich_click as click

from os import makedirs
from datetime import datetime
from phyloblitz.__about__ import __version__
from phyloblitz.utils import check_dependencies


click.rich_click.OPTION_GROUPS = {
    "phyloblitz": [
        {
            "name": "Input",
            "options": ["db", "dbindex", "reads", "num_reads", "platform"],
        },
        {
            "name": "Output",
            "options": [
                "prefix",
                "outdir",
                "report",
                "keeptmp",
                "log",
                "write_cluster_alns",
            ],
        },
        {
            "name": "Run parameters",
            "options": [
                "threads",
                "cluster_tool",
                "align_minlen",
                "summary_taxlevel",
                "min_clust_size",
                "resume",
                "debug",
            ],
        },
        {
            "name": "ava + mcl options",
            "help": "Only used if --cluster_tool mcl is specified",
            "options": ["dv_max", "dv_max_auto", "inflation"],
        },
        {
            "name": "Experimental",
            "options": ["twopass", "flanking"],
        },
    ]
}

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(
    context_settings=CONTEXT_SETTINGS,
    help="Rapid SSU rRNA marker gene screening of long read metagenomes.",
)
@click.option(
    "--db",
    "-d",
    help="Path to preprocessed SILVA database fasta file",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--dbindex",
    help="Path to minimap2 index of database file (optional)",
    type=click.Path(exists=True),
    default=None,
)
@click.option(
    "--reads",
    "-r",
    help="Fastq or Fastq.gz read file to screen",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--num_reads",
    help="Use subsample of N reads; default is to use all reads",
    default=None,
    type=int,
)
@click.option(
    "--platform",
    help="Sequencing platform used, argument passed to minimap2 -x option",
    type=click.Choice(["map-ont", "map-pb", "lr:hq", "map-hifi"]),
    default="lr:hq",
    show_default=True,
)
@click.option(
    "--prefix", "-p", help="Output filename prefix", default="pbz", show_default=True
)
@click.option(
    "--outdir",
    "-o",
    help="Output folder path",
    default="pbz_test",
    show_default=True,
    type=click.Path(),
)
@click.option(
    "--report/--noreport",
    help="Generate report file",
    default=True,
    is_flag=True,
    show_default=True,
)
@click.option("--keeptmp", help="Do not delete temp files", default=False, is_flag=True)
@click.option("--log", help="Write logging messages to this file", type=click.Path())
@click.option(
    "--write_cluster_alns",
    help="Write cluster alignments to files",
    default=False,
    is_flag=True,
    show_default=True,
)
@click.option(
    "--threads",
    "-t",
    help="Number of parallel threads",
    default=12,
    type=int,
    show_default=True,
)
@click.option(
    "--cluster_tool",
    help="Tool(s) to use for sequence clustering",
    type=click.Choice(["mcl", "isonclust3"]),
    default="isonclust3",
    show_default=True,
)
@click.option(
    "--align_minlen",
    help="Minimum length of aligned segment",
    default=1200,
    type=int,
    show_default=True,
)
@click.option(
    "--summary_taxlevel",
    help="Depth of taxonomy string for summary in report",
    default=4,
    type=int,
    show_default=True,
)
@click.option(
    "--min_clust_size",
    help="Minimum cluster size to assemble a consensus sequence",
    default=5,
    type=int,
    show_default=True,
)
@click.option(
    "--resume",
    help="Resume partially completed run based on expected filenames",
    default=False,
    is_flag=True,
    show_default=True,
)
@click.option(
    "--debug",
    help="Display logging DEBUG level messages to console",
    default=False,
    is_flag=True,
)
@click.option(
    "--dv_max",
    help="Maximum pairwise sequence divergence in minimap2 all-vs-all mapping to retain for clustering",
    default=0.03,
    type=float,
    show_default=True,
)
@click.option(
    "--dv_max_auto",
    help="Set dv_max parameter automatically at max(0.001, 2 * median of all-vs-all divergence value)",
    default=False,
    is_flag=True,
    show_default=True,
)
@click.option(
    "--inflation", help="Inflation parameter for MCL", default=2, show_default=True
)
@click.option(
    "--twopass",
    help="[EXPERIMENTAL] Extract read segments and map again to reference",
    default=False,
    is_flag=True,
)
@click.option(
    "--flanking",
    help="[EXPERIMENTAL] Sequence flanking the mapped hits on query reads to extract",
    default=0,
    type=int,
)
@click.version_option(version=__version__)
def main(
    db,
    dbindex,
    reads,
    num_reads,
    platform,
    prefix,
    outdir,
    report,
    keeptmp,
    log,
    write_cluster_alns,
    threads,
    cluster_tool,
    align_minlen,
    summary_taxlevel,
    min_clust_size,
    resume,
    debug,
    dv_max,
    dv_max_auto,
    inflation,
    twopass,
    flanking,
):
    logging.basicConfig(level=logging.DEBUG)
    root_logger = logging.getLogger()
    root_logger.handlers.clear()  # avoid duplicate handlers
    logger = logging.getLogger(__name__)  # Logger for this module

    # args = init_args()
    args = dict(
        zip(
            [
                "db",
                "dbindex",
                "reads",
                "num_reads",
                "platform",
                "prefix",
                "outdir",
                "report",
                "keeptmp",
                "log",
                "write_cluster_alns",
                "threads",
                "cluster_tool",
                "align_minlen",
                "summary_taxlevel",
                "min_clust_size",
                "resume",
                "debug",
                "dv_max",
                "dv_max_auto",
                "inflation",
                "twopass",
                "flanking",
            ],
            [
                db,
                dbindex,
                reads,
                num_reads,
                platform,
                prefix,
                outdir,
                report,
                keeptmp,
                log,
                write_cluster_alns,
                threads,
                cluster_tool,
                align_minlen,
                summary_taxlevel,
                min_clust_size,
                resume,
                debug,
                dv_max,
                dv_max_auto,
                inflation,
                twopass,
                flanking,
            ],
        )
    )

    console_handler = logging.StreamHandler(sys.stderr)
    if debug:
        console_handler.setLevel(logging.DEBUG)
    else:
        console_handler.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "%(levelname)s : %(module)s : %(asctime)s : %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    if log:
        logfile_handler = logging.FileHandler(log)
        logfile_handler.setFormatter(formatter)
        logfile_handler.setLevel(logging.DEBUG)
        root_logger.addHandler(logfile_handler)

    logger.debug("Arguments:")
    for i in args:
        logger.debug(f" {i} : {str(args[i])}")

    logger.debug("Dependencies:")
    deps = check_dependencies()
    for d in deps:
        logger.debug(f"  {d} : {deps[d]}")

    logger.info("Starting phyloblitz run ... ")

    logger.info(f"Creating output folder {outdir}")
    try:
        makedirs(outdir, exist_ok=False)
    except FileExistsError:
        if not resume:
            logger.error(
                f"Output folder {outdir} already exists, and option --resume not used"
            )
            sys.exit(1)
        else:
            logger.error(f"Output folder {outdir} already exists, resuming run")

    p = pipeline.Pipeline(args)
    p.run_minimap(threads=threads, mode=platform, sample=num_reads)

    if twopass:
        logger.info("[EXPERIMENTAL] Applying two-pass mode")
        p.twopass_extract_read_intervals(minlen=align_minlen)
        p.run_minimap_secondmap(threads=threads, mode="map-" + platform)

    # Extract reads for clustering
    p.extract_reads_for_ava(
        twopass=twopass, align_minlen=align_minlen, flanking=flanking
    )

    # All-vs-all mapping
    p.ava_map(mode=platform, threads=threads)
    p.paf_file_filter_overhangs(max_overhang_frac=0.05)
    p.paf_get_dvs()

    # Clustering
    if cluster_tool == "mcl":
        p.pymcl_cluster(dv_max=dv_max, dv_max_auto=dv_max_auto, inflation=inflation)
    elif cluster_tool == "isonclust3":
        p.isonclust3_cluster()
    p.assemble_clusters(
        cluster_tool=cluster_tool,
        threads=threads,
        keeptmp=keeptmp,
        min_clust_size=min_clust_size,
    )

    # Cluster flanking sequences
    if flanking > 500:
        p.cluster_flanking_isonclust3()

    # Taxonomy summary
    p.db_taxonomy()
    p.summarize_initial_mapping_taxonomy(
        twopass=twopass,
        minlen=align_minlen,
        taxlevel=summary_taxlevel,
    )
    p.cluster_asm_tophits(threads=threads)
    p.summarize_tophit_paf()

    # Write reports
    p.write_report_json()
    if write_cluster_alns:
        p.write_cluster_alns()

    if report:
        p.write_report_histogram()
        p.write_report_markdown()
        p.write_report_html()

    logger.info("-------------- phyloblitz run complete --------------")

    if log:
        logfile_handler.close()
