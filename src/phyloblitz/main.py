import logging
import sys

import rich_click as click
from rich.logging import RichHandler

from phyloblitz import downloads, pipeline
from phyloblitz.__about__ import __version__
from phyloblitz.utils import check_dependencies, check_outdir

click.rich_click.OPTION_GROUPS = {
    "phyloblitz download": [],
    "phyloblitz run": [
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
                "flanking",
                "no_supplementary",
                "resume",
                "debug",
            ],
        },
        {
            "name": "ava + mcl options",
            "help": "Only used if --cluster_tool mcl is specified",
            "options": ["dv_max", "dv_max_auto", "inflation"],
        },
    ],
}

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
MIN_FLANKING = 500


@click.group(
    context_settings=CONTEXT_SETTINGS,
    help="Rapid rRNA marker gene screening of long read metagenomes.",
)
@click.version_option(version=__version__)
def main() -> None:
    """Top level command."""


@main.command(
    context_settings=CONTEXT_SETTINGS,
    help="Download reference databases.",
    no_args_is_help=True,
)
@click.option(
    "--list_versions",
    "-l",
    help="List available database versions and exit",
    default=False,
    is_flag=True,
)
@click.option(
    "--which_db",
    help="Which database to download",
    type=click.Choice(["SSU", "LSU"]),
    default="SSU",
    show_default=True,
)
@click.option(
    "--db_version",
    help="Version of database to download, use --list to see available versions; if not specified, will download latest version",
    default="latest",
    show_default=True,
)
@click.option(
    "--outdir",
    "-o",
    help="Output folder path",
    default="pbz_db",
    show_default=True,
    type=click.Path(),
)
@click.option(
    "--dryrun",
    "-n",
    help="Only print download URL without downloading the file",
    default=False,
    is_flag=True,
)
@click.option(
    "--overwrite",
    help="Overwrite existing file if it exists",
    default=False,
    is_flag=True,
)
@click.option(
    "--debug",
    help="Display logging DEBUG level messages to console",
    default=False,
    is_flag=True,
)
@click.option("--log", help="Write logging messages to this file", type=click.Path())
def download(
    list_versions, which_db, db_version, outdir, dryrun, overwrite, debug, log
) -> None:
    """Command line interface to download reference databases from Zenodo."""
    logging.basicConfig(level=logging.DEBUG)
    root_logger = logging.getLogger()
    root_logger.handlers.clear()  # avoid duplicate handlers
    logger = logging.getLogger(__name__)  # Logger for this module

    loglevel = logging.DEBUG if debug else logging.INFO
    root_logger.addHandler(RichHandler(level=loglevel))

    formatter = logging.Formatter(
        "%(levelname)s : %(module)s : %(asctime)s : %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    if log:
        logfile_handler = logging.FileHandler(log)
        logfile_handler.setFormatter(formatter)
        logfile_handler.setLevel(logging.DEBUG)
        root_logger.addHandler(logfile_handler)

    versions, latest = downloads.list_versions()
    logger.debug("Latest version is %s", latest)
    if list_versions:
        for v, meta in versions.items():
            report = f"Version: {v}, DOI: {meta['doi']}, Created: {meta['created']}"
            if v == str(latest):
                print("* " + report)
            else:
                print("  " + report)
            for marker, filedata in meta["files"].items():
                print(
                    f"    Marker: {marker}, Filename: {filedata['filename']}, Size: {filedata['size']} bytes, Checksum: {filedata['checksum']}"
                )
        sys.exit(0)
    if db_version == "latest":
        logger.info("Using latest version %s of database", str(latest))
        db_version = str(latest)

    logger.info("Creating output folder %s", outdir)
    try:
        check_outdir(outdir, resume=True)
    except Exception as e:
        logger.error(str(e))
        sys.exit(1)

    try:
        filepath = downloads.get_file(
            versions, which_db, db_version, outdir, dryrun, overwrite
        )
    except Exception as e:
        logger.error(str(e))
        sys.exit(1)

    if dryrun:
        logger.info("Dry run complete, no file downloaded")
    else:
        logger.info("Database file downloaded to %s", filepath)

    try:
        checksum_ok = downloads.check_md5sum_file(
            versions, which_db, db_version, filepath
        )
    except Exception as e:
        logger.error(str(e))
        sys.exit(1)
    if checksum_ok:
        logger.info("MD5 checksum OK")
        logger.info("------------ phyloblitz download complete ------------")
    else:
        logger.error("MD5 checksum does not match expected value!")
        sys.exit(1)

    if log:
        logfile_handler.close()


@main.command(
    context_settings=CONTEXT_SETTINGS,
    help="Run phyloblitz pipeline.",
    no_args_is_help=True,
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
    "--prefix",
    "-p",
    help="Output filename prefix",
    default="pbz",
    show_default=True,
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
    "--flanking",
    help="Sequence flanking the mapped hits on query reads to extract",
    default=1000,
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
    "--no_supplementary",
    help="Ignore supplementary alignments, only parse one marker sequence per read",
    default=False,
    is_flag=True,
)
def run(
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
    flanking,
    no_supplementary,
):
    """Command line interface to run phyloblitz pipeline."""
    logging.basicConfig(level=logging.DEBUG)
    root_logger = logging.getLogger()
    root_logger.handlers.clear()  # avoid duplicate handlers
    logger = logging.getLogger(__name__)  # Logger for this module

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
                "flanking",
                "no_supplementary",
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
                flanking,
                no_supplementary,
            ],
            strict=False,
        ),
    )

    loglevel = logging.DEBUG if debug else logging.INFO
    root_logger.addHandler(RichHandler(level=loglevel))

    formatter = logging.Formatter(
        "%(levelname)s : %(module)s : %(asctime)s : %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    if log:
        logfile_handler = logging.FileHandler(log)
        logfile_handler.setFormatter(formatter)
        logfile_handler.setLevel(logging.DEBUG)
        root_logger.addHandler(logfile_handler)

    logger.debug("Arguments:")
    for arg, val in args.items():
        logger.debug(" %s : %s", str(arg), str(val))

    logger.debug("Dependencies:")
    deps = check_dependencies()
    for dep, ver in deps.items():
        logger.debug("  %s : %s", dep, ver)

    logger.info("Starting phyloblitz run ... ")

    logger.info("Creating output folder %s", outdir)
    try:
        check_outdir(outdir, resume=resume)
    # catch any exceptions and log them, then exit with error code 1
    except Exception as e:
        logger.error(str(e))
        sys.exit(1)

    p = pipeline.Pipeline(args)
    p.run_minimap(threads=threads, mode=platform, sample=num_reads)

    # Extract reads for clustering
    p.extract_reads_for_ava(
        align_minlen=align_minlen,
        no_supplementary=no_supplementary,
        flanking=flanking,
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
    if flanking < MIN_FLANKING:
        logger.info("Flanking sequence length must be >=500 bp; setting to 500 bp")
        flanking = MIN_FLANKING
    p.cluster_flanking_isonclust3()
    p.cluster_flanking_kmercount(k=11, minlen=500)
    p.cluster_flanking_kmercount_plot(min_clust_size=min_clust_size)

    # Taxonomy summary
    p.db_taxonomy()
    p.summarize_initial_mapping_taxonomy(
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
