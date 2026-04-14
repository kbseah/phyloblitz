"""Main command line functions for phyloblitz."""

import logging
import sys

import rich_click as click
from rich.logging import RichHandler

from phyloblitz import downloads
from phyloblitz.__about__ import __version__
from phyloblitz.compare import Compare
from phyloblitz.pipeline import Run
from phyloblitz.utils import check_dependencies, check_outdir

logging.basicConfig(level=logging.DEBUG)
root_logger = logging.getLogger()
formatter = logging.Formatter(
    "%(levelname)s : %(module)s : %(asctime)s : %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

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
                "max_clust_size",
                "rseed",
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
    "phyloblitz compare": [
        {
            "name": "Input",
            "options": ["db", "dbindex", "input_table"],
        },
        {
            "name": "Output",
            "options": ["outdir", "prefix", "log"],
        },
        {
            "name": "Run parameters",
            "options": [
                "threads",
                "ignore_db_mismatch",
                "min_clust_size",
                "max_clust_size",
                "rseed",
                "debug",
            ],
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
@click.pass_context
def download(ctx, **kwargs) -> None:
    """Command line interface to download reference databases from Zenodo."""
    root_logger.handlers.clear()  # avoid duplicate handlers
    logger = logging.getLogger(__name__)  # Logger for this module

    loglevel = logging.DEBUG if ctx.params["debug"] else logging.INFO
    root_logger.addHandler(RichHandler(level=loglevel))

    if ctx.params["log"]:
        logfile_handler = logging.FileHandler(ctx.params["log"])
        logfile_handler.setFormatter(formatter)
        logfile_handler.setLevel(logging.DEBUG)
        root_logger.addHandler(logfile_handler)

    versions, latest = downloads.list_versions()
    logger.debug("Latest version is %s", latest)
    if ctx.params["list_versions"]:
        for v, meta in versions.items():
            report = f"Version: {v}, DOI: {meta['doi']}, Created: {meta['created']}"
            if v == str(latest):
                print("* " + report)
            else:
                print("  " + report)
            for marker, filedata in meta["files"].items():
                print(
                    f"    Marker: {marker}, Filename: {filedata['filename']}, Size: {filedata['size']} bytes, Checksum: {filedata['checksum']}",
                )
        sys.exit(0)
    if ctx.params["db_version"] == "latest":
        logger.info("Using latest version %s of database", str(latest))
        db_version = str(latest)

    logger.info("Creating output folder %s", ctx.params["outdir"])
    try:
        check_outdir(ctx.params["outdir"], resume=True)
    except Exception as e:
        logger.error(str(e))
        sys.exit(1)

    try:
        filepath = downloads.get_file(
            versions,
            ctx.params["which_db"],
            db_version,
            ctx.params["outdir"],
            ctx.params["dryrun"],
            ctx.params["overwrite"],
        )
    except Exception as e:
        logger.error(str(e))
        sys.exit(1)

    if ctx.params["dryrun"]:
        logger.info("Dry run complete, no file downloaded")
    else:
        logger.info("Database file downloaded to %s", filepath)

    try:
        checksum_ok = downloads.check_md5sum_file(
            versions,
            ctx.params["which_db"],
            db_version,
            filepath,
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

    if ctx.params["log"]:
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
    "--max_clust_size",
    help="Clusters above this size will be downsampled for consensus assembly",
    default=500,
    type=int,
    show_default=True,
)
@click.option(
    "--rseed",
    help="Random seed for subsampling reads and downsampling clusters",
    default=12345,
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
    "--inflation",
    help="Inflation parameter for MCL",
    default=2,
    show_default=True,
)
@click.option(
    "--no_supplementary",
    help="Ignore supplementary alignments, only parse one marker sequence per read",
    default=False,
    is_flag=True,
)
@click.pass_context
def run(ctx, **kwargs):
    """Command line interface to run phyloblitz pipeline."""
    root_logger.handlers.clear()  # avoid duplicate handlers
    logger = logging.getLogger(__name__)  # Logger for this module

    loglevel = logging.DEBUG if ctx.params["debug"] else logging.INFO
    root_logger.addHandler(RichHandler(level=loglevel))

    if ctx.params["log"]:
        logfile_handler = logging.FileHandler(ctx.params["log"])
        logfile_handler.setFormatter(formatter)
        logfile_handler.setLevel(logging.DEBUG)
        root_logger.addHandler(logfile_handler)

    logger.debug("Arguments:")
    for arg, val in ctx.params.items():
        logger.debug(" %s : %s", str(arg), str(val))

    logger.debug("Dependencies:")
    deps = check_dependencies()
    for dep, ver in deps.items():
        logger.debug("  %s : %s", dep, ver)

    logger.info("Starting phyloblitz run ... ")

    logger.info("Creating output folder %s", ctx.params["outdir"])
    try:
        check_outdir(ctx.params["outdir"], resume=ctx.params["resume"])
    # catch any exceptions and log them, then exit with error code 1
    except Exception as e:
        logger.error(str(e))
        sys.exit(1)

    p = Run(ctx.params)
    p.run_minimap(
        threads=ctx.params["threads"],
        mode=ctx.params["platform"],
        sample=ctx.params["num_reads"],
        keeptmp=ctx.params["keeptmp"],
    )

    # Extract reads for clustering
    p.extract_reads_for_ava(
        align_minlen=ctx.params["align_minlen"],
        no_supp=ctx.params["no_supplementary"],
        flanking=ctx.params["flanking"],
    )

    # All-vs-all mapping
    p.ava_map(mode=ctx.params["platform"], threads=ctx.params["threads"])
    p.paf_file_filter_overhangs(max_overhang_frac=0.05)
    p.paf_get_dvs()

    # Clustering
    if ctx.params["cluster_tool"] == "mcl":
        p.pymcl_cluster(
            dv_max=ctx.params["dv_max"],
            dv_max_auto=ctx.params["dv_max_auto"],
            inflation=ctx.params["inflation"],
        )
    elif ctx.params["cluster_tool"] == "isonclust3":
        p.isonclust3_cluster()
    p.assemble_clusters(
        cluster_tool=ctx.params["cluster_tool"],
        threads=ctx.params["threads"],
        keeptmp=ctx.params["keeptmp"],
        min_clust_size=ctx.params["min_clust_size"],
        max_clust_size=ctx.params["max_clust_size"],
        rseed=ctx.params["rseed"],
    )

    # Cluster flanking sequences
    if ctx.params["flanking"] < MIN_FLANKING:
        logger.info("Flanking sequence length must be >=500 bp; setting to 500 bp")
        ctx.params["flanking"] = MIN_FLANKING
    p.cluster_flanking_isonclust3()
    p.cluster_flanking_kmercount(k=11, minlen=500)
    p.cluster_flanking_kmercount_plot(min_clust_size=ctx.params["min_clust_size"])

    # Taxonomy summary
    p.db_taxonomy()
    p.summarize_initial_mapping_taxonomy(
        minlen=ctx.params["align_minlen"],
        taxlevel=ctx.params["summary_taxlevel"],
    )
    p.cluster_asm_tophits(threads=ctx.params["threads"])
    p.summarize_tophit_paf()

    # Write reports
    p.write_report_json()
    if ctx.params["write_cluster_alns"]:
        p.write_cluster_alns()

    if ctx.params["report"]:
        p.write_report_histogram()
        p.write_report_markdown()
        p.write_report_html()

    logger.info("-------------- phyloblitz run complete --------------")

    if ctx.params["log"]:
        logfile_handler.close()


@main.command(
    context_settings=CONTEXT_SETTINGS,
    help="Compare phyloblitz runs.",
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
    "--input_table",
    help="Path to table of sample names and report JSON file paths to compare, in TSV format with columns 'sample' and 'report'",
    type=click.Path(exists=True),
    default=None,
    required=True,
)
@click.option(
    "--outdir",
    help="Output folder path",
    default="pbz_compare",
    show_default=True,
    type=click.Path(),
)
@click.option(
    "--prefix",
    help="Output filename prefix",
    default="pbz_compare",
    show_default=True,
)
@click.option(
    "--ignore_db_mismatch",
    help="Continue even if different databases were used to generate reports",
    default=False,
    is_flag=True,
)
@click.option(
    "--threads",
    help="Number of parallel threads",
    default=12,
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
    "--max_clust_size",
    help="Clusters above this size will be downsampled for consensus assembly",
    default=500,
    type=int,
    show_default=True,
)
@click.option(
    "--rseed",
    help="Random seed for subsampling reads and downsampling clusters",
    default=12345,
    type=int,
    show_default=True,
)
@click.option(
    "--log",
    help="Path to write log file",
    type=click.Path(exists=True),
    default=None,
)
@click.option(
    "--debug",
    help="Display logging DEBUG level messages to console",
    default=False,
    is_flag=True,
)
@click.pass_context
def compare(ctx, **kwargs) -> None:
    """Compare phyloblitz runs.

    Takes JSON results files from phyloblitz runs. Runs should be produced by
    the same sequencing platform, and processed by the same phyloblitz version
    and reference database. Compares the different libraries by re-clustering
    the pooled reads and comparing which samples are represented in each
    cluster.
    """
    root_logger.handlers.clear()  # avoid duplicate handlers
    logger = logging.getLogger(__name__)  # Logger for this module

    loglevel = logging.DEBUG if ctx.params["debug"] else logging.INFO
    root_logger.addHandler(RichHandler(level=loglevel))

    if ctx.params["log"]:
        logfile_handler = logging.FileHandler(ctx.params["log"])
        logfile_handler.setFormatter(formatter)
        logfile_handler.setLevel(logging.DEBUG)
        root_logger.addHandler(logfile_handler)

    c = Compare(args=ctx.params)

    logger.info("Creating output folder %s", ctx.params["outdir"])
    try:
        check_outdir(ctx.params["outdir"], resume=True)
    except Exception as e:
        logger.error(str(e))
        sys.exit(1)

    # Check for various potential issues
    if not c.check_database_checksums():
        if ctx.params["ignore_db_mismatch"]:
            logger.info(
                "Continuing because --ignore_db_mismatch specified. I assume you know what you are doing...",
            )
        else:
            sys.exit(1)

    if not c.check_sequencing_platforms():
        sys.exit(1)

    if not c.segment_by_sample():
        logger.error("Repeated read segment names were found")
        sys.exit(1)

    c.write_segments_to_fastq()

    c.cluster_segments()

    c.assemble_clusters(
        threads=ctx.params["threads"],
        rseed=ctx.params["rseed"],
        min_clust_size=ctx.params["min_clust_size"],
        max_clust_size=ctx.params["max_clust_size"],
    )

    logger.info("Map cluster memberships to samples ...")
    c.cluster_memberships()

    c.write_report_json()
    c.write_reports()
