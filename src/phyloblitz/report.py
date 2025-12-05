#!/usr/bin/env python3

import re
import logging
import pysam

import numpy as np
import matplotlib.pyplot as plt

from phyloblitz.__about__ import __version__
import phyloblitz.utils as utils
from mistune.renderers.markdown import MarkdownRenderer
from mistune import create_markdown, html
from collections import defaultdict

logger = logging.getLogger(__name__)


def dict2markdowntable(
    d,
    keys=None,
    col1="Parameter",
    col2="Value",
    order_by_value=False,
    order_descending=True,
):
    """Convert dict to a markdown table

    Key will be cast as row names

    :param d: Dict to convert to table
    :param keys: Keys to include in table; if None, all keys will be used
    :param col1: Column header for first column
    :param col2: Column header for second column
    :returns: Text table in markdown format
    :rtype: str
    """
    if order_by_value:
        ordered_keys = (
            sorted(d, key=lambda x: -d[x])
            if order_descending
            else sorted(d, key=lambda x: d[x])
        )
    else:
        ordered_keys = list(d.keys())
    if keys is not None:
        keys = [i for i in ordered_keys if i in keys]
    else:
        keys = ordered_keys
    out = []
    out.append(f"| {col1} | {col2} |")
    out.append(f"| :----- | :----- |")
    for k in keys:
        try:
            out.append(f"| {str(k)} | {str(d[k])} |")
        except KeyError:
            raise Exception(f"Key {k} not defined in dict")
    return "\n".join(out)


def dod2markdowntable(d, keys, col1="Name"):
    """Convert dict of dicts to a markdown table

    The first key will be cast as row names, second key as column names.

    :param d: Dict of dicts
    :param keys: Keys to use as columns, must be keys of the second-level dicts
    :param col1: Column header to use for the first column
    :returns: Text table in markdown format
    :rtype: str
    """
    out = []
    out.append(f"| {col1} | " + " | ".join(keys) + " |")
    out.append("| :----: | " + " | ".join([":----:"] * len(keys)) + " |")
    for c in d:
        out.append(
            f"| {str(c)} | "
            + " | ".join([str(d[c][k]) if k in d[c] else "-" for k in keys])
            + " |"
        )
    return "\n".join(out)


def generate_report_plots(stats, args):
    """Generate plots for report file

    :param stats: `stats` dict produced in phyloblitz.main.main
    :param args: Command line arguments from ArgumentParser.parse_args()
    """

    if not utils.check_run_file(args, "report_dvs_hist"):
        dvs = [min(stats["dvs"][r]) for r in stats["dvs"]]  # TODO
        stats["runstats"]["ava min dvs median"] = "{:.4f}".format(np.median(dvs))
        fig, axs = plt.subplots(1, figsize=(3, 2))
        axs.hist(dvs, bins="auto", density=True)
        axs.axvline(args.dv_max, color="red")
        axs.set_title("Histogram of ava min dvs")
        fig.tight_layout()
        fig.savefig(utils.pathto(args, "report_dvs_hist"))


def generate_report_md(stats, args):
    """Generate markdown report from stats collected during phyloblitz run

    :param stats: `stats` dict produced in phyloblitz.main.main
    :param args: Command line arguments from ArgumentParser.parse_args()
    :returns: Report in markdown format
    :rtype: str
    """
    cluster_table_fields = [
        "numseq",
        "tophit",
        "higher taxonomy",
        "tophit species",
        "align %id",
        "alnlen",
        "query %aln",
        "target %aln",
        "seq match",
        "seq mismatch",
        "insertion",
        "deletion",
    ]

    raw = f"""# phyloblitz run report

* Run started: {str(stats["starttime"])}
* Run ended: {str(stats["endtime"])}
* phyloblitz version: {__version__}

phyloblitz [homepage](https://github.com/kbseah/phyloblitz)


## Input parameters

phyloblitz was called with the following parameters:

{dict2markdowntable(stats["args"])}


## Run statistics

{dict2markdowntable(stats["runstats"])}

![]({utils.pathto(args, 'report_dvs_hist', basename_only=True)})


## Taxonomy summary from initial mapping

{dict2markdowntable(stats['initial_taxonomy'], order_by_value=True)}


## Assembled sequence clusters

Summary of assembled sequence clusters and their top hits in the reference
database. Cluster sequences with few underlying reads and many
mismatches/indels to the reference hits are lower quality because of
insufficient coverage, and should not be used for phylogenetics or probe
design.

{dod2markdowntable(per_cluster_summarize(stats), cluster_table_fields, col1='Cluster ID')}

"""

    format_md = create_markdown(renderer=MarkdownRenderer())
    return format_md(raw)


def generate_report_html(stats, args):
    """Generate HTML report from stats collected during phyloblitz run"""

    return html(generate_report_md(stats, args))


def per_cluster_summarize(stats):
    out = defaultdict(lambda: defaultdict(dict))
    for c in stats["cluster2seq"]:
        out["cluster_" + str(c)]["numseq"] = len(stats["cluster2seq"][c])
    for c in stats["cluster_tophits"]:
        out[c].update(stats["cluster_tophits"][c])
    return out
