#!/usr/bin/env python3

import logging
from collections import defaultdict

import matplotlib.pyplot as plt
from mistune import create_markdown, html
from mistune.renderers.markdown import MarkdownRenderer

from phyloblitz.__about__ import __version__

logger = logging.getLogger(__name__)
# mute verbose debug messages from matplotlib and PngImagePlugin
plt.set_loglevel(level="warning")
pil_logger = logging.getLogger("PIL")
pil_logger.setLevel(logging.INFO)

HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8"/>
    <meta name="viewport" content="width=device-width, initial-scale=1"/>
    <title>{{title}}</title>
    <link rel="stylesheet" href="https://cdn.simplecss.org/simple.min.css"/>
</head>
<body style="max-width:1200px">
{{markdown_report}}
</body>
</html>
"""


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
    :param order_by_value: Order rows by dict values
    :param order_descending: Descending order if true
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
    out.append("| :----- | :----- |")
    for k in keys:
        try:
            out.append(f"| {k!s} | {d[k]!s} |")
        except KeyError:
            raise Exception(f"Key {k} not defined in dict")
    return "\n".join(out)


def dod2markdowntable(
    d, keys, order_by_value=True, order_descending=True, order_by="numseq", col1="Name"
):
    """Convert dict of dicts to a markdown table

    The first key will be cast as row names, second key as column names.

    :param d: Dict of dicts
    :param keys: Keys to use as columns, must be keys of the second-level dicts
    :param col1: Column header to use for the first column
    :param order_by_value: Order rows by dict values
    :param order_descending: Descending order if true
    :param order_by: Second-order key to use for ordering (i.e. which column)
    :returns: Text table in markdown format
    :rtype: str
    """
    assert order_by in keys
    out = []
    out.append(f"| {col1} | " + " | ".join(keys) + " |")
    out.append("| :----: | " + " | ".join([":----:"] * len(keys)) + " |")
    sign = -1 if order_descending else 1
    dd = (
        sorted(d.items(), key=lambda cv: sign * cv[1][order_by])
        if order_by_value
        else d
    )
    for c, v in dd:
        out.append(
            f"| {c!s} | "
            + " | ".join([str(v[k]) if k in v else "-" for k in keys])
            + " |"
        )
    return "\n".join(out)


def generate_histogram(vals, vline, title, outfile, figsize=(3, 2)):
    """Generate plots for report file

    :param vals: Iterable of values to generate histogram for
    :param vline: x-value to draw vertical red line
    :param title: Title for the plot
    :param outfile: Path to write output file
    :param figsize: Dimensions of figure
    """

    fig, axs = plt.subplots(1, figsize=figsize)
    axs.hist(vals, bins="auto", density=True)
    if vline is not None:
        axs.axvline(vline, color="red")
    axs.set_title(title)
    fig.tight_layout()
    fig.savefig(outfile)


def generate_report_md(stats, histogram_file_path, kmercount_plot_path):
    """Generate markdown report from stats collected during phyloblitz run

    :param stats: `stats` dict produced in phyloblitz.main.main
    :param histogram_file_path: Path to histogram image file, relative to report file path
    :param kmercount_plot_path: Path to k-mer count plot file, relative to report file path
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
    if "cluster flanking numclust" in stats:
        cluster_table_fields.append("flanking seq clusters")

    raw = f"""# phyloblitz run report

* Run started: {stats["starttime"]!s}
* Run ended: {stats["endtime"]!s}
* phyloblitz version: {__version__}

phyloblitz [homepage](https://github.com/kbseah/phyloblitz)


## Input parameters

phyloblitz was called with the following parameters:

{dict2markdowntable(stats["args"])}


## Run statistics

{dict2markdowntable(stats["runstats"])}

<figure>

![]({histogram_file_path})

<figcaption>Histogram of the lowest non-self divergence score per read from
all-vs-all mapping of extracted rRNA sequences.</figcaption>
</figure>

## Taxonomy summary from initial mapping

For each read in the initial mapping, the consensus taxonomy of top hits in
reference database was taken; these are summarized here at the requested taxon
level {stats["args"]["summary_taxlevel"]!s}. This should not be interpreted
as a direct measure of abundance, because the number of copies of rRNA genes
per genome is variable between species. Eukaryotes especially often have high
copy rRNA copy numbers.

{dict2markdowntable(stats['initial_taxonomy'], order_by_value=True, col1='Taxon', col2='Read count')}

## Assembled marker sequence clusters

Summary of assembled marker sequence clusters and their top hits in the
reference database. Cluster sequences with few underlying reads and many
mismatches/indels to the reference hits are lower quality because of
insufficient coverage, and should not be used for phylogenetics or probe
design. Clusters with a high number of underlying reads but low alignment
identity to a reference sequence may represent novel taxa, but should be
checked for sequence chimerism or misassembly.

{dod2markdowntable(per_cluster_summarize(stats), cluster_table_fields, col1='Cluster ID')}

<figure>

![]({kmercount_plot_path})

<figcaption>K-mer multiplicity plots for flanking sequences of each marker sequence cluster.</figcaption>
</figure>


---

phyloblitz depends on the following tools; please cite them:
[`minimap2`](https://github.com/lh3/minimap2) ([Li, 2018](https://doi.org/10.1093/bioinformatics/bty191)),
[`isONclust3`](https://github.com/aljpetri/isONclust3) ([Petri & Sahlin, 2025](https://doi.org/10.1093/bioinformatics/btaf207)),
[`pymarkovclustering`](https://github.com/moshi4/pyMarkovClustering),
[`mcl`](https://micans.org/mcl/) ([van Dongen, 2008](http://link.aip.org/link/?SJMAEL/30/121/1)),
[`pyfastx`](https://pyfastx.readthedocs.io/) ([Du, et al., 2020](https://doi.org/10.1093/bib/bbaa368)),
[`samtools`](https://www.htslib.org/) ([Li, Handsaker, et al., 2009](https://doi.org/10.1093/bioinformatics/btp352); [Bonfield, Marshall, Danecek, et al., 2021](https://doi.org/10.1093/gigascience/giab007)),
[`pysam`](https://github.com/pysam-developers/pysam),
[`spoa`](https://github.com/rvaser/spoa).

Please cite the [SILVA](https://www.arb-silva.de/) reference database
([Chuvochina, Gerken, et al., 2026](https://doi.org/10.1093/nar/gkaf1247)) if
you use it.

If you use `phyloblitz` in published research, please cite the GitHub
repository URL and software version.

HTML styling with [Simple.css](https://simplecss.org/)
"""

    format_md = create_markdown(renderer=MarkdownRenderer())
    return format_md(raw)


def generate_report_html(stats, histogram_file_path, kmercount_plot_path):
    """Generate HTML report from stats collected during phyloblitz run

    Same parameters as `generate_report_md`

    :returns: Report in HTML format
    :rtype: str
    """

    return HTML_TEMPLATE.replace(
        "{{markdown_report}}",
        html(generate_report_md(stats, histogram_file_path, kmercount_plot_path)),
    ).replace("{{title}}", "phyloblitz run report")


def per_cluster_summarize(stats):
    out = defaultdict(lambda: defaultdict(dict))
    for c in stats["cluster2seq"]:
        out["cluster_" + str(c)]["numseq"] = len(stats["cluster2seq"][c])
    for c in stats["cluster_tophits"]:
        out[c].update(stats["cluster_tophits"][c])
    if "cluster flanking numclust" in stats:
        for c in stats["cluster flanking numclust"]:
            out["cluster_" + str(c)]["flanking seq clusters"] = stats[
                "cluster flanking numclust"
            ][c]
    return out
