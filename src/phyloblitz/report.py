#!/usr/bin/env python3

import re
import logging

import numpy as np
import matplotlib.pyplot as plt

from phyloblitz.__about__ import __version__
import phyloblitz.main as main
from mistune.renderers.markdown import MarkdownRenderer
from mistune import create_markdown, html
from collections import defaultdict

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="%(levelname)s : %(module)s : %(asctime)s : %(message)s",
    level=logging.DEBUG,
    datefmt="%Y-%m-%d %H:%M:%S",
)

CIGAROPS = {
    "M": "match",
    "I": "insertion",  # to reference
    "D": "deletion",  # to reference
    "N": "skipped",  # region to reference
    "S": "soft clipping",
    "H": "hard clipping",
    "P": "padding",
    "=": "seq match",
    "X": "seq mismatch",
}


def dict2markdowntable(d, keys=None, col1="Parameter", col2="Value"):
    """Convert dict to a markdown table

    Key will be cast as row names

    :param d: Dict to convert to table
    :param keys: Keys to include in table; if None, all keys will be used
    :param col1: Column header for first column
    :param col2: Column header for second column
    :returns: Text table in markdown format
    :rtype: str
    """
    if keys is None:
        keys = list(d.keys())
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

    if not main.check_run_file(args, "report_dvs_hist"):
        dvs = [min(stats["dvs"][r]) for r in stats["dvs"]]
        stats["runstats"]["ava min dvs median"] = "{:.4f}".format(np.median(dvs))
        fig, axs = plt.subplots(1, figsize=(3, 2))
        axs.hist(dvs, bins="auto", density=True)
        axs.axvline(args.dv_max, color="red")
        axs.set_title("Histogram of ava min dvs")
        fig.tight_layout()
        fig.savefig(main.pathto(args, "report_dvs_hist"))


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

![]({main.pathto(args, 'report_dvs_hist', basename_only=True)})

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


def db_taxonomy(silva_fasta):
    """Get taxonomy string from SILVA headers in database Fasta file

    :param silva_fasta: Path to SILVA reference database in Fasta format; headers must follow SILVA format
    :returns: Dict of taxonomy strings keyed by SILVA accession
    :rtype: dict
    """
    acc2tax = {}
    with open(silva_fasta, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                spl = line.lstrip(">").rstrip().split(" ")
                acc = spl[0]
                taxstring = " ".join(spl[1:]).split(";")
                acc2tax[acc] = taxstring
    return acc2tax


def parse_cigar_ops(cigar):
    """Summarize operations in a CIGAR string"""
    summary = defaultdict(int)
    cigar = re.findall(r"(\d+)(\D)", cigar)
    for num, op in cigar:
        summary[op] += int(num)
    return summary


def lists_common_prefix(lol):
    """Get common prefix in a list of lists

    :param lol: list of lists of strings
    :returns: list of the common prefix
    :rtype: list
    """
    out = []
    for j in range(min([len(l) for l in lol])):
        s = set([i[j] for i in lol])
        if len(s) == 1:
            out.append(s.pop())
        else:
            break
    return out


def summarize_tophit_paf(paf_file, silva_fasta, acc2tax):
    """Summarize top hits of assembled seqs mapped to SILVA database by minimap2

    :param paf_file: Path to PAF output from minimap2, must have CIGAR string
    :param silva_fasta: Path to Fasta formatted SILVA database
    :param acc2tax: Dict of SILVA taxonomy strings keyed by ref db sequence id
    :returns: Dict of summary stats for each hit, keyed by query sequence name
    :rtype: dict
    """
    out = {}

    with open(paf_file, "r") as fh:
        for line in fh:
            spl = line.rstrip().split("\t")
            hits = dict(
                zip(
                    [
                        "qname",
                        "qlen",
                        "qstart",
                        "qend",
                        "strand",
                        "tname",
                        "tlen",
                        "tstart",
                        "tend",
                        "alnmatch",
                        "alnlen",
                    ],
                    spl[0:11],
                )
            )
            cigar = [i for i in spl if i.startswith("cg:Z:")][0]
            # Calculate derived metrics from PAF fields
            cigar_summary = parse_cigar_ops(cigar[5:])
            hits.update({CIGAROPS[c]: cigar_summary[c] for c in cigar_summary})
            hits["align %id"] = "{:.2%}".format(
                int(hits["alnmatch"]) / int(hits["alnlen"])
            ).rstrip(
                "%"
            )  # remove redundant % sign for display
            hits["query %aln"] = "{:.2%}".format(
                (int(hits["qend"]) - int(hits["qstart"])) / int(hits["qlen"])
            ).rstrip("%")
            hits["target %aln"] = "{:.2%}".format(
                (int(hits["tend"]) - int(hits["tstart"])) / int(hits["tlen"])
            ).rstrip("%")
            out[spl[0]] = hits

    # Taxonomy of hit targets
    for c in out:
        try:
            # hyperlink to ENA record
            out[c]["tophit"] = (
                "["
                + out[c]["tname"]
                + "](https://www.ebi.ac.uk/ena/browser/view/"
                + out[c]["tname"].split(".")[0]
                + ")"
            )
            out[c]["tophit taxonomy"] = ";".join(acc2tax[out[c]["tname"]])
            out[c]["tophit species"] = acc2tax[out[c]["tname"]][-1]
            # Higher taxonomy to class level, except for chloroplast and mitochondria
            # Assumes SILVA taxonomy is in use
            if out[c]["tophit taxonomy"].startswith(
                "Bacteria;Cyanobacteria;Cyanobacteriia;Chloroplast"
            ):
                out[c]["higher taxonomy"] = "[Chloroplast]"
            elif out[c]["tophit taxonomy"].startswith(
                "Bacteria;Proteobacteria;Alphaproteobacteria;Rickettsiales;Mitochondria"
            ):
                out[c]["higher taxonomy"] = "[Mitochondria]"
            else:
                out[c]["higher taxonomy"] = ";".join(acc2tax[out[c]["tname"]][0:3])
        except KeyError:
            KeyError(f"Accession {out[c]['tname']} not found in database?")
    return out


def per_cluster_summarize(stats):
    out = defaultdict(lambda: defaultdict(dict))
    for c in stats["cluster2seq"]:
        out["cluster_" + str(c)]["numseq"] = len(stats["cluster2seq"][c])
    for c in stats["cluster_tophits"]:
        out[c].update(stats["cluster_tophits"][c])
    return out
