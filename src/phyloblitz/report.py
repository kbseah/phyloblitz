#!/usr/bin/env python3

import re
import logging

from phyloblitz.__about__ import __version__
from mistune.renderers.markdown import MarkdownRenderer
from mistune import create_markdown, html
from collections import defaultdict

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="%(levelname)s: %(asctime)s : %(message)s",
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


def generate_report_md(stats):
    """Generate markdown report from stats collected during phyloblitz run"""

    out = []
    out.append("# phyloblitz run report\n")
    out.append("* Run started: " + str(stats["starttime"]))
    out.append("* Run ended: " + str(stats["endtime"]))
    out.append("* phyloblitz version: " + __version__)
    out.append("")
    out.append("phyloblitz [homepage](https://github.com/kbseah/phyloblitz)")
    out.append("")
    out.append("## Argument string\n")
    out.append("```")
    out.append(
        " ".join([" ".join(["--" + i, str(stats["args"][i])]) for i in stats["args"]])
    )
    out.append("```")
    out.append("")
    out.append("## Cluster stats\n")
    cluster_table_dict = per_cluster_summarize(stats)
    cluster_table_fields = [
        "numseq",
        "qstart",
        "qend",
        "tname",
        "tophit taxonomy",
        "tstart",
        "tend",
        "seq match",
        "seq mismatch",
        "insertion",
        "deletion",
    ]
    out.append("| cluster name | " + " | ".join(cluster_table_fields) + " |")
    out.append(
        "| ------ | " + " | ".join(["------" for f in cluster_table_fields]) + " |"
    )
    for c in cluster_table_dict:
        out.append(
            "| "
            + c
            + " | "
            + " | ".join(
                [
                    str(cluster_table_dict[c][k]) if k in cluster_table_dict[c] else "-"
                    for k in cluster_table_fields
                ]
            )
            + " |"
        )
    out.append("")

    raw = "\n".join(out)
    format_md = create_markdown(renderer=MarkdownRenderer())
    return format_md(raw)


def generate_report_html(stats):
    """Generate HTML report from stats collected during phyloblitz run"""

    return html(generate_report_md(stats))


def db_taxonomy(silva_fasta):
    """Get taxonomy string from SILVA headers in database Fasta file"""
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


def summarize_tophit_paf(paf_file, silva_fasta):
    out = {}

    with open(paf_file, "r") as fh:
        for line in fh:
            spl = line.rstrip().split("\t")
            hits = dict(
                zip(
                    ["qname", "qstart", "qend", "tname", "tstart", "tend"],
                    [spl[0], spl[2], spl[3], spl[5], spl[7], spl[8]],
                )
            )
            cigar = [i for i in spl if i.startswith("cg:Z:")][0]
            cigar_summary = parse_cigar_ops(cigar[5:])
            hits.update({CIGAROPS[c]: cigar_summary[c] for c in cigar_summary})
            out[spl[0]] = hits

    logger.info("Reading taxonomy from SILVA database file")
    acc2tax = db_taxonomy(silva_fasta)
    logger.debug(f" Accessions read: {str(len(acc2tax))}")
    for c in out:
        try:
            out[c]["tophit taxonomy"] = ";".join(acc2tax[out[c]["tname"]])
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
