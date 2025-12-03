#!/usr/bin/env python3

from phyloblitz.__about__ import __version__
from mistune.renderers.markdown import MarkdownRenderer
from mistune import create_markdown, html


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
    out.append("| Cluster name | Number of reads |")
    out.append("| ------------ | --------------- |")
    for c in stats["cluster2seq"]:
        out.append(
            "| "
            + "cluster_"
            + str(c)
            + " | "
            + str(len(stats["cluster2seq"][c]))
            + " |"
        )
    out.append("")

    raw = "\n".join(out)
    format_md = create_markdown(renderer=MarkdownRenderer())
    return format_md(raw)


def generate_report_html(stats):
    """Generate HTML report from stats collected during phyloblitz run"""

    return html(generate_report_md(stats))
