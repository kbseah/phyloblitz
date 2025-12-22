#!/usr/bin/env python3

import logging
import re
import os.path

from collections import defaultdict
from subprocess import Popen, PIPE, STDOUT

logger = logging.getLogger(__name__)

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


def filter_paf_overhang(line: str, max_overhang_frac: float = 0.05):
    """Filter out PAF alignments with incompatible overhangs

    If two aligned reads have overhangs that do not align, and the overhangs
    are on the same side of the alignment, this means that the underlying
    sequences have a conserved homologous region (e.g. a repeat or conserved
    homolog), flanked by at least one non-homologous region.

    :param line: Single line from PAF entry
    :param max_overhang_frac: Max fraction of read length that same-side
        overhang is allowed to be
    :returns: Input `line` if alignment does not have overhangs on same side(s)
        of alignment.
    """
    [qlen, qstart, qend, tlen, tstart, tend] = [
        int(line.rstrip().split("\t")[i]) for i in [1, 2, 3, 6, 7, 8]
    ]
    s_q = qstart / qlen
    s_t = tstart / tlen
    t_q = (qlen - qend) / qlen
    t_t = (tlen - tend) / tlen
    if (s_q > max_overhang_frac and s_t > max_overhang_frac) or (
        t_q > max_overhang_frac and t_t > max_overhang_frac
    ):
        # incompatible overlap
        pass
    else:
        return line


def check_dependencies():
    vers = {}
    from sys import version as python_version
    from pysam import __version__ as pysam_version
    from pyfastx import __version__ as pyfastx_version
    from mistune import __version__ as mistune_version
    from numpy import __version__ as numpy_version
    from matplotlib import __version__ as matplotlib_version

    vers["python"] = python_version
    vers["pysam"] = pysam_version
    vers["pyfastx"] = pyfastx_version
    vers["mistune"] = mistune_version
    vers["numpy"] = numpy_version
    vers["matplotlib"] = matplotlib_version

    for tool in ["minimap2", "spoa", "mcl", "isONclust3"]:
        p = Popen([tool, "--version"], stdout=PIPE, stderr=STDOUT, text=True)
        vers[tool] = (
            p.communicate()[0].split("\n")[0].rstrip()
        )  # mcl has multiline output
    return vers
