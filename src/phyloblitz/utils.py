#!/usr/bin/env python3

import logging
import re
import os.path

from collections import defaultdict

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

