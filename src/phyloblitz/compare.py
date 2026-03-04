"""Compare phyloblitz runs."""

import json
import logging

from collections import defaultdict
from pathlib import Path

logger = logging.getLogger(__name__)


def read_tsv(filepath) -> dict:
    out = defaultdict(list)
    with Path.open(filepath, "r") as fh:
        header = next(fh)
        header = header.rstrip().split("\t")
        for line in fh:
            for key, val in zip(header, line.rstrip().split("\t"), strict=True):
                out[key].append(val)
    return out


class Compare:
    def __init__(self, infile: str | Path) -> None:
        df = read_tsv(infile)
        samples = df["sample"]
        reports = df["report"]
        self._reports = {}
        try:
            for sample, report in zip(samples, reports, strict=True):
                with Path.open(report, "r") as fh:
                    self._reports[sample] = json.load(fh)
            logger.info("%d reports read", len(samples))
        except ValueError as e:
            e.add_note("Number of sample names and report files do not agree")
            raise

    def check_database_checksums(self) -> None:
        try:
            if len({self._reports[s]["db_md5"] for s in self._reports}) > 1:
                logger.info("Mismatched database checksums found")
                for s in self._reports:
                    logger.debug("Sample %s checksum %s", s, self._reports[s]["db_md5"])
                return False
            else:
                logger.info("Checksums ok")
                return True
        except KeyError as e:
            e.add_note("Reports must be produced by phyloblitz v1.0.0 or higher")
            raise
