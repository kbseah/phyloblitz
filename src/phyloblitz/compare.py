"""Compare phyloblitz runs."""

import json
import logging
from collections import defaultdict
from pathlib import Path

from phyloblitz.utils import run_md5

logger = logging.getLogger(__name__)


def read_tsv(filepath:str|Path) -> dict[list]:
    """Read TSV file to dict oriented by column.

    :param filepath: Path to TSV file to parse.
    :returns: dict with column headers as keys and lists of column data as values.
    :rtype: dict
    """
    out = defaultdict(list)
    with Path.open(filepath, "r") as fh:
        header = next(fh)
        header = header.rstrip().split("\t")
        for line in fh:
            for key, val in zip(header, line.rstrip().split("\t"), strict=True):
                out[key].append(val)
    return out


class Compare:
    # def __init__(self, db: str | Path, infile: str | Path, outdir: str|Path, prefix:str) -> None:
    def __init__(self, db: str | Path, infile: str | Path) -> None:
        df = read_tsv(infile)
        samples = df["sample"]
        reports = df["report"]
        self._ref = db
        self._ref_md5 = run_md5(db)
        self._reports = {}
        self._segment2sample = {}
        logger.debug("Database checksum: %s", self._ref_md5)
        try:
            for sample, report in zip(samples, reports, strict=True):
                with Path.open(report, "r") as fh:
                    self._reports[sample] = json.load(fh)
            logger.info("%d reports read", len(samples))
        except ValueError as e:
            e.add_note("Number of sample names and report files do not agree")
            raise

    def check_database_checksums(self) -> bool:
        """Check whether the same database was used for runs to be compared.

        Warn user if different databases were used in the phyloblitz runs to be
        compared. Also check if the database for classifying the assembled
        clusters downstream is the same as the one used to extract the read
        segments. Can be overridden if this is intentional.

        :returns: True if checks are OK, else False
        :rtype: bool
        """
        try:
            reports_md5 = {self._reports[s]["db_md5"] for s in self._reports}
            if len(reports_md5) > 1:
                logger.info("Mismatched database checksums found")
                for s in self._reports:
                    logger.debug("Sample %s checksum %s", s, self._reports[s]["db_md5"])
                return False
            logger.info(
                "Results derived from phyloblitz runs against the same database",
            )
            if next(iter(reports_md5)) != self._ref_md5:
                logger.warning(
                    "Database used to generate phyloblitz runs has different checksum %s from database used as reference for phyloblitz compare %s",
                    next(iter(reports_md5)),
                    self._ref_md5,
                )
                return False
            logger.info(
                "Checksums of database %s and that used to generate reports match",
                self._ref,
            )
            return True
        except KeyError as e:
            e.add_note("Reports must be produced by phyloblitz v1.0.0 or higher")
            raise
