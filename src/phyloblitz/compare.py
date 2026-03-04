"""Compare phyloblitz runs."""

import json
import logging

from collections import defaultdict
from hashlib import md5
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
    def __init__(self, db: str | Path, infile: str | Path) -> None:
        df = read_tsv(infile)
        samples = df["sample"]
        reports = df["report"]
        self._ref = db
        self._ref_md5 = ""
        self._reports = {}
        try:
            for sample, report in zip(samples, reports, strict=True):
                with Path.open(report, "r") as fh:
                    self._reports[sample] = json.load(fh)
            logger.info("%d reports read", len(samples))
        except ValueError as e:
            e.add_note("Number of sample names and report files do not agree")
            raise

    def checksum_db(self) -> None:
        """Checksum database file for validation when comparing samples.

        Get path to reference database from Pipeline._ref attribute. Write
        checksum in Compare._ref_md5
        """
        md5_hash = md5()
        with Path.open(self._ref, "rb") as f:
            for byte_block in iter(lambda: f.read(4096), b""):
                md5_hash.update(byte_block)
        checksum = md5_hash.hexdigest()
        logger.debug("Database file %s has MD5 checksum %s", self._ref, checksum)
        self._ref_md5 = checksum

    def check_database_checksums(self) -> None:
        try:
            reports_md5 = {self._reports[s]["db_md5"] for s in self._reports}
            if len(reports_md5) > 1:
                logger.info("Mismatched database checksums found")
                for s in self._reports:
                    logger.debug("Sample %s checksum %s", s, self._reports[s]["db_md5"])
                return False
            else:
                logger.info(
                    "Results were derived from phyloblitz runs against the same database"
                )
                if list(reports_md5)[0] != self._ref_md5:
                    logger.warning(
                        "Database used to generate phyloblitz runs has different checksum %s from database used as reference for phyloblitz compare %s",
                        list(reports_md5)[0],
                        self._ref_md5,
                    )
                else:
                    logger.info(
                        "Checksums of database %s and that used to generate reports match",
                        self._ref,
                    )
                    return True
        except KeyError as e:
            e.add_note("Reports must be produced by phyloblitz v1.0.0 or higher")
            raise
