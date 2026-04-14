"""Compare phyloblitz runs."""

import json
import logging
from collections import defaultdict
from pathlib import Path

from phyloblitz.utils import cluster_seqs_from_isonclust3, run_isonclust3, run_md5

logger = logging.getLogger(__name__)


def read_tsv(filepath: str | Path) -> dict[list]:
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
    def __init__(self, args: dict) -> None:
        df = read_tsv(args["input_table"])
        samples = df["sample"]
        reports = df["report"]
        self._ref = args["db"]
        self._refindex = args["dbindex"]
        self._ref_md5 = run_md5(args["db"])
        self._reports = {}
        self._segment2sample = {}
        self._outdir = args["outdir"]
        self._prefix = args["prefix"]
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
                logger.error("Mismatched database checksums found")
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

    def check_sequencing_platforms(self) -> bool:
        """Check whether the same sequencing platform was used for runs to be compared.

        Warn user if different sequencing platforms were used in the phyloblitz
        runs to be compared. Clustering won't work properly if different
        platforms were used, so this is a hard check.

        :returns: True if checks are OK, else False
        :rtype: bool
        """
        try:
            platforms = {self._reports[s]["args"]["platform"] for s in self._reports}
            if len(platforms) > 1:
                logger.error("Mismatched sequencing platforms found")
                for s in self._reports:
                    logger.debug(
                        "Sample %s platform %s",
                        s,
                        self._reports[s]["args"]["platform"],
                    )
                return False
            self._platform = next(iter(platforms))
            logger.info(
                "Results derived from same sequencing platform %s",
                self._platform,
            )
            return True
        except KeyError as e:
            e.add_note("Reports must be produced by phyloblitz v1.0.0 or higher")
            raise

    def segment_by_sample(self) -> bool:
        """Create mapping of read segments to samples.

        Needed for clustering and for comparing the clusters across samples.
        Assumes that the same read segment will not be seen in more than one
        sample, else reports an error.
        """
        returned = True
        for sample, report in self._reports.items():
            for segment in report["segments"]:
                if segment in self._segment2sample:
                    logger.debug("Segment %s seen more than once", segment)
                    returned = False
                self._segment2sample[segment] = sample
        return returned

    def write_segments_to_fastq(self) -> None:
        """Write read segments to Fastq file for clustering."""
        fastq_path = Path(self._outdir) / Path(self._prefix + "_segments.fastq")
        # If file already exists, remove it to avoid appending to an old file
        if fastq_path.exists():
            logger.warning(
                "Fastq file %s already exists, it will be overwritten",
                fastq_path,
            )
            fastq_path.unlink()
        with Path.open(fastq_path, "a") as fh:
            for report in self._reports.values():
                for seg_name, seg_dict in report["segments"].items():
                    fastq_rec = (
                        "@"
                        + seg_name
                        + "\n"
                        + seg_dict["seq"]
                        + "\n+\n"
                        + seg_dict["quals"]
                        + "\n"
                    )
                    fh.write(fastq_rec)

    def cluster_segments(self) -> int:
        """Cluster read segments with isONclust3."""
        if self._platform in ["lr:hq", "map-ont"]:
            mode = "ont"
        elif self._platform in ["map-hifi", "map-pb"]:
            mode = "pacbio"
        fastq_path = Path(self._outdir) / Path(self._prefix + "_segments.fastq")
        outfolder = Path(self._outdir) / Path(self._prefix + "_clusters")
        return run_isonclust3(fastq_path, mode, outfolder)

    def cluster_memberships(self) -> None:
        """Map cluster memberships to samples."""
        fastq_path = Path(self._outdir) / Path(self._prefix + "_segments.fastq")
        isonclust3_out = (
            Path(self._outdir)
            / Path(self._prefix + "_clusters")
            / Path("clustering")
            / Path("final_clusters.tsv")
        )
        _fastq_handles, cluster2seq = cluster_seqs_from_isonclust3(
            isonclust3_out,
            fastq_path,
            keeptmp=False,
            min_clust_size=5,
        )
        # Combine segment2sample and cluster2seq
        self._cluster2sample = {}
        for cluster, seqs in cluster2seq.items():
            samples = defaultdict(list)
            for seq in seqs:
                if seq in self._segment2sample:
                    samples[self._segment2sample[seq]].append(seq)
                else:
                    logger.warning(
                        "Sequence %s not found in segment2sample mapping",
                        seq,
                    )
            self._cluster2sample[cluster] = samples
