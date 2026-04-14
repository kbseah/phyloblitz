"""Compare phyloblitz runs."""

import json
import logging
from collections import defaultdict
from datetime import datetime
from multiprocessing import Pool
from pathlib import Path

from phyloblitz.__about__ import __version__
from phyloblitz.utils import (
    Pipeline,
    cluster_seqs_from_isonclust3,
    count_spoa_aln_persite_vars,
    count_spoa_aln_vars,
    parse_spoa_r2,
    run_isonclust3,
    run_md5,
    spoa_assemble_fasta,
)

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


class Compare(Pipeline):
    def __init__(self, args: dict) -> None:
        df = read_tsv(args["input_table"])
        samples = df["sample"]
        reports = df["report"]
        self._ref = args["db"]
        self._refindex = args["dbindex"]
        self._ref_md5 = run_md5(args["db"])
        self._reports = {}
        self._outdir = args["outdir"]
        self._prefix = args["prefix"]
        self._stats = {
            "version": __version__,
            "args": args,
            "runstats": {},
            "starttime": str(datetime.now()),
            "cluster2seq": {},
            "segment2sample": {},
            "cluster2sample": {},
        }
        logger.debug("Database checksum: %s", self._ref_md5)
        try:
            for sample, report in zip(samples, reports, strict=True):
                with Path.open(report, "r") as fh:
                    self._reports[sample] = json.load(fh)
            logger.info("%d reports read", len(samples))
        except ValueError as e:
            e.add_note("Number of sample names and report files do not agree")
            raise

    OUTFILE_SUFFIX = {
        "pooled_segments": "_segments.fastq",
        "isonclust3_cluster": "_isonclust3_out/clustering/final_clusters.tsv",
        "cluster_asm": "_final.fasta",
        "cluster_tophits": "_final_tophits.paf",
        "report_json": "_report.json",
        "report_md": "_report.md",
        "report_html": "_report.html",
    }

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
                if segment in self._stats["segment2sample"]:
                    logger.debug("Segment %s seen more than once", segment)
                    returned = False
                self._stats["segment2sample"][segment] = sample
        return returned

    def write_segments_to_fastq(self) -> None:
        """Write read segments to Fastq file for clustering."""
        fastq_path = self.pathto("pooled_segments")
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
        return run_isonclust3(
            self.pathto("pooled_segments"),
            mode,
            self.pathto("isonclust3_cluster").parent.parent,
        )

    def assemble_clusters(
        self,
        threads: int = 12,
        rseed: int = 12345,
        keeptmp: bool = False,
        min_clust_size: int = 5,
        max_clust_size: int = 500,
    ) -> None:
        """Extract cluster sequences and assemble with spoa.

        :param cluster_tool: Clustering method used, either "mcl" or "isonclust3".
        :param threads: Number of parallel assembly jobs to run.
        :param rseed: Random seed for downsampling reads in large clusters.
        :param keeptmp: If True, do not delete Fastq files with extracted reads.
        :param min_clust_size: Only assemble clusters containing at least this
            number of reads.
        :param max_clust_size: Downsample reads for clusters above this size.
        """
        fastq_handles, cluster2seq = cluster_seqs_from_isonclust3(
            self.pathto("isonclust3_cluster"),
            self.pathto("pooled_segments"),
            keeptmp=False,
            min_clust_size=min_clust_size,
            max_clust_size=max_clust_size,
            rseed=rseed,
        )
        self._stats["runstats"].update(
            {
                "number of clusters": len(cluster2seq),
                "number of clusters > 5 reads": len(
                    [i for i in cluster2seq if len(cluster2seq[i]) > 5],
                ),
                "total reads in clusters": sum(
                    [len(cluster2seq[c]) for c in cluster2seq],
                ),
            },
        )
        self._stats["cluster2seq"] = cluster2seq
        logger.info("Assemble consensus from clustered sequences with spoa")
        with Pool(threads) as pool:
            cluster_cons_tuples = pool.map(
                spoa_assemble_fasta,
                [(c, handle.name) for c, handle in fastq_handles.items()],
            )
        # Close NamedTemporaryFile handles
        for handle in fastq_handles.values():
            handle.close()

        cluster_cons = dict(cluster_cons_tuples)

        cluster_cons_parsed = {i: parse_spoa_r2(cluster_cons[i]) for i in cluster_cons}
        cluster_variant_counts = {
            i: count_spoa_aln_vars(cluster_cons_parsed[i]) for i in cluster_cons_parsed
        }
        cluster_persite_variant_counts = {  # WIP
            i: count_spoa_aln_persite_vars(cluster_cons_parsed[i])
            for i in cluster_cons_parsed
        }
        self._stats.update(
            {
                "cluster variant counts": cluster_variant_counts,
                "cluster persite variant counts": cluster_persite_variant_counts,  # WIP
                "cluster cons parsed": cluster_cons_parsed,
            },
        )
        with Path.open(self.pathto("cluster_asm"), "w") as fh:
            fh.writelines(
                f">cluster_{c!s} Consensus\n"
                + cluster_cons_parsed[c]["Consensus"].replace("-", "")
                + "\n"
                for c in cluster_cons_parsed
            )
            logger.info("Assembled sequences written to %s", self.pathto("cluster_asm"))

    def cluster_memberships(self) -> None:
        """Map cluster memberships to samples."""
        # Combine segment2sample and cluster2seq
        for cluster, seqs in self._stats["cluster2seq"].items():
            samples = defaultdict(list)
            for seq in seqs:
                if seq in self._stats["segment2sample"]:
                    samples[self._stats["segment2sample"][seq]].append(seq)
                else:
                    logger.warning(
                        "Sequence %s not found in segment2sample mapping",
                        seq,
                    )
            self._stats["cluster2sample"][cluster] = samples
