"""Compare phyloblitz runs."""

import json
import logging
from collections import defaultdict
from datetime import datetime
from multiprocessing import Pool
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from mistune import create_markdown, html
from mistune.renderers.markdown import MarkdownRenderer
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

from phyloblitz.__about__ import __version__
from phyloblitz.report import HTML_TEMPLATE, dict2markdowntable, dod2markdowntable
from phyloblitz.utils import (
    Pipeline,
    check_stage_file,
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


# TODO: User specify clustering and distance methods
def clustermap_rows_cols(
    X,
    *,
    row_labels: list | None,
    col_labels: list | None,
    row_metric: str = "euclidean",
    col_metric: str = "euclidean",
    row_method: str = "ward",
    col_method: str = "ward",
    cmap: str = "viridis",
    figsize: tuple = (10, 10),
    xtick_rotation: float = 90,
    xtick_fontsize: float = 8,
    ytick_fontsize: float = 8,
) -> tuple:
    """Plot heatmap and dendrograms of 2-d array.

    Hierarchical clustering of rows and columns of X, then plot heatmap with
    row/col dendrograms (scipy + matplotlib only), with row/col labels.

    :param X: (nrow, ncol) array-like
    :param row_labels: Labels for rows (length nrow). If None, uses 0..nrow-1.
    :param col_labels: Labels for cols (length ncol). If None, uses 0..ncol-1.
    :param row_metric: Distance metric for clustering rows (passed to scipy.pdist).
    :param col_metric: Distance metric for clustering columns (passed to scipy.pdist).
    :param row_method: Linkage method for clustering rows (passed to scipy.linkage).
    :param col_method: Linkage method for clustering columns (passed to scipy.linkage).
    :param cmap: Colormap for heatmap (passed to matplotlib.pyplot.imshow).
    :param figsize: Figure size (passed to matplotlib.pyplot.figure).
    :param xtick_rotation: Rotation angle for x-axis tick labels (column labels).
    :param xtick_fontsize: Font size for x-axis tick labels (column labels).
    :param ytick_fontsize: Font size for y-axis tick labels (row labels).
    :returns: Figure and dict with data.
    """
    X = np.asarray(X)
    nrow, ncol = X.shape

    if row_labels is None:
        row_labels = [str(i) for i in range(nrow)]
    if col_labels is None:
        col_labels = [str(j) for j in range(ncol)]
    if len(row_labels) != nrow:
        raise ValueError(f"row_labels must have length {nrow}, got {len(row_labels)}")
    if len(col_labels) != ncol:
        raise ValueError(f"col_labels must have length {ncol}, got {len(col_labels)}")

    # cluster rows
    row_dist = pdist(X, metric=row_metric)
    row_link = linkage(row_dist, method=row_method)
    row_den = dendrogram(row_link, no_plot=True)
    row_order = row_den["leaves"]

    # cluster columns (cluster the transpose)
    col_dist = pdist(X.T, metric=col_metric)
    col_link = linkage(col_dist, method=col_method)
    col_den = dendrogram(col_link, no_plot=True)
    col_order = col_den["leaves"]

    # reorder matrix and labels
    Xo = X[np.ix_(row_order, col_order)]
    row_labels_o = [row_labels[i] for i in row_order]
    col_labels_o = [col_labels[j] for j in col_order]

    # layout: row dendrogram (left), col dendrogram (top), heatmap (center), colorbar (right)
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(
        nrows=2,
        ncols=3,
        width_ratios=(1.2, 6.0, 0.3),
        height_ratios=(1.2, 6.0),
        wspace=0.05,
        hspace=0.05,
    )

    ax_col = fig.add_subplot(gs[0, 1])
    ax_row = fig.add_subplot(gs[1, 0])
    ax_hm = fig.add_subplot(gs[1, 1])
    ax_cb = fig.add_subplot(gs[1, 2])

    # Column dendrogram (top)
    dendrogram(
        col_link,
        ax=ax_col,
        orientation="top",
        no_labels=True,
        color_threshold=0,
        above_threshold_color="black",
    )
    ax_col.set_xticks([])
    ax_col.set_yticks([])
    for sp in ax_col.spines.values():
        sp.set_visible(False)

    # Row dendrogram (left)
    dendrogram(
        row_link,
        ax=ax_row,
        orientation="left",
        no_labels=True,
        color_threshold=0,
        above_threshold_color="black",
    )
    ax_row.set_xticks([])
    ax_row.set_yticks([])
    for sp in ax_row.spines.values():
        sp.set_visible(False)

    # Heatmap
    im = ax_hm.imshow(
        Xo,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        origin="upper",
    )

    # Set ticks at cell centers
    ax_hm.set_xticks(np.arange(ncol))
    ax_hm.set_yticks(np.arange(nrow))

    ax_hm.set_xticklabels(
        col_labels_o,
        rotation=xtick_rotation,
        ha="right",
        fontsize=xtick_fontsize,
    )
    ax_hm.set_yticklabels(row_labels_o, fontsize=ytick_fontsize)

    # Put x labels at bottom only
    ax_hm.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)

    # Optional: small tick marks off for a cleaner heatmap
    ax_hm.tick_params(axis="both", which="both", length=0)

    # Clean spines
    for sp in ax_hm.spines.values():
        sp.set_visible(False)

    # Colorbar
    fig.colorbar(im, cax=ax_cb)

    return fig, {
        "row_linkage": row_link,
        "col_linkage": col_link,
        "row_order": row_order,
        "col_order": col_order,
        "X_reordered": Xo,
        "row_labels_reordered": row_labels_o,
        "col_labels_reordered": col_labels_o,
    }


class Compare(Pipeline):
    """phyloblitz compare pipeline object."""

    def __init__(self, args: dict) -> None:
        """Construct Run object."""
        df = read_tsv(args["input_table"])
        self._samples = df["sample"]
        self._report_files = df["report"]
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
            for sample, report in zip(self._samples, self._report_files, strict=True):
                with Path.open(report, "r") as fh:
                    self._reports[sample] = json.load(fh)
            logger.info("%d reports read", len(self._samples))
            self._stats["runstats"]["number of samples compared"] = len(self._samples)
        except ValueError as e:
            e.add_note("Number of sample names and report files do not agree")
            raise

    OUTFILE_SUFFIX = {
        "pooled_segments": "_segments.fastq",
        "isonclust3_cluster": "_isonclust3_out/clustering/final_clusters.tsv",
        "cluster_asm": "_final.fasta",
        "cluster_tophits": "_final_tophits.paf",
        "cluster_membership_heatmap": "_cluster_membership_heatmap.png",
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

    @check_stage_file(
        stage="pooled_segments",
        message="Pool read segments from all samples for clustering",
    )
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
                    fh.write(
                        f"@{seg_name!s}\n{seg_dict['seq']!s}\n+\n{seg_dict['quals']!s}\n",
                    )

    @check_stage_file(
        stage="isonclust3_cluster",
        message="Clustering with isonclust3",
    )
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

    @check_stage_file(
        stage="cluster_asm",
        message="Extract cluster sequences and assemble with spoa",
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

    # TODO: Classify assembled cluster sequences against reference database

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

    @check_stage_file(
        stage="report_json",
        message="Writing report stats in JSON format",
    )
    def write_report_json(self) -> None:
        """Dump run stats file in JSON format.

        Dump the Run._stats attribute to a JSON file for troubleshooting
        and comparison of different phyloblitz runs.
        """
        self._stats.update({"endtime": str(datetime.now())})
        with Path.open(self.pathto("report_json"), "w") as fh:
            json.dump(self._stats, fh, indent=4)

    @check_stage_file(
        stage="report_md",
        message="Writing report in Markdown format",
    )
    @check_stage_file(
        stage="report_html",
        message="Writing report in HTML format",
    )
    def write_reports(self) -> None:
        """Write report in Markdown and HTML formats."""
        # Input table of sample and report file paths
        input_table_md = dict2markdowntable(
            dict(zip(self._samples, self._report_files, strict=True)),
            col1="Sample",
            col2="Report file",
        )

        # Count reads per cluster per sample
        cluster2sample_counts = {
            c: {
                s: len(self._stats["cluster2sample"][c][s])
                for s in self._stats["cluster2sample"][c]
            }
            for c in self._stats["cluster2sample"]
        }
        counts_md = dod2markdowntable(
            cluster2sample_counts,
            keys=self._samples,
            order_by=self._samples[0],
            order_by_value=False,
            col1="Cluster",
            fill_empty="0",
        )
        # TODO: Use percentages instead of counts for visualization
        # Generate dendrogram and heatmap of cluster memberships across samples
        # Rows: Clusters, Columns: Samples
        cluster2sample_array = np.array(
            [
                [cluster2sample_counts[c].get(s, 0) for s in self._samples]
                for c in cluster2sample_counts
            ],
        )
        # Normalize by column (sample) sums to get relative abundances
        cluster2sample_array = cluster2sample_array / cluster2sample_array.sum(axis=0, keepdims=True)
        fig, out = clustermap_rows_cols(
            cluster2sample_array,
            row_labels=[f"Cluster {c}" for c in cluster2sample_counts],
            col_labels=self._samples,
            row_method="ward",
            col_method="ward",
            cmap="viridis",
            figsize=(10, max(6, len(cluster2sample_counts) * 0.5)),
        )
        fig.savefig(self.pathto("cluster_membership_heatmap"), bbox_inches="tight")

        raw = f"""# phyloblitz compare report

* Compare started: {self._stats["starttime"]!s}
* phyloblitz version: {__version__}

phyloblitz [homepage](https://github.com/kbseah/phyloblitz)


## Input parameters

`phyloblitz compare` was called with the following parameters:

{dict2markdowntable(self._stats["args"])}


### Input files

{input_table_md}


## Run statistics

{dict2markdowntable(self._stats["runstats"])}


## Co-assembled marker sequence clusters

<figure>

![]({self.pathto("cluster_membership_heatmap", basename_only=True)!s})

<figcaption>Cluster membership heatmap. Rows are clusters, columns are samples,
and values are relative read abundance from each sample per sequence cluster.
Clusters and samples are ordered by hierarchical clustering with Ward's
linkage and Euclidean distance. Clusters below minimum cluster size for
assembly are not shown.</figcaption>
</figure>


### Read counts per sequence cluster per sample

{counts_md}
"""

        format_md = create_markdown(renderer=MarkdownRenderer())
        out_md = format_md(raw)

        with Path.open(self.pathto("report_md"), "w") as fh:
            fh.write(out_md)

        out_html = HTML_TEMPLATE.replace(
            "{{markdown_report}}",
            html(out_md),
        ).replace("{{title}}", "phyloblitz compare report")

        with Path.open(self.pathto("report_html"), "w") as fh:
            fh.write(out_html)
