"""Shared utility functions for phyloblitz."""

import json
import logging
import re
import sys
from collections import defaultdict
from datetime import datetime
from functools import wraps
from hashlib import md5
from multiprocessing import Pool
from os import W_OK, access
from pathlib import Path
from random import sample, seed
from subprocess import PIPE, STDOUT, Popen
from tempfile import NamedTemporaryFile

import numpy as np
import pyfastx
from Bio import AlignIO
from Bio import __version__ as biopython_version
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import write as phylo_write
from matplotlib import __version__ as matplotlib_version
from mistune import __version__ as mistune_version
from pymarkovclustering import __version__ as pymcl_version
from pysam import __version__ as pysam_version

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


def parse_cigar_ops(cigar: str) -> dict:
    """Summarize operations in a CIGAR string."""
    summary = defaultdict(int)
    cigar = re.findall(r"(\d+)(\D)", cigar)
    for num, op in cigar:
        summary[op] += int(num)
    return summary


def lists_common_prefix(lol: list[list]) -> list:
    """Get common prefix in a list of lists.

    :param lol: list of lists of strings
    :returns: list of the common prefix
    :rtype: list
    """
    out = []
    for j in range(min([len(l) for l in lol])):
        s = {i[j] for i in lol}
        if len(s) == 1:
            out.append(s.pop())
        else:
            break
    return out


def filter_paf_overhang(line: str, max_overhang_frac: float = 0.05) -> str | None:
    """Filter out PAF alignments with incompatible overhangs.

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
        return None
    return line


def check_dependencies() -> dict:
    """Check if depdendencies present and get versions.

    Report dependency versions in log files for reproducibility and to
    encourage users to cite them.

    :return: dict of dependency versions keyed by name
    :rtype: dict
    """
    vers = dict(
        zip(
            [
                "python",
                "pysam",
                "pyfastx",
                "mistune",
                "numpy",
                "matplotlib",
                "pymarkovclustering",
            ],
            [
                sys.version,
                pysam_version,
                pyfastx.__version__,
                mistune_version,
                np.__version__,
                matplotlib_version,
                pymcl_version,
            ],
            strict=True,
        ),
    )
    for tool in ["minimap2", "spoa", "isONclust3"]:
        p = Popen([tool, "--version"], stdout=PIPE, stderr=STDOUT, text=True)
        # Split in case version message is multiline, e.g. mcl
        vers[tool] = p.communicate()[0].split("\n")[0].rstrip()
    return vers


def check_outdir(outdir: str | Path, resume: bool = True) -> None:
    """Check if output directory exists, and create it if it doesn't.

    Raise exceptions if output path exists but is not a directory, or if it is
    not writable. If output directory already exists and resume is False, raise
    an exception. Otherwise create output directory recursively if it doesn't
    exist, or simply carry on if it exists and resume is True.

    :param outdir: Path to output directory
    :param resume: If True, allow existing output directory to be used for resuming
    """
    outdir = Path(outdir)
    if not Path.exists(outdir):
        Path.mkdir(outdir, parents=True, exist_ok=False)
        logger.debug("Created output directory: %s", outdir)
    elif not Path.is_dir(outdir):
        msg = f"Output path {outdir!s} exists but is not a directory."
        raise NotADirectoryError(msg)
    # outdir exists but is not writable
    elif not access(outdir, W_OK):
        msg = f"Output directory {outdir!s} is not writable."
        raise PermissionError(msg)
    elif resume:
        logger.info("Output directory %s already exists, resuming.", outdir)
    else:
        msg = f"Output directory {outdir!s} already exists, but resume is False."
        raise FileExistsError(msg)


def run_md5(file: str | Path) -> str:
    """Calculate MD5 hash for a file.

    :param file: Path to file.
    :returns: MD5 hash hex digest.
    :rtype: str
    """
    md5_hash = md5()
    with Path.open(file, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            md5_hash.update(byte_block)
    return md5_hash.hexdigest()


def cluster_seqs_from_isonclust3(
    isonclust3_out: Path,
    reads: Path,
    keeptmp: bool,
    min_clust_size: int = 5,
    max_clust_size: int = 1000,
    rseed: int = 12345,
) -> tuple:
    """Extract sequences from each isONclust3 cluster to Fastq files.

    The numbering of clusters by isONclust3 itself is not consistent
    between the final_clusters.tsv file and the output fastq files, so we
    extract the sequences ourselves.

    :param isonclust3_out: Path to isONclust3 output file
    :param reads: Path to reads from extract_reads_for_ava step
    :param keeptmp: Keep temporary files?
    :param min_clust_size: Minimum cluster size to return Fastq file
    :param max_clust_size: Downsample clusters if above this cluster size
    :param rseed: Random seed for downsampling clusters above max_clust_size
    :returns: Dict of file handles to each Fastq file keyed by cluster ID;
        excludes clusters below min_clust_size, and reads for clusters above
        max_clust_size are downsampled to max_clust_size.
    :returns: Dict of sequence IDs keyed by cluster ID; includes all sequences
        regardless of cluster size.
    :rtype: tuple
    """
    fastq_handles = {}
    seq2cluster = {}
    with Path.open(isonclust3_out) as fh:
        for line in fh:
            clust, seqname = line.rstrip().split("\t")
            seq2cluster[seqname] = clust
    # Cluster memberships for all read segments regardless of cluster size
    cluster2seq = defaultdict(list)
    for seqname in seq2cluster:
        cluster2seq[seq2cluster[seqname]].append(seqname)
    selectedseqs = []
    for cluster, seqnames in cluster2seq.items():
        logger.info("Cluster %s comprises %d sequences", cluster, len(seqnames))
        if len(seqnames) < min_clust_size:
            logger.debug("Cluster %s below size cutoff", cluster)
            continue
        if len(seqnames) > max_clust_size:
            logger.debug(
                "Cluster %s has %d reads, downsampling with random seed %d...",
                cluster,
                len(seqnames),
                rseed,
            )
            seed(rseed)
            selected_idx = sorted(sample(range(len(seqnames)), k=max_clust_size))
            selectedseqs.extend([seqnames[i] for i in selected_idx])
        else:
            selectedseqs.extend(seqnames)
        # Do not use a context manager here because we need file later
        fastq_handles[cluster] = NamedTemporaryFile(
            suffix=".fastq",
            mode="w",
            delete=(not keeptmp),
            delete_on_close=False,
        )
    # Write fastq only for clusters that are above min_clust_size and
    # downsampled if above max_clust_size
    for seqname, seq, qual in pyfastx.Fastx(reads):
        if seqname in selectedseqs and seq2cluster[seqname] in fastq_handles:
            fastq_rec = f"@{seqname!s}\n{seq!s}\n+\n{qual!s}\n"
            fastq_handles[seq2cluster[seqname]].write(fastq_rec)
    for handle in fastq_handles.values():
        handle.close()
    return fastq_handles, cluster2seq


def cluster_seqs_from_mcl(
    mcl_out: Path,
    reads: Path,
    keeptmp: bool,
    min_clust_size: int = 5,
    max_clust_size: int = 1000,
    rseed: int = 12345,
) -> tuple:
    """Extract sequences from each MCL cluster to Fastq files.

    :param mcl_out: Path to MCL output file
    :param reads: Path to reads from extract_reads_for_ava step
    :param keeptmp: Keep temporary files?
    :param min_clust_size: Minimum cluster size to return Fastq file
    :returns: Dict of file handles to each Fastq file keyed by cluster ID;
        excludes clusters below min_clust_size, and reads for clusters above
        max_clust_size are downsampled to max_clust_size.
    :returns: Dict of sequence IDs keyed by cluster ID; includes all sequences
        regardless of cluster size.
    :rtype: tuple
    """
    fastq_handles = {}
    seq2cluster = {}
    selectedseqs = []
    with Path.open(mcl_out) as fh:
        for counter, line in enumerate(fh):
            seqnames = line.rstrip().split("\t")
            logger.info("Cluster %d comprises %d sequences", counter, len(seqnames))
            for seqname in seqnames:
                seq2cluster[seqname] = counter  # assume each seq in only one cluster
            if len(seqnames) < min_clust_size:
                logger.debug("Cluster %d below size cutoff", counter)
                continue
            if len(seqnames) > max_clust_size:
                logger.debug(
                    "Cluster %d has %d reads, downsampling with random seed %d...",
                    counter,
                    len(seqnames),
                    rseed,
                )
                seed(rseed)
                selected_idx = sorted(sample(range(len(seqnames)), k=max_clust_size))
                selectedseqs.extend([seqnames[i] for i in selected_idx])
            else:
                selectedseqs.extend(seqnames)
            # Don't use context manager here because we need file later
            fastq_handles[counter] = NamedTemporaryFile(
                suffix=".fastq",
                mode="w",
                delete=(not keeptmp),
                delete_on_close=False,
            )
    # Write fastq only for clusters that are above min_clust_size and
    # downsampled if above max_clust_size
    for seqname, seq, qual in pyfastx.Fastx(reads):
        if seqname in selectedseqs and seq2cluster[seqname] in fastq_handles:
            fastq_rec = f"@{seqname!s}\n{seq!s}\n+\n{qual!s}\n"
            fastq_handles[seq2cluster[seqname]].write(fastq_rec)
    # Cluster memberships for all read segments regardless of cluster size
    cluster2seq = defaultdict(list)
    for seqname, cluster in seq2cluster.items():
        cluster2seq[cluster].append(seqname)
    for handle in fastq_handles.values():
        handle.close()
    return fastq_handles, cluster2seq


def spoa_assemble_fasta(label_fastq: tuple) -> tuple[str, str]:
    """Run spoa assembly on a Fastq input file.

    :param label_fastq: Tuple of input label (str) and path to Fastq file with
        reads to assemble
    :returns: Tuple of input label and alignment of consensus and input
        sequences in Fasta format (stdout from spoa -r 2)
    :rtype: tuple
    """
    label, fastq = label_fastq
    cmd = ["spoa", "-r", "2", fastq]
    logger.debug("spoa command: %s", " ".join([str(i) for i in cmd]))
    proc = Popen(cmd, stdout=PIPE, text=True)
    # TODO directing stderr to PIPE and logger prevents mp.Pool from closing?
    # ignoring stderr from spoa for now
    out = proc.communicate()[0]
    logger.info("Assembly complete for cluster %s", str(label))
    return label, out


def parse_spoa_r2(fasta: str) -> dict:
    """Parse spoa gapped alignment + consensus in Fasta format (-r 2 output).

    :param fasta: Assembly from spoa as list of strings
    :returns: dict of sequences (multiple lines concatenated) keyed by headers.
        The conensus sequence assembled by spoa must have key "Consensus".
    :rtype: dict
    """
    seqs = {}
    prev_hdr = ""
    prev_seq = ""
    for line in fasta.split("\n"):
        if line.startswith(">"):
            if len(prev_hdr) == 0 and len(prev_seq) == 0:
                prev_hdr = line.rstrip()[1:]
            elif len(prev_hdr) > 0 and len(prev_seq) > 0:
                seqs[prev_hdr] = prev_seq
                prev_seq = ""
                prev_hdr = line.rstrip()[1:]
        else:
            prev_seq += line.rstrip()
    # catch last
    if len(prev_hdr) > 0 and len(prev_seq) > 0:
        seqs[prev_hdr] = prev_seq
    return seqs


def count_spoa_aln_vars(seqs: dict) -> dict:
    """Count mismatches and gaps vs consensus for each sequence in a cluster.

    :param seqs: Dict of sequences parsed by parse_spoa_r2; the consensus
        sequence must have key "Consensus"
    :returns: Count of base matches and mismatches, gaps relative to query and
        consensus, leading and trailing query gaps, for each sequence in the
        dict relative to consensus
    :rtype: dict
    """
    var = defaultdict(lambda: defaultdict(int))
    cons = seqs["Consensus"]
    for hdr, seq in seqs.items():
        if hdr == "Consensus":
            continue
        # get coordinates without trailing and leading gaps
        span = re.search(r"^-*([^-].*[^-])-*$", seq).span(1)
        var[hdr]["query_lead_gap"] = span[0]
        var[hdr]["query_trail_gap"] = len(seq) - span[1]
        for i in range(span[0], span[1]):
            if seq[i] == cons[i] and seq[i] != "-":
                # Matching base but not common gap
                var[hdr]["match"] += 1
            elif seq[i] == "-" and cons[i] != "-":
                var[hdr]["query_gap"] += 1
            elif seq[i] != "-" and cons[i] == "-":
                var[hdr]["cons_gap"] += 1
            elif seq[i] != "-" and cons[i] != "-" and seq[i] != cons[i]:
                var[hdr]["mismatch"] += 1
    return var


def count_spoa_aln_persite_vars(seqs: dict) -> dict:
    """[WIP] Calculate entropy per alignment position for clustered sequences vs consensus.

    If mismatches/gaps between sequences and the consensus are solely due to
    technical sequencing error, rather than erroneous clustering of divergent
    underlying reads, we expect fraction of variants per site to be roughly the
    sequencing error. If for example two sequences are misclustered, then we
    should see a secondary peak of ~50% variant coverage.

    :param seqs: Dict of aligned sequences parsed by parse_spoa_r2; the
        consensus sequence must have key "Consensus"
    :returns: dict keyed by alignment position, value is the sequence entropy
        per site in bits.
    """
    var = {}
    for i in range(len(seqs["Consensus"])):
        column = [seqs[hdr][i] for hdr in seqs if hdr != "Consensus"]
        _values, counts = np.unique(column, return_counts=True)
        counts_norm = counts / counts.sum()
        var[i] = -(counts_norm * np.log(counts_norm) / np.log(2)).sum()
    return var


def check_stage_file(stage: str, message: str):
    """Check whether outputs for each stage of Pipeline already exist.

    Output files are checked by their expected filenames. If the expected
    outputs already exist, the stage is skipped entirely. Use this function
    as a decorator for individual stage methods. Enables resuming the
    pipeline from partial output.

    :param stage: Name of run stage, must be a key of self.OUTFILE_SUFFIX
    :param message: Logging message to emit on starting each stage
    """

    # TODO: check required input files
    def check_stage_decorator(func):
        @wraps(func)
        def wrapped_function(self, *args, **kwargs):
            if self.check_run_file(stage):
                if self._resume:
                    logger.info(
                        "Stage %s file output already present, skipping",
                        stage,
                    )
                    return None
                logger.error(
                    "Stage %s file output already present and option --resume not used, exiting",
                    stage,
                )
                sys.exit()
            else:
                logger.info(message)
                return func(self, *args, **kwargs)

        return wrapped_function

    return check_stage_decorator


class Pipeline:
    """Generic pipeline class to manage intermediate file paths.

    The Pipeline class is intended to be subclassed by specific pipelines,
    which will define run stages and expected output files in a dict
    self.OUTFILE_SUFFIX.
    """

    def check_run_file(self, stage: str) -> bool:
        """Check if intermediate output file has been created.

        :param stage: Name of run stage, must be a key of self.OUTFILE_SUFFIX
        :returns: True if file already exists at expected path
        :rtype: bool
        """
        return Path.is_file(self.pathto(stage))

    def pathto(self, stage: str, basename_only: bool = False) -> Path:
        """Combine output directory prefix and filenames to intermediate file path.

        :param stage: Name of run stage, must be a key of self.OUTFILE_SUFFIX
        :param basename_only: Only report the base filename if True
        :returns: Expected path to intermediate output file
        :rtype: Path
        """
        try:
            if basename_only:
                return Path(Path(self._prefix + self.OUTFILE_SUFFIX[stage]).name)
            return Path(self._outdir) / Path(self._prefix + self.OUTFILE_SUFFIX[stage])
        except KeyError as e:
            e.add_note(f"Unknown intermediate file {stage}")
            raise

    def db_taxonomy(self) -> None:
        """Get taxonomy string from SILVA headers in database Fasta file.

        Parse SILVA-style headers to get dict of taxonomy strings keyed by
        accession, stored in the Run._acc2tax attribute.
        """
        logger.info("Reading taxonomy from SILVA database file")
        self._acc2tax = {}
        with Path.open(self._ref) as fh:
            for line in fh:
                if line.startswith(">"):
                    spl = line.lstrip(">").rstrip().split(" ")
                    acc = spl[0]
                    taxstring = " ".join(spl[1:]).split(";")
                    self._acc2tax[acc] = taxstring
        logger.debug(" Accessions read: %d", len(self._acc2tax))

    def isonclust3_cluster(self, outfolder: str | Path, reads: str | Path) -> int:
        """Cluster marker read segments with isONclust3.

        isONclust3 was originally designed for long-read transcriptome
        datasets, and uses minimizers with iterative cluster merging. The
        default preset for Nanopore, "ont", is applied here for both legacy and
        Q20+ Nanopore reads, while the "pacbio" preset is used for both PacBio
        CLR and HiFi. Cluster output is a directory with multiple files.

        :param outfolder: Path to folder to write results
        :param reads: Path to reads in Fastq format
        :returns: Return code of the isONclust3 process.
        :rtype: int
        """
        if self._platform in ["lr:hq", "map-ont"]:
            mode = "ont"
        elif self._platform in ["map-hifi", "map-pb"]:
            mode = "pacbio"
        cmd = [
            "isONclust3",
            "--no-fastq",
            "--fastq",
            reads,
            "--mode",
            mode,
            "--outfolder",
            outfolder,
        ]
        if mode == "ont":
            cmd.append("--post-cluster")
        logger.debug("isonclust3 command: %s", " ".join([str(i) for i in cmd]))
        proc = Popen(cmd, stdout=PIPE)
        for l in proc.stdout:
            logger.debug("  isonclust3 log: %s", l.decode().rstrip())
        return proc.wait()

    def assemble_clusters(
        self,
        cluster_out: str | Path,
        reads: str | Path,
        cluster_asm: str | Path,
        cluster_tool: str = "isonclust3",
        threads: int = 12,
        rseed: int = 12345,
        keeptmp: bool = False,
        min_clust_size: int = 5,
        max_clust_size: int = 500,
    ) -> None:
        """Extract cluster sequences and assemble with spoa.

        Cluster reads with specified clustering tool, extract read segments for
        each cluster to separate Fastq files, and assemble consensus sequence
        for each with Spoa. Assembly step is embarrassingly parallelized.

        Updates Pipeline._stats["cluster2seq"] with lists of read names keyed
        by cluster id (number prefixed with `cluster_`). Some summary stats on
        clusters are written to Pipeline._stats["runstats"].

        Other per-cluster summaries in Pipeline._stats:
          * "cluster variant counts" -- number of variant sites by type for
            each read vs. its cluster consensus.
          * "cluster persite variant counts" -- per-site entropy vs. consensus
            [WIP].
          * "cluster cons parsed" -- cluster alignment including consensus,
            produced by spoa.

        :param cluster_out: Path to output file from clustering step; expected
            format depends on cluster_tool
        :param reads: Path to Fastq file with read segments to be clustered and
            assembled
        :param cluster_asm: Path to write assembled cluster consensus sequences
            in Fasta format
        :param cluster_tool: Clustering method used, either "mcl" or "isonclust3".
        :param threads: Number of parallel assembly jobs to run.
        :param rseed: Random seed for downsampling reads in large clusters.
        :param keeptmp: If True, do not delete Fastq files with extracted reads.
        :param min_clust_size: Only assemble clusters containing at least this
            number of reads.
        :param max_clust_size: Downsample reads for clusters above this size.
        """
        if cluster_tool == "mcl":
            fastq_handles, cluster2seq = cluster_seqs_from_mcl(
                cluster_out,
                reads,
                keeptmp=keeptmp,
                min_clust_size=min_clust_size,
                max_clust_size=max_clust_size,
                rseed=rseed,
            )
        elif cluster_tool == "isonclust3":
            fastq_handles, cluster2seq = cluster_seqs_from_isonclust3(
                cluster_out,
                reads,
                keeptmp=keeptmp,
                min_clust_size=min_clust_size,
                max_clust_size=max_clust_size,
                rseed=rseed,
            )
        self._stats.update({"cluster2seq": cluster2seq})
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
        with Path.open(cluster_asm, "w") as fh:
            fh.writelines(
                f">cluster_{c!s} Consensus\n"
                + cluster_cons_parsed[c]["Consensus"].replace("-", "")
                + "\n"
                for c in cluster_cons_parsed
            )
            logger.info("Assembled sequences written to %s", cluster_asm)

    def cluster_cons_mafft_aln(
        self,
        cluster_asm: str | Path,
        cluster_cons_aln: str | Path,
        threads: int = 12,
    ) -> int:
        """Multiple sequence alignment of cluster consensus sequences with MAFFT.

        Alignment will be used to generate a phylogenetic tree of assembled
        marker sequences to include in report. Use MAFFT `linsi` algorithm
        because this performs best with sequences with mixture of conserved and
        variable regions with variable lengths, as expected for rRNA gene.

        :param cluster_asm: Path to assembled cluster consensus sequences in Fasta format
        :param cluster_cons_aln: Path to write aligned cluster consensus sequences in Fasta format
        :param threads: Number of threads for MAFFT to use
        :returns: Return code of the MAFFT process
        """
        cmd = [
            "linsi",
            "--thread",
            str(threads),
            cluster_asm,
        ]
        logger.debug("MAFFT command: %s", " ".join([str(i) for i in cmd]))
        with Path.open(cluster_cons_aln, "w") as fh:
            proc = Popen(cmd, stdout=fh, stderr=PIPE, text=True)
            for l in proc.stderr:
                logger.debug("  MAFFT log: %s", l.rstrip())
            return proc.wait()

    def cluster_cons_distance_tree(
        self,
        cluster_cons_aln: str | Path,
        cluster_cons_tree: str | Path,
    ) -> None:
        """Generate distance tree from aligned cluster consensus sequences.

        Use Bio.Phylo.TreeConstruction to generate a distance tree from the
        aligned cluster consensus sequences. Simple tree construction method is
        adequate for visualization and comparison of sample composition.

        Bio.Phylo.Tree object is written to internal self._constree attribute
        and also to file in Newick format.

        :param cluster_cons_aln: Path to aligned cluster consensus sequences in Fasta format
        :param cluster_cons_tree: Path to write distance tree in Newick format
        """
        aln = AlignIO.read(cluster_cons_aln, "fasta")
        calculator = DistanceCalculator("identity")
        dm = calculator.get_distance(aln)
        constructor = DistanceTreeConstructor()
        self._constree = constructor.nj(dm)
        self._constree.root_at_midpoint()
        self._constree.ladderize()
        phylo_write(self._constree, cluster_cons_tree, "newick")
        logger.info("Distance tree written to %s", cluster_cons_tree)

    def cluster_asm_tophits(
        self,
        tophits: str | Path,
        asm: str | Path,
        threads: int = 12,
    ) -> int:
        """Map sequences to reference database with minimap2 and get top hits.

        Find best reference match for cluster consensus sequences by aligning
        with minimap2. Alignment is written to file.

        :param tophits: Path to write minimap2 alignment output in PAF format
        :param asm: Path to assembled cluster consensus sequences in Fasta format
        :param threads: Number of threads for minimap2 to use
        :returns: Return code of the minimap2 process
        """
        cmd = [
            "minimap2",
            "-x",
            "asm5",
            "-c",
            "--eqx",
            "--secondary=no",
            "-t",
            str(threads),
            "-o",
            tophits,
        ]
        (
            cmd.extend([self._refindex, asm])
            if self._refindex is not None
            else cmd.extend([self._ref, asm])
        )
        logger.debug("minimap command: %s", " ".join([str(i) for i in cmd]))
        proc = Popen(cmd, stderr=PIPE)
        for l in proc.stderr:
            logger.debug("  minimap log: %s", l.decode().rstrip())
        return proc.wait()

    def summarize_tophit_paf(self, tophits: str | Path) -> None:
        """Summarize top hits of sequences mapped to SILVA database by minimap2.

        Update Pipeline._stats with a dict of summary stats "cluster_tophits"
        for each hit, keyed by query sequence name.

        :param tophits: Path to minimap2 alignment output in PAF format
        """
        out = {}

        with Path.open(tophits) as fh:
            for line in fh:
                spl = line.rstrip().split("\t")
                hits = dict(
                    zip(
                        [
                            "qname",
                            "qlen",
                            "qstart",
                            "qend",
                            "strand",
                            "tname",
                            "tlen",
                            "tstart",
                            "tend",
                            "alnmatch",
                            "alnlen",
                        ],
                        spl[0:11],
                        strict=False,
                    ),
                )
                cigar = next(i for i in spl if i.startswith("cg:Z:"))
                # Calculate derived metrics from PAF fields
                cigar_summary = parse_cigar_ops(cigar[5:])
                hits.update({CIGAROPS[c]: cigar_summary[c] for c in cigar_summary})
                # remove redundant % sign for display
                hits["align %id"] = "{:.2%}".format(
                    int(hits["alnmatch"]) / int(hits["alnlen"]),
                ).rstrip("%")
                hits["query %aln"] = "{:.2%}".format(
                    (int(hits["qend"]) - int(hits["qstart"])) / int(hits["qlen"]),
                ).rstrip("%")
                hits["target %aln"] = "{:.2%}".format(
                    (int(hits["tend"]) - int(hits["tstart"])) / int(hits["tlen"]),
                ).rstrip("%")
                out[spl[0]] = hits

        # Taxonomy of hit targets
        for rec in out.values():
            try:
                # hyperlink to ENA record
                # TODO: This only works for SILVA where accessions are derived
                # from ENA accessions. May not work with other databases e.g.
                # Greengenes
                rec["tophit"] = (
                    f"[{rec['tname']!s}](https://www.ebi.ac.uk/ena/browser/view/{rec['tname'].split('.')[0]!s})"
                )
                rec["tophit taxonomy"] = ";".join(self._acc2tax[rec["tname"]])
                rec["tophit species"] = self._acc2tax[rec["tname"]][-1]
                # Higher taxonomy to class level, except for chloroplast and mitochondria
                # Assumes SILVA taxonomy is in use
                if rec["tophit taxonomy"].startswith(
                    "Bacteria;Cyanobacteria;Cyanobacteriia;Chloroplast",
                ):
                    rec["higher taxonomy"] = "[Eukaryotic organelle Chloroplast]"
                elif rec["tophit taxonomy"].startswith(
                    "Bacteria;Proteobacteria;Alphaproteobacteria;Rickettsiales;Mitochondria",
                ):
                    rec["higher taxonomy"] = "[Eukaryotic organelle Mitochondria]"
                else:
                    rec["higher taxonomy"] = "; ".join(
                        self._acc2tax[rec["tname"]][0:-1],
                    )
            except KeyError as e:
                e.add_note(f"Accession {rec['tname']} not found in database?")
                raise
        self._stats.update({"cluster_tophits": out})

    def write_report_json(self, out: str | Path) -> None:
        """Dump run stats file in JSON format.

        Dump the Pipeline._stats attribute to a JSON file for troubleshooting
        and comparison of different phyloblitz runs.

        :param out: Path to write JSON file
        """
        self._stats.update({"endtime": str(datetime.now())})
        with Path.open(out, "w") as fh:
            json.dump(self._stats, fh, indent=4)
