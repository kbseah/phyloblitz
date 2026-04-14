"""Shared utility functions for phyloblitz."""

import logging
import re
import sys
from collections import defaultdict
from functools import wraps
from hashlib import md5
from os import W_OK, access
from pathlib import Path
from random import sample, seed
from subprocess import PIPE, STDOUT, Popen
from sys import version as python_version
from tempfile import NamedTemporaryFile

import numpy as np
import pyfastx
from matplotlib import __version__ as matplotlib_version
from mistune import __version__ as mistune_version
from numpy import __version__ as numpy_version
from pyfastx import __version__ as pyfastx_version
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
        pass
    else:
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
                python_version,
                pysam_version,
                pyfastx_version,
                mistune_version,
                numpy_version,
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


def run_isonclust3(reads: str | Path, mode: str, outfolder: str | Path) -> int:
    """Run isONclust3.

    :param reads: Path to reads in Fastq format
    :param mode: Either "ont" or "pacbio"
    :param outfolder: Path to folder to write results
    :returns: Return code of the isONclust3 process.
    :rtype: int
    """
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
            fastq_rec = "@" + seqname + "\n" + seq + "\n+\n" + qual + "\n"
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
            fastq_rec = "@" + seqname + "\n" + seq + "\n+\n" + qual + "\n"
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

    # TODO specify more than one output file; check required input files
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
            e.add_not(f"Unknown intermediate file {stage}")
            raise
