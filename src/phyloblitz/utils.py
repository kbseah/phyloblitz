"""Shared utility functions for phyloblitz."""

import logging
import pyfastx
import re
from collections import defaultdict
from hashlib import md5
from os import W_OK, access
from pathlib import Path
from subprocess import PIPE, STDOUT, Popen
from sys import version as python_version
from matplotlib import __version__ as matplotlib_version
from mistune import __version__ as mistune_version
from numpy import __version__ as numpy_version
from pyfastx import __version__ as pyfastx_version
from pymarkovclustering import __version__ as pymcl_version
from pysam import __version__ as pysam_version
from tempfile import NamedTemporaryFile

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
    """Summarize operations in a CIGAR string."""
    summary = defaultdict(int)
    cigar = re.findall(r"(\d+)(\D)", cigar)
    for num, op in cigar:
        summary[op] += int(num)
    return summary


def lists_common_prefix(lol):
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


def filter_paf_overhang(line: str, max_overhang_frac: float = 0.05):
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
        ),
    )
    for tool in ["minimap2", "spoa", "isONclust3"]:
        p = Popen([tool, "--version"], stdout=PIPE, stderr=STDOUT, text=True)
        # Split in case version message is multiline, e.g. mcl
        vers[tool] = p.communicate()[0].split("\n")[0].rstrip()
    return vers


def check_outdir(outdir, resume=True):
    """Check if output directory exists, and create it if it doesn't

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
    isonclust3_out: Path, reads: Path, keeptmp: bool, min_clust_size: int = 5
) -> tuple:
    """Extract sequences from each isONclust3 cluster to Fastq files.

    The numbering of clusters by isONclust3 itself is not consistent
    between the final_clusters.tsv file and the output fastq files; best to
    extract the sequences ourselves.

    :param isonclust3_out: Path to isONclust3 output file
    :param reads: Path to reads from extract_reads_for_ava step
    :param keeptmp: Keep temporary files?
    :param min_clust_size: Minimum cluster size to return Fastq file
    :returns: Dict of file handles to each Fastq file keyed by cluster ID
    :returns: Dict of sequence IDs keyed by cluster ID
    """
    fastq_handles = {}
    seq2cluster = {}
    with Path.open(isonclust3_out) as fh:
        for line in fh:
            clust, seqname = line.rstrip().split("\t")
            seq2cluster[seqname] = clust
    cluster2seq = defaultdict(list)
    for seqname in seq2cluster:
        cluster2seq[seq2cluster[seqname]].append(seqname)
    for cluster in cluster2seq:
        logger.info(
            "Cluster %s comprises %d sequences", cluster, len(cluster2seq[cluster])
        )
        if len(cluster2seq[cluster]) >= min_clust_size:
            fastq_handles[cluster] = NamedTemporaryFile(
                suffix=".fastq",
                mode="w",
                delete=(not keeptmp),
                delete_on_close=False,
            )
        else:
            logger.debug("Cluster %s below size cutoff", cluster)
    for seqname, seq, qual in pyfastx.Fastx(reads):
        if seqname in seq2cluster and seq2cluster[seqname] in fastq_handles:
            fastq_rec = "@" + seqname + "\n" + seq + "\n+\n" + qual + "\n"
            fastq_handles[seq2cluster[seqname]].write(fastq_rec)
    for handle in fastq_handles.values():
        handle.close()
    return fastq_handles, cluster2seq


def cluster_seqs_from_mcl(
    mcl_out: Path, reads: Path, keeptmp: bool, min_clust_size: int = 5
) -> tuple:
    """Extract sequences from each MCL cluster to Fastq files.

    :param mcl_out: Path to MCL output file
    :param reads: Path to reads from extract_reads_for_ava step
    :param keeptmp: Keep temporary files?
    :param min_clust_size: Minimum cluster size to return Fastq file
    :returns: Dict of file handles to each Fastq file keyed by cluster ID
    :returns: Dict of sequence IDs keyed by cluster ID
    """
    counter = 0
    fastq_handles = {}
    seq2cluster = {}
    with Path.open(mcl_out) as fh:
        for line in fh:
            seqs = line.rstrip().split("\t")
            logger.info("Cluster %d comprises %d sequences", counter, len(seqs))
            for seq in seqs:
                seq2cluster[seq] = counter  # assume each seq in only one cluster
            if len(seqs) >= min_clust_size:
                fastq_handles[counter] = NamedTemporaryFile(
                    suffix=".fastq",
                    mode="w",
                    delete=(not keeptmp),
                    delete_on_close=False,
                )
            else:
                logger.debug("Cluster %d below size cutoff", counter)
            counter += 1
    for name, seq, qual in pyfastx.Fastx(reads):
        if name in seq2cluster and seq2cluster[name] in fastq_handles:
            fastq_rec = "@" + name + "\n" + seq + "\n+\n" + qual + "\n"
            fastq_handles[seq2cluster[name]].write(fastq_rec)
    cluster2seq = defaultdict(list)
    for seq, cluster in seq2cluster.items():
        cluster2seq[cluster].append(seq)
    for handle in fastq_handles.values():
        handle.close()
    return fastq_handles, cluster2seq
