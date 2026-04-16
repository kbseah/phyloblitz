# Run pipeline


## Options


### Input
 * `db` - Path to preprocessed SILVA database fasta file
 * `dbindex` - Path to minimap2 index of database file (optional, not
    provided). This may have a slight speed advantage by avoiding repeated
    indexing of the Fasta file.
 * `reads` - Read file in FASTQ or gzipped FASTQ format to screen.
 * `num_reads` - Use subsample of N reads; default is to use all reads.
 * `platform` - Sequencing platform used, argument passed to `minimap2 -x`
    option; one of "map-ont", "map-pb", "lr:hq", "map-hifi".


### Output
 * `prefix` - Output file name prefix
 * `outdir` - Output folder path. Will be created if it does not exist. If it
    does, the pipeline will stop unless `--resume` is specified, in which case
    it assumes this was an interrupted run and attempts to continue (at least
    the mapping step must be complete.)
 * `report` - Generate a report file; True by default, use `--noreport` to turn
    off reporting.
 * `keeptmp` - Do not delete temporary files.
 * `log` - Write logging messages to this file.
 * `write_cluster_alns` - Write cluster alignments to files.


### Run parameters
 * `threads` - Number of parallel threads; for mapping and cluster assembly.
 * `cluster_tool` - Tool to use for sequence clustering. Default is
    `isonclust3`. The other option using `mcl` is experimental.
 * `align_minlen` - Minimum length of aligned segment.
 * `summary_taxlevel` - Depth of taxonomy string for summary in report. For
    example, `4` (default) will summarize at taxonomic order.
 * `min_clust_size` - Minimum cluster size to assemble a consensus sequence.
 * `max_clust_size` - Clusters above this size will be downsampled for
    consensus assembly.
 * `rseed` - Random seed for subsampling reads and downsampling clusters.
 * `flanking` - Sequence flanking the mapped hits on query reads to extract, on
    both sides. Default is 1000 bp. If too short, the flanking sequence will be
    less informative. If too long, probability of finding reads of that length
    goes down substantially.
 * `no_supplementary` - Ignore supplementary alignments, only parse one marker
    sequence per read. May be useful if the genomes of interest contain high
    copy-number tandem rRNA gene arrays.
 * `resume` - Resume partially completed run based on expected filenames.
 * `debug` - Display logging DEBUG level messages to console. (Only affects
    console, all logging levels will be written to log file if `--log` is
    specified.)


### Options for MCL clustering

Only used if `--cluster_tool mcl` is specified. These are experimental for now.

 * `dv_max` - Maximum pairwise sequence divergence in minimap2 all-vs-all
    mapping to retain for clustering.
 * `dv_max_auto` - Set dv_max parameter automatically at max(0.001, 2 * median
    of all-vs-all divergence value).
 * `inflation` - Inflation parameter for MCL.


## Run steps


### Initial mapping

Input reads are mapped to the reference database with `minimap2`. Secondary
mappings and unmapped reads are not retained. Supplementary mappings are
retained because they correspond to multiple marker sequences present on the
same read.

The taxonomy of mapped read segments is summarized for all reads at the taxon
level specified with option `--summary_taxlevel`. This assumes that the
reference database is using an internally consistent set of taxon ranks in its
taxon strings, like SILVA. If a read has more than one marker sequence, the
consensus for that read is first taken, so each read is only counted once.


### Extract aligned segments and flanking sequence

The aligned segments of mapped reads are extracted (and reverse complemented if
necessary). If an aligned segment overlaps with another segment that has
already been extracted, only the first (primary) alignment is retained. This
sometimes happens with supplementary alignments to other reference sequences
that only partially overlap with the primary alignment. To avoid
double-counting, such overlapping segments should be skipped.

Flanking sequence context is also extracted, to evaluate potential
strain diversity. For example, if `--flanking 1000` is specified, 1 kbp
upstream and 1 kbp downstream of each aligned segment is extracted. For
downstream analyses, the aligned marker segment and the flanking sequences are
processed separately.


### Estimate sequencing error from all-vs-all mapping

Aligned segments are aligned against each other with `minimap2` in all-vs-all
("ava") mode. For Nanopore and PacBio CLR reads, the `minimap2` presets
`ava-ont` and `ava-pb` respectively are applied. For Nanopore high quality Q20+
reads and PacBio HiFi (CCS) reads, the "ava" presets were modified with
settings for those sequencing modes that reflect the lower sequence error.
Filter "ava" alignments to remove those that are not end-to-end (allow maximum
5% soft-clipped unaligned overhangs).


### Cluster marker segments

Cluster extracted segments with `isonclust3`. `--mode ont` setting used for
both legacy and Q20+ Nanopore reads, while `--mode pacbio` for both CLR and
HiFi PacBio reads. The `--post-cluster` option is applied for Nanopore reads.


### Assemble cluster consensus sequences

Clusters comprising total number of reads above the cutoff specified by
`--min_clust_size` are assembled with `spoa` to a consensus sequence for each
cluster. If there are too few reads, the cluster consensus will be fragmented
and contain more errors. If the cluster size is above `--max_clust_size`, the
reads are downsampled to `--max_clust_size` before assembly. Beyond about 500
reads, assembly becomes noticeably slower without much improvement in quality.


### Cluster flanking sequences of each cluster


### Map cluster consensus sequences to reference database
