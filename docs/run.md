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
 * `min_clust_size` - Minimum cluster size to assemble a consensus sequence. If
    there are too few reads, the cluster consensus will be fragmented and
    contain more errors.
 * `max_clust_size` - Clusters above this size will be downsampled for
    consensus assembly. Beyond about 500 reads, assembly becomes noticeably
    slower without much improvement in quality.
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
