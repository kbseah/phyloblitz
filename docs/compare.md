# Compare pipeline

The `phyloblitz compare` subcommand compares multiple datasets using their
individual `phyloblitz run` results files. The extracted marker sequence
segments from all samples are pooled, clustered, and co-assembled into
consensus sequences. The relative abundances of the clusters in each sample
(cluster memberships) are then used to cluster the samples.


## Command line options


### Input

 * `--db` - Path to preprocessed SILVA database fasta file
 * `--dbindex` - Path to minimap2 index of database file (optional)
 * `--input_table` - Path to table of sample names and report JSON file paths
    to compare, in TSV format with columns 'sample' and 'report'

### Output

 * `--outdir` - Output folder path
 * `--prefix` - Output filename prefix
 * `--log` - Path to write log file

### Run parameters

 * `--ignore_db_mismatch` - Continue even if different databases were used to
    generate reports
 * `--threads` - Number of parallel threads
 * `--min_clust_size` - Minimum cluster size to assemble a consensus sequence
 * `--max_clust_size` - Clusters above this size will be downsampled for
    consensus assembly
 * `--rseed` - Random seed for subsampling reads and downsampling clusters
 * `--cluster_method` - Method to use for comparing clusters between samples,
    passed to scipy.cluster.hierarchy.linkage method parameter; note that
    'ward' can only be used with euclidean distance. Choice of "ward",
    "single", "complete", "average", "weighted", "centroid", "median"
 * `--cluster_metric` - Distance metric to use for comparing clusters between
    samples. Choice of: "euclidean", "minkowski", "cityblock", "cosine",
    "jaccard"
 * `--debug` - Display logging DEBUG level messages to console


## Compare steps


### Check input files

Three checks are performed:

(1) Check that the same reference database was used for all the phyloblitz runs
being compared, and that it is the same as the one specified now to `--db`. The
reads extracted depend on the underlying database, so a comparison of runs
where different databases were used would not be fair.

(2) Check that read names are not repeated across samples.

(3) Check that the same sequencing platform was used for all samples. Different
platforms have substantially different error profiles so the clustering and
assembly of reads pooled from them may not be reliable.


### Pooled clustering of read segments and assembly of consensus

Read segments containing the marker of interest are extracted in individual
`phyloblitz` runs and stored in the JSON report file. Sequences are pooled,
clustered with `isonclust3`, and assembled with `spoa` with the same parameters
as for `phyloblitz run`.


### Count sequence cluster memberships per sample, cluster samples

Cluster consensus sequences are mapped to the reference database with
`minimap2` to find the top hits as in the `run` pipeline. The numbers of read
segments from each sample underlying each sequence cluster (cluster
memberships) are used to compare samples with each other, to draw a heatmap and
hierarchical clustering dendrogram of the samples.
