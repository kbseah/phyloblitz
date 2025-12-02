# phyloblitz 🌳🌩️

-----

The name is a reference to precursor tool
[phyloFlash](https://hrgv.github.io/phyloFlash/) and to
[bioblitz](https://en.wikipedia.org/wiki/BioBlitz) events.

## Usage

Dependencies are managed with [pixi](https://pixi.sh/). Run without arguments
or with the `--help` parameter to view help message. When run for the first
time, pixi will resolve and install dependencies.

```console
pixi run phyloblitz --help
```

Required inputs are a preprocessed SILVA reference database and the long read
files in fastq(.gz) format.

## Reference database

Download the database files from https://doi.org/10.5281/zenodo.7892522

These have been formatted and indexed for use with
[phyloFlash](https://hrgv.github.io/phyloFlash/). However for `phyloblitz`,
only the Fasta file with trimmed, dereplicated sequences is necessary:
`SILVA_SSU.noLSU.masked.trimmed.NR99.fixed.fasta`.

## License

`phyloblitz` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.

## Pipeline

 * Map to SILVA with minimap2 and output SAM
 * Taxonomic summary per read with primary + secondary alignments
 * Extract aligned portion of mapped reads (optionally with flanking
   sequence?), use primary alignment only otherwise we get duplicate read
   segments with slightly differing boundaries
 * All-vs-all mapping of extracted reads with minimap2
 * Define read clusters: Filter out all-vs-all mappings where per-base sequence
   divergence reported by minimap (dv: tag) is less than expected sequence
   divergence from sequencing platform; use mcl to generate clusters (may not
   be necessary)
 * Get consensus sequence per read cluster with spoa
 * Generate metrics per consensus cluster for diagnostics: expect similar error
   rate per read vs. consensus, identify clusters with too few reads


## TODO

 * Investigate effect of dv cutoff values and clustering methods
 * Can we set thresholds more naturally by bootstrapping read metrics from mapping steps?

## Citations

Please cite the following dependencies:

 * [`spoa`](https://github.com/rvaser/spoa)
 * [`mcl`](https://micans.org/mcl/) Stijn van Dongen, Graph Clustering Via a Discrete Uncoupling Process, SIAM Journal on Matrix Analysis and Applications, 30(1):121-141, 2008. http://link.aip.org/link/?SJMAEL/30/121/1
 * [`minimap2`](https://github.com/lh3/minimap2) Heng Li, Minimap2: pairwise alignment for nucleotide sequences, Bioinformatics, 34(18):3094-3100, 2018. https://doi.org/10.1093/bioinformatics/bty191
 * [`pyfastx`](https://pyfastx.readthedocs.io/)
 * [`samtools`](https://www.htslib.org/)
 * [`pysam`](https://github.com/pysam-developers/pysam)

Please cite the SILVA reference database if you use it:

 * [SILVA](https://www.arb-silva.de/) Maria Chuvochina, Jan Gerken, et al. SILVA in 2026: a global core biodata resource for rRNA within the DSMZ digital diversity. Nucleic Acids Research, gkaf1247, 2026. https://doi.org/10.1093/nar/gkaf1247
