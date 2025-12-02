# phyloblitz

-----

The name is a reference to precursor tool
[phyloFlash](https://hrgv.github.io/phyloFlash/) and to
[bioblitz](https://en.wikipedia.org/wiki/BioBlitz) events.

## Usage

Dependencies are managed with [pixi](https://pixi.sh/). Run without arguments
to view help message.

```console
pixi run phyloblitz
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
