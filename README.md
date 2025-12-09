# phyloblitz 🌳🌩️

Rapid SSU rRNA marker gene screening of long read metagenomes

-----

`phyloblitz` is a tool to rapidly screen long read metagenomes with ribosomal
RNA genes, generating a taxonomic summary and targeted assemblies of
full-length SSU rRNA sequences.

 * Quickly check the composition and complexity of metagenomes before
     committing time and resources to de novo assembly.
 * Leverage the broad phylogenetic coverage of the SILVA rRNA database to
     detect taxa underrepresented in genome databases.
 * Use assembled gene sequences for phylogenetics and probe design.

The rRNA-targeted screening approach for metagenomes was originally implemented
for Illumina reads in the software tool
[phyloFlash](https://hrgv.github.io/phyloFlash/); the name `phyloblitz`
references phyloFlash and [bioblitzes](https://en.wikipedia.org/wiki/BioBlitz).


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
[phyloFlash](https://hrgv.github.io/phyloFlash/). For `phyloblitz`, only the
Fasta file with trimmed, dereplicated sequences is necessary:
`SILVA_SSU.noLSU.masked.trimmed.NR99.fixed.fasta`.

If you build your own database, low complexity regions should be masked (e.g.
with `bbmask.sh`) otherwise mapping will be slow. Refer to the original
phyloFlash paper ([Gruber-Vodicka, Seah & Pruesse,
2021](https://doi.org/10.1128/mSystems.00920-20)) for details.


### Test datasets

To try out `phyloblitz`, you can use published datasets of the [ZymoBIOMICS Gut
Microbiome
Standard](https://www.zymoresearch.com/products/zymobiomics-gut-microbiome-standard):

* Nanopore GridION R9.4.1 flowcells, Kit 9 and Q20+ chemistry: SRR23926923
  ([Liu et al., 2022](https://doi.org/10.1186/s40168-022-01415-8))
* PacBio Sequel II, Standard Input library: SRR13128014

A subsample of ~100k reads is sufficient for a quick overview.


## Pipeline overview

 * Map reads to SILVA with minimap2 and retain only mapped reads
 * Summarize taxonomy across mapped reads, using consensus taxonomy of each read
 * Extract aligned portion of mapped reads, use primary alignment only
     otherwise we get duplicate read segments with slightly differing
     boundaries
 * All-vs-all mapping of extracted read segments with minimap2
 * Filter out all-vs-all mappings where per-base sequence divergence reported
     by minimap (`dv:` tag) is too low; this parameter should be adjusted,
     depending on the expected read accuracy of the sequencing platform
 * Cluster sequences with mcl
 * Assemble consensus sequence per read cluster with spoa
 * Generate metrics per cluster for diagnostics: expect similar error rate per
     read vs. consensus, identify clusters with too few reads


## TODO

v0.1.0 targets:

 - [x] Taxonomy summary from initial mapping
 - [x] Reads in input, number mapped, number used for cluster assembly
 - [x] Choose dv_max threshold from observed values in all-vs-all mapping
 - [x] Object orientation 🫠
 - [x] CSS stylesheet for HTML report
 - [x] Check if tool works on PacBio data
 - [x] Option to map only a sample of reads
 - [x] Fix dv_max at a minimal value if observed median is zero
 - [x] Report metrics into a user-friendly file like phyloFlash
     - [x] Numbers of reads mapped and taxonomic summary
     - [x] Assembled sequences, top hits, and respective read coverage and
         cluster metrics (plus heuristic assembly quality score)
     - [x] Link SILVA record by accession
     - [x] Alignment identity (using Blast-like % id definition)

Future steps:

 - [ ] Embed graphics as PNGs into HTML with markdown-embedimages
 - [ ] Divergence of reads in each cluster from consensus to detect potential chimeras
 - [ ] Does extracting flanking sequence context improve strain resolution?
 - [ ] Investigate effect of dv cutoff values and clustering methods
 - [ ] Benchmarking against defined test datasets
 - [ ] Check if lr:hq mode is better for Q20 ONT reads than map:ont
 - [ ] MultiQC integration
 - [ ] Detailed documentation


## Citations

Please cite the following dependencies:

 * [`spoa`](https://github.com/rvaser/spoa)
 * [`mcl`](https://micans.org/mcl/) Stijn van Dongen, Graph Clustering Via a
     Discrete Uncoupling Process, SIAM Journal on Matrix Analysis and
     Applications, 30(1):121-141, 2008.
     http://link.aip.org/link/?SJMAEL/30/121/1
 * [`minimap2`](https://github.com/lh3/minimap2) Heng Li, Minimap2: pairwise
     alignment for nucleotide sequences, Bioinformatics, 34(18):3094-3100,
     2018.  https://doi.org/10.1093/bioinformatics/bty191
 * [`pyfastx`](https://pyfastx.readthedocs.io/)
 * [`samtools`](https://www.htslib.org/)
 * [`pysam`](https://github.com/pysam-developers/pysam)

Please cite the SILVA reference database if you use it:

 * [SILVA](https://www.arb-silva.de/) Maria Chuvochina, Jan Gerken, et al.
     SILVA in 2026: a global core biodata resource for rRNA within the DSMZ
     digital diversity. Nucleic Acids Research, gkaf1247, 2026.
     https://doi.org/10.1093/nar/gkaf1247


## License

`phyloblitz` is distributed under the terms of the
[MIT](https://spdx.org/licenses/MIT.html) license.
