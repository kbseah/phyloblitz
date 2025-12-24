# phyloblitz 🌳🌩️

**Rapid SSU rRNA marker gene screening of long read metagenomes.**

Note: This tool is still under development and the interface may change without
warning. We appreciate any bug reports or feedback; please [create an issue on
GitHub](https://github.com/kbseah/phyloblitz/issues/new). Please understand
that, because of limited resources, we cannot guarantee a quick response.

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

Dependencies are managed with [pixi](https://pixi.sh/). When `pixi shell` is
run for the first time in the folder containing the `pixi.toml` configuration
file, pixi will resolve and install dependencies. Run `phyloblitz` without
arguments or with the `--help` parameter to view help message.

```console
git clone git@github.com:kbseah/phyloblitz.git
cd phyloblitz
pixi shell # set up workspace and start pixi shell session
phyloblitz --help
```

Use `exit` or `Ctrl-D` to exit pixi shell session.

Required inputs are a preprocessed SILVA reference database and the long read
files in fastq(.gz) format.


## Reference database

Download the database files from https://doi.org/10.5281/zenodo.7892522

These have been formatted and indexed for use with
[phyloFlash](https://hrgv.github.io/phyloFlash/). For `phyloblitz`, only the
Fasta file with trimmed, dereplicated sequences is necessary:
`SILVA_SSU.noLSU.masked.trimmed.NR96.fixed.fasta`.

You can use `pixi run download` to execute the download command.

If you build your own database, low complexity regions should be masked (e.g.
with `bbmask.sh`) otherwise mapping will be slow. Refer to the original
phyloFlash paper ([Gruber-Vodicka, Seah & Pruesse,
2021](https://doi.org/10.1128/mSystems.00920-20)) for details.


### Test datasets

After downloading the database file with `pixi run download` (or manually to
the subfolder `138.1/`) you can test `phyloblitz` on the test data (a sample of
SSU rRNA-containing Nanopore reads from SRR17913200) with `pixi run run`. The
output will be in `runtest/run`.

Try out `phyloblitz` on other published datasets that have sequenced the
[ZymoBIOMICS Gut Microbiome
Standard](https://www.zymoresearch.com/products/zymobiomics-gut-microbiome-standard):

* Nanopore GridION R9.4.1 flowcells, Kit 9 and Q20+ chemistry: SRR23926923
  ([Liu et al., 2022](https://doi.org/10.1186/s40168-022-01415-8))
* PacBio Sequel II, Standard Input library: SRR13128014

A subsample of ~100k reads (`--num_reads 100000`) from a metagenome dataset is
sufficient for a quick overview.


## Pipeline overview

 * Map reads to SILVA with minimap2 and retain only mapped reads
 * Summarize taxonomy across mapped reads, using consensus taxonomy of each read
 * Extract aligned portion of mapped reads, use primary alignment only
     otherwise we get duplicate read segments with slightly differing
     boundaries
 * Sequence clustering, either with isONclust3 or minimap2 + mcl
 * Assemble consensus sequence per read cluster with spoa
 * Generate metrics per cluster for diagnostics: expect similar error rate per
     read vs. consensus, identify clusters with too few reads (TODO)


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

v0.2.0 targets:

 - [x] Check if lr:hq mode is better for Q20 ONT reads than map:ont
 - [x] Option to cluster with isONclust3 instead of minimap2 ava
 - [ ] Option to search final assembled cluster sequences with vsearch instead
       of minimap2; better handling of divergent sequences: Example, the
       Aestuariibacter sequence in SRR28830816
 - [ ] Divergence of reads in each cluster from consensus to detect potential
       chimeras
 - [ ] Sequence diversity in reads flanking the SSU sequence per cluster to
       identify meta-consensus sequences
 - [ ] Detailed documentation
 - [ ] Multisample comparison and co-assembly

Future plans:

 - [ ] Does extracting flanking sequence context improve strain resolution?
     - [ ] Filter out all-vs-all hits with overhangs
     - [ ] Investigate effect of dv cutoff values and clustering methods
 - [ ] Embed graphics as PNGs into HTML with markdown-embedimages
 - [ ] Benchmarking against defined test datasets
 - [ ] MultiQC integration
 - [ ] Fix the CSS stylesheet


## Citations

Please cite the following dependencies:

 * [`minimap2`](https://github.com/lh3/minimap2) Heng Li, Minimap2: pairwise
     alignment for nucleotide sequences, Bioinformatics, 34(18):3094-3100,
     2018.  https://doi.org/10.1093/bioinformatics/bty191
 * [`isONclust3`](https://github.com/aljpetri/isONclust3) Alexander J Petri,
     Kristoffer Sahlin, De novo clustering of large long-read transcriptome
     datasets with isONclust3, Bioinformatics, 41(5):batf207, 2025.
     https://doi.org/10.1093/bioinformatics/btaf207
 * [`pymarkovclustering`](https://github.com/moshi4/pyMarkovClustering),
     implementation of MCL in Python
 * [`mcl`](https://micans.org/mcl/) Stijn van Dongen, Graph Clustering Via a
     Discrete Uncoupling Process, SIAM Journal on Matrix Analysis and
     Applications, 30(1):121-141, 2008.
     http://link.aip.org/link/?SJMAEL/30/121/1
 * [`pyfastx`](https://pyfastx.readthedocs.io/) Lianming Du, et al., Pyfastx: a
     robust Python package for fast random access to sequences from plain and
     gzipped FASTA/Q files, Briefings in Bioinformatics, 22(4):bbaa368, 2020.
     https://doi.org/10.1093/bib/bbaa368
 * [`samtools`](https://www.htslib.org/) Heng Li, Bob Handsaker, et al., The
     Sequence Alignment/Map format and SAMtools, Bioinformatics,
     25(16):2078-2079, 2009. https://doi.org/10.1093/bioinformatics/btp352
     James K Bonfield, John Marshall, Petr Danecek, et al., HTSlib: C library
     for reading/writing high-throughput sequencing data, GigaScience,
     10(2):giab007, 2021. https://doi.org/10.1093/gigascience/giab007
 * [`pysam`](https://github.com/pysam-developers/pysam)
 * [`spoa`](https://github.com/rvaser/spoa)

Please cite the SILVA reference database if you use it:

 * [SILVA](https://www.arb-silva.de/) Maria Chuvochina, Jan Gerken, et al.
     SILVA in 2026: a global core biodata resource for rRNA within the DSMZ
     digital diversity. Nucleic Acids Research, gkaf1247, 2026.
     https://doi.org/10.1093/nar/gkaf1247

If you use `phyloblitz` in published research, please cite the GitHub
repository URL and software version.


## License

`phyloblitz` is distributed under the terms of the
[MIT](https://spdx.org/licenses/MIT.html) license.
