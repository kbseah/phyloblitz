# phyloblitz 🌳🌩️

**Rapid SSU rRNA marker gene screening of long read metagenomes.**

![GitHub Tag](https://img.shields.io/github/v/tag/kbseah/phyloblitz)
![GitHub License](https://img.shields.io/github/license/kbseah/phyloblitz)
![Conda Version](https://img.shields.io/conda/vn/bioconda/phyloblitz)


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


## Installation


### Development version

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
files in fastq(.gz) format. After navigating to the phyloblitz folder, fetch
the (SSU rRNA) database file from Zenodo on the command line with:

```console
pixi run download
```

After downloading the database file to the subfolder `138.2/` you can test
`phyloblitz` on the test data (a sample of SSU rRNA-containing Nanopore reads
from SRR17913200) with `pixi run run`. The output will be in `runtest/run`.


### Reference database

Download the database files from https://doi.org/10.5281/zenodo.18627380

These are derived from SILVA database files, with some additional filtering and
trimming of contaminant sequences and masking of low-complexity repeats to
avoid excessive false-positive matches. Use the appropriate file for either SSU
or LSU rRNA genes. Refer to the [database processing
pipeline](https://github.com/kbseah/phyloblitz-db) for the actual commands
used to prepare the database.

If you build your own database, low complexity regions should be masked (e.g.
with `bbmask.sh`) otherwise mapping will be slow. Refer to the original
phyloFlash paper ([Gruber-Vodicka, Seah & Pruesse,
2021](https://doi.org/10.1128/mSystems.00920-20)) for details.


## Usage

Refer to the help message for details on the phyloblitz run parameters:

```console
phyloblitz --help
```


### Test datasets

Try out `phyloblitz` on published datasets that have sequenced the
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


## Citations

Please cite the following dependencies:

 * [`minimap2`](https://github.com/lh3/minimap2) Heng Li, Minimap2: pairwise
     alignment for nucleotide sequences, Bioinformatics, 34(18):3094-3100,
     2018.  https://doi.org/10.1093/bioinformatics/bty191
 * [`isONclust3`](https://github.com/aljpetri/isONclust3) Alexander J Petri,
     Kristoffer Sahlin, De novo clustering of large long-read transcriptome
     datasets with isONclust3, Bioinformatics, 41(5):batf207, 2025.
     https://doi.org/10.1093/bioinformatics/btaf207
 * [`mcl`](https://micans.org/mcl/) Stijn van Dongen, Graph Clustering Via a
     Discrete Uncoupling Process, SIAM Journal on Matrix Analysis and
     Applications, 30(1):121-141, 2008.
     http://link.aip.org/link/?SJMAEL/30/121/1
 * [`pymarkovclustering`](https://github.com/moshi4/pyMarkovClustering),
     implementation of MCL in Python
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
 * [`oxli`](https://github.com/oxli-bio/oxli)

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
