# phyloblitz 🌳🌩️

**Rapid SSU rRNA marker gene screening of long read metagenomes.**

[![GitHub Tag](https://img.shields.io/github/v/tag/kbseah/phyloblitz)](https://github.com/kbseah/phyloblitz/tags)
[![GitHub License](https://img.shields.io/github/license/kbseah/phyloblitz)](https://github.com/kbseah/phyloblitz?tab=MIT-1-ov-file)
[![Conda Version](https://img.shields.io/conda/vn/bioconda/phyloblitz)](https://anaconda.org/channels/bioconda/packages/phyloblitz/overview)


Note: This tool is still under development and the interface may change without
warning. We appreciate any bug reports or feedback; please [create an issue on
GitHub](https://github.com/kbseah/phyloblitz/issues/new). Please understand
that, because of limited resources, we cannot guarantee a quick response.

-----

`phyloblitz` is a tool to rapidly screen long read metagenomes with ribosomal
RNA genes, generating a taxonomic summary and targeted assemblies of
full-length SSU or LSU rRNA sequences.

 * Quickly check the composition and complexity of metagenomes before
     committing time and resources to de novo assembly.
 * Leverage the broad phylogenetic coverage of the SILVA rRNA database to
     detect taxa underrepresented in genome databases.
 * Use assembled gene sequences for phylogenetics and probe design.

The rRNA-targeted screening approach for metagenomes was originally implemented
for Illumina reads in the software tool
[phyloFlash](https://hrgv.github.io/phyloFlash/); the name `phyloblitz`
references phyloFlash and [bioblitzes](https://en.wikipedia.org/wiki/BioBlitz).


## Quick start


### Install

Using Conda/Mamba or Pixi, add the `conda-forge` and `bioconda` channels, and
install the `phyloblitz` package from Bioconda.

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n phyloblitz_env phyloblitz
```


### Download reference database

`phyloblitz` has a `download` subcommand to fetch the preformatted reference
databases for either SSU or LSU rRNA.

```bash
phyloblitz download --which_db SSU --db_version latest --outdir pbz_db
```

If you are unable to access the internet from the command line, you can
download the files with a web browser from
[Zenodo](https://doi.org/10.5281/zenodo.18627380) and transfer them later.


### Usage

`phyloblitz` has three subcommands:
 * `download` downloads the preformatted reference database
 * `run` applies the pipeline to a single read library
 * `compare` compares the run results between different libraries

Refer to the help messages for details on the phyloblitz run parameters:

```bash
phyloblitz --help
phyloblitz run --help
```

The following command runs phyloblitz on a Nanopore Q20+ dataset
`my_reads.fastq.gz`, and writes output files to the folder `pbz_out/`:

```bash
phyloblitz run --threads 12 --platform lr:hq \
--db pbz_db/SILVA_SSURef_NR99.nocontam.masked.trimmed.NR99.fasta \
--reads my_reads.fastq.gz --outdir pbz_out
```

A human-readable report is written in HTML and Markdown formats, along with a
JSON file containing run data that can be parsed by the `phyloblitz compare`
subcommand to compare composition of different samples.


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


### Run pipeline

 * Map reads to SILVA with minimap2 and retain only mapped reads
 * Summarize taxonomy across mapped reads, using consensus taxonomy of each read
 * Extract aligned portion of mapped reads, including supplementary alignments
   to account for multiple marker genes in one read. Also extract flanking
   sequence context to identify potential strain diversity.
 * Cluster extracted segments, either with isONclust3 or minimap2 + mcl
 * Assemble a consensus sequence per read cluster with spoa
 * Generate metrics per cluster for diagnostics: expect similar error rate per
     read vs. consensus, identify clusters with too few reads (TODO)


### Compare pipeline

 * Pool read segments containing marker of interest that were extracted from
   individual phyloblitz runs.
 * Cluster and assemble each cluster as above.
 * Report closest taxon hits for each assembled marker sequence, and the number
   of reads per sample represented in that cluster.


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
