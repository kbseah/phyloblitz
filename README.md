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

phyloblitz is distributed via [Bioconda](https://bioconda.github.io/). Read
[this explainer](https://conda.org/blog/2024-08-14-conda-ecosystem-explained/)
if you are new to the Conda ecosystem. I recommend either installing
Conda/Mamba with [Miniforge](https://conda-forge.org/download/) or using
[pixi](https://pixi.prefix.dev/latest/).


### Install with Conda/Mamba

Set up your Conda/Mamba configuration as recommended for Bioconda, if you have
not already done so:

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

Install phyloblitz to a new environment (here named `phyloblitz_env`):

```bash
conda create -n phyloblitz_env phyloblitz
```

Activate environment and view help message:

```bash
conda activate phyloblitz_env
phyloblitz --help
```


### Install with Pixi

Configure [pixi](https://pixi.prefix.dev/latest/) to use Bioconda channel in
addition to conda-forge:

```bash
pixi config set default-channels '["conda-forge", "bioconda"]'
```

Create a new pixi workspace in a folder named `phyloblitz_workspace` and
install phyloblitz there:

```bash
mkdir phyloblitz_workspace
cd phyloblitz_workspace
pixi init
pixi add phyloblitz
```

Start a pixi shell session and view help message:

```bash
pixi shell
phyloblitz --help
```

Use Ctrl-D or `exit` to exit the pixi shell session.


### Install as a container from BioContainers

Bioconda packages are automatically containerized and uploaded to the
[BioContainers](https://biocontainers.pro/tools/phyloblitz) registry, so you
could simply pull the container with either [Docker](https://www.docker.com/),
[Singularity](https://sylabs.io/docs/), or [Apptainer](https://apptainer.org/)
(my preference for HPC environments):

```bash
docker pull quay.io/biocontainers/phyloblitz:0.2.0--pyhdfd78af_0
```

```bash
singularity pull phyloblitz.sif docker://quay.io/biocontainers/phyloblitz:0.2.0--pyhdfd78af_0
```


### Development version

**Note: Currently not accepting pull requests. Please create an issue if you
want to contribute.**

Build and run time dependencies are managed with [pixi](https://pixi.sh/). When
`pixi shell` is run for the first time in the folder containing the `pixi.toml`
configuration file, pixi will resolve and install dependencies. Run
`phyloblitz` without arguments or with the `--help` parameter to view help
message.

```bash
git clone git@github.com:kbseah/phyloblitz.git
cd phyloblitz
pixi shell # set up workspace and start pixi shell session
phyloblitz --help
```

Use `exit` or `Ctrl-D` to exit pixi shell session.

Required inputs are a preprocessed SILVA reference database and the long read
files in fastq(.gz) format. After navigating to the phyloblitz folder, fetch
the (SSU rRNA) database file from Zenodo on the command line with:

```bash
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

```bash
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
 * Extract aligned portion of mapped reads, including supplementary alignments
   to account for multiple marker genes in one read. Also extract flanking
   sequence context to identify potential strain diversity.
 * Cluster extracted segments, either with isONclust3 or minimap2 + mcl
 * Assemble a consensus sequence per read cluster with spoa
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
