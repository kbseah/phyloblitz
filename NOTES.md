## Pipeline design notes


### Generating read clusters for assembly

I was hoping that miniasm would generate unitigs directly from minimap2 ava
mappings, but generally I get zero reads remaining after the initial
prefiltering step. The parameters probably need to be adjusted because miniasm
assumes we want to generate a genomic assembly whereas here we aim for a
targeted local assembly with ~1.5 kbp contig/unitig lengths.

Filter the AVA PAF to remove alignments that do not overlap over most of length
(relative to the shorter sequence of a pair). Especially for the 'flanking'
mode.


### Defining read clusters

How to define read clusters? Maybe use MCL on the minimap paf output: decide
which field to use as the mapping 'score'.

I observe that if we don't split clusters cleanly, they may contain more than
one species. SPOA appears to favour the majority, so even if the cluster
members are a mixture of different species the resulting consensus will be a
good match for the majority sequence, rather than a chimeric mixture.

If nominal sequence accuracy is Q20 (99%), then the expected sequence
divergence from pairwise mapping of two reads will be about 1 - 0.99*0.99 =
0.02; we therefore cut off clusters at about 2 to 3% divergence.


### Distinguishing strains

If we extract aligned read segments from the initial mapping plus some flanking
sequence, we may be able to distinguish strains because we expect there to be
more variation in the flanking intergenic sequence than in the SSU rRNA
sequence itself. This might work for bacteria where there is ~1 copy of the
rRNA genes per genome, but would not work so well for eukaryotes and species
with multiple copies

See recent survey of rRNA copy number and strain variation in prokaryotes: 10.1128/aem.02108-22

Why not extract ITS (prokaryotes) and ITS1, ITS2 (eukaryotes) for clustering too?

Run phyloblitz with SSU and LSU databases, identify reads with adjacent SSU and
LSU fragments to extract ITS. Assemble SSU and LSU separately because (i)
reference DBs are different, (ii) arrangement of the gene cluster differs
between prokaryotes, eukaryotes, organelles, and has exceptions.


### Secondary and supplementary mappings

phyloblitz uses minimap2, set to skip reporting of secondary alignments, but to
report supplementary alignments with soft clipping instead of hard clipping.

In conventional read mapping, metagenomic reads (query) are mapped to a
reference genome (subject), so secondary mappings represent alternative
positions in the genome that a read may map to, e.g., paralogs or repeat
sequences. Supplementary mappings represent parts of the read that may map to
locations in the reference that are discontinuous with the primary mapping; the
SAM specification refers to them as 'chimeric alignments'.

In the context of phyloblitz, long reads (query) are mapped to reference marker
gene database (subject) where the subjects are on average shorter than the
queries, so secondary mappings represent alternative reference sequences that a
read may map to (e.g. close relatives), while supplementary mappings represent
multiple paralogs or copies of the marker gene on the same read. With long
reads in the ~kbp range, it is possible that a single read may span two copies
of the SSU rRNA gene (~1.8 kbp) especially in eukaryotes, where rRNA operons
are often present as tandem repeat arrays.

TODO: Skip supplementary mappings that are also secondary to other
supplementary mappings, i.e. they overlap with another supplementary mapping.


## Notes

phyloblitz is essentially running nanopore amplicon pipeline, but first mapping
out the reads from a metagenome.

existing pipelines for working on amplicons -- either reference-based taxonomic
profiling or so-called LACA (longread amplicon consensus assembly)

we could extract reads and then just run them through existing pipeline

but phyloblitz provides an opportunity to extract sequence flanking the SSU
rRNA to resolve strains, a step short of complete de novo assembly

if we do not extract the flanking sequence, then the consensus sequences
assembled are a meta-consensus of related strains; unable to distinguish easily
because of read error > SSU rRNA sequence conservation

But read coverage per strain in a metagenome is expected to be low. The SSU
rRNA will always be an "island of conservation" in a variable genomic context.
With noisy reads we expect to seldom have enough coverage to assemble flanking
sequence accurately for most strains. What we can do instead: cluster the SSU
extracted sequences and assemble a consensus SSU per cluster, then estimate the
strain diversity represented by that consensus sequence from the flanking
regions through some kind of quick k-mer diversity approach. Do not attempt to
assemble the flanking regions.

How about iterative fishing assembler?


## Comparison with other tools

From metagenomes:

PenguiN assembler https://github.com/soedinglab/plass designed as a
strain-resolved de novo metagenomic assembler for viruses and microbial 16S
rRNA; limited to short reads (for now), closest to the phyloFlash targeted
marker assembly approach

SingleM https://wwood.github.io/singlem/ uses protein coding markers, so far
only applied to prokaryotes and dsDNA viruses. Produces taxonomy summary but
the actual marker sequences themselves are not directly accessible.

MetaMaps https://github.com/DiltheyLab/MetaMaps
uses MashMap to approximately map long reads to reference DB then estimates
sequence composition with an EM-based method.

MetaPhlAn https://github.com/biobakery/MetaPhlAn for shotgun metagenomes,
appears to be limited to short reads

K-mer hashing tools, such as sendsketch.sh from BBtools (see its [user
guide](https://bbmap.org/docs/guides/BBSketchGuide.md)): Much quicker lookup
compared to conventional read mapping, powerful tools for searching large
sequence collections. However results are usually reported as a summary of
references matched, reads with matches are not reported. Most sketch databases
are for whole genomes (although JGI does have a sketch server built from SILVA
database).

For working with amplicon libraries:

Nanoclust https://github.com/genomicsITER/NanoCLUST for 16S amplicon clustering
and assembly, no longer maintained and does not work with Nextflow >22.0
https://www.reddit.com/r/bioinformatics/comments/1j3860d/pipelines_for_metagenomics_nanopore_data/

Emu https://github.com/treangenlab/emu for long reads, implementing EM
algorithm for read classification, only outputs a summary table, may
over-classify

phyloblitz statement of novelty: works with long reads, lightweight database
and fast execution, works with eukaryotes and organelles, better with divergent
lineages because SSU is more conserved (flip side is that it is less sensitive
to strain differences), still can give usable results with noisy data. With
full-length assembled consensus, user can perform more detailed phylogenetic
analysis, not just stop at taxonomy summary. Also useful for lower-coverage
components of the metagenome.

Intended users: Working with low- to moderate-diversity metagenomes, aim to do
QC to determine metagenome composition before de novo assembly and binning,
interested in assembly of specific members of the metagenome, or may be
checking for contaminants in a single-species genome sequencing project.
