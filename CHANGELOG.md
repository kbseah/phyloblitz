# CHANGES

## v0.1.0

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

## v0.2.0

 - [x] Check if lr:hq mode is better for Q20 ONT reads than map:ont
 - [x] Option to cluster with isONclust3 instead of minimap2 ava
 - [x] Sequence diversity in reads flanking the SSU sequence per cluster to
       identify meta-consensus sequences
 - [x] Separately cluster flanking sequences as measure of diversity
 - [x] Consensus sequence length in table
 - [x] Reference database on Zenodo for SSU and LSU
 - [x] Account for multiple SSUs in one read; see supplementary alignments


## v0.3.0

 - [x] Pretty log messages
 - [x] Deprecate two-pass mode
 - [x] Subcommands for run and database download
 - [x] Database download from Zenodo with checksums
 - [x] Parse supplementary alignments by default
 - [x] Skip duplicate segments extracted from overlapping alignments
 - [x] Downsample clusters with too many reads for assembly (~500 should be ok)
 - [x] Multisample comparison and co-assembly
 - [x] Embed graphics as PNGs into HTML
 - [x] Detailed documentation


## v1.0.0 targets

 - [ ] Phylogenetic tree of assembled sequences
 - [ ] Compare samples by Unifrac distances from tree of assembled sequences
 - [ ] Interactive graphics, deprecate Markdown report
 - [ ] Fix the CSS stylesheet
 - [ ] Greengenes2 database


## Future plans/ideas

 - [ ] MultiQC integration
 - [ ] Option to search final assembled cluster sequences with vsearch instead
       of minimap2; better handling of divergent sequences: Example, the
       Aestuariibacter sequence in SRR28830816
 - [ ] Divergence of reads vs consensus in cluster to detect chimeras
 - [ ] Does extracting flanking sequence context improve strain resolution?
     - [ ] Filter out all-vs-all hits with overhangs
     - [ ] Investigate effect of dv cutoff values and clustering methods
 - [ ] Benchmarking against defined test datasets
 - [ ] Extract the ITS too?
