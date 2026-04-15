# Reference database

The database files for `phyloblitz` are derived from the [SILVA ribosomal RNA
database](https://www.arb-silva.de/), with some additional processing (see
below). Sequences in SILVA have been curated and classified with a consistent
taxonomy for all domains of life. It is a commonly used reference for probe
design and research in molecular ecology and phylogenetics.

## Download

`phyloblitz` has a `download` subcommand to fetch the preformatted reference
databases for either SSU or LSU rRNA.

```bash
phyloblitz download --which_db SSU --db_version latest --outdir pbz_db
```

If you are unable to access the internet from the command line, you can
download the files with a web browser from
[Zenodo](https://doi.org/10.5281/zenodo.18627380) and transfer them later.


## Preprocessing

The original SILVA files were processed with additional filtering, trimming of
contaminant sequences, and masking of low-complexity repeats to avoid excessive
false-positive matches. Use the appropriate file for either SSU or LSU rRNA
genes. Refer to the [database processing
pipeline](https://github.com/kbseah/phyloblitz-db) for the actual commands used
to prepare the database.

If you build your own database, low complexity regions should be masked (e.g.
with `bbmask.sh`) otherwise mapping will be slow. The database processing steps
were originally developed for phyloFlash; refer to the phyloFlash paper
([Gruber-Vodicka, Seah & Pruesse,
2021](https://doi.org/10.1128/mSystems.00920-20)) for details and discussion.


## Format

Database files are in FASTA format.

Header lines are in the standard SILVA format:

 * Accession comprising original ENA accession and start/end coordinates of the
   rRNA marker sequence within the ENA record, separated by `.`
 * Single space character
 * Taxonomy string from higher to lower taxon names, separated by `;` without
   spaces, however taxon names may contain spaces.

Sequences use DNA bases (T instead of U); masked bases are lower-case `n`.
