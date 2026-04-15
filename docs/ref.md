# Reference database

`phyloblitz` has a `download` subcommand to fetch the preformatted reference
databases for either SSU or LSU rRNA.

```bash
phyloblitz download --which_db SSU --db_version latest --outdir pbz_db
```

If you are unable to access the internet from the command line, you can
download the files with a web browser from
[Zenodo](https://doi.org/10.5281/zenodo.18627380) and transfer them later.

The database files are derived from the SILVA database, with additional
filtering, trimming of contaminant sequences, and masking of low-complexity
repeats to avoid excessive false-positive matches. Use the appropriate file for
either SSU or LSU rRNA genes. Refer to the [database processing
pipeline](https://github.com/kbseah/phyloblitz-db) for the actual commands used
to prepare the database.

If you build your own database, low complexity regions should be masked (e.g.
with `bbmask.sh`) otherwise mapping will be slow. The database processing steps
were originally developed for phyloFlash; refer to the phyloFlash paper
([Gruber-Vodicka, Seah & Pruesse,
2021](https://doi.org/10.1128/mSystems.00920-20)) for details and discussion.
