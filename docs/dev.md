# Development

**Note: Currently not accepting pull requests. Please create an issue if you
want to contribute.**

Build and run time dependencies are managed with [pixi](https://pixi.sh/). When
`pixi shell` is run for the first time in the folder containing the `pixi.toml`
configuration file, pixi will resolve and install dependencies. Run
`phyloblitz` without arguments or with the `--help` parameter to view the help
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


## Testing

Test commands are managed as pixi run commands in the `pixi.toml` configuration
file.

Unit tests (under `tests/`) have only been implemented for a few key functions,
because the majority of the code handles external tools and file formats. Run
unit tests and code formatter (`black`) with `pixi run test`.

Small test datasets for Nanopore and PacBio (subsetted from published data on
SRA) are included in the repo at `runtest/data`. After downloading the database
file to the subfolder `138.2/` you can test `phyloblitz` on the test data with
`pixi run run` (Nanopore) and `pixi run run_pb` (PacBio). The output will be in
`runtest/run` or `runtest/run_pb`.

