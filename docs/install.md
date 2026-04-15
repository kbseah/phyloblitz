# Installation

phyloblitz is distributed via [Bioconda](https://bioconda.github.io/). Read
[this explainer](https://conda.org/blog/2024-08-14-conda-ecosystem-explained/)
if you are new to the Conda ecosystem. I recommend either installing
Conda/Mamba with [Miniforge](https://conda-forge.org/download/) or using
[pixi](https://pixi.prefix.dev/latest/).


## Install with Conda/Mamba

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


## Install with Pixi

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


## Install as a container from BioContainers

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
