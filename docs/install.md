# Installation

## Overview

MicroHapulator depends on several Python packages and bioinformatics tools.
Using conda to install and manage these packages and their dependencies is recommended.
See [this page](https://bioconda.github.io/user/install.html#set-up-channels) for instructions on configuring conda to enable installing packages from the bioconda channel.

```
conda create --name microhapulator -y python=3.11 microhapulator
conda activate microhapulator
mhpl8r --help
```

To test whether installation was successful, pytest is recommended.

```
conda install -y pytest
pytest --pyargs microhapulator
```


## Development quick start

If you're setting up an environment for developing MicroHapulator, you probably want to skip the procedure outlined above and use the following instead, after [installing pixi](https://pixi.sh/latest/installation/).

```bash
git clone https://github.com/bioforensics/MicroHapulator.git
cd MicroHapulator/
pixi install
pixi install -e test
pixi task list
pixi run hooks
```
