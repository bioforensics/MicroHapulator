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


## Dependencies

MicroHapulator depends on several Python packages.
These are listed in `setup.py` in the main source code distribution.
If performing a non-standard installation with `pip`, these dependencies will automatically be installed from the Python Package Index (PyPI).

Preparing Illumina reads for analysis and interpretation with MicroHapulator also depends on several bioinformatics tools, including [FLASH](https://ccb.jhu.edu/software/FLASH/), [Minimap2](https://lh3.github.io/minimap2/minimap2.html) (or an alternative short read aligner), [SAMtools](http://www.htslib.org/), and [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
These software tools are installed automatically when performing a standard installation with Conda.
Otherwise, they must be installed manually.


## Development quick start

If you're setting up an environment for developing MicroHapulator, you may want to skip the procedure outlined above and use the following instead.

```bash
conda create --new microhapulator python=3.11 flash minimap2 samtools fastqc
conda activate microhapulator
git clone https://github.com/bioforensics/MicroHapulator.git
cd MicroHapulator/
pip install -e .  # Install the package and its Python dependencies
make devdeps      # Install development packages
make devhooks     # Register pre-commit hooks for development
```
