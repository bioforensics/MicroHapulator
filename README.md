[![MicroHapulator build status][travisbadge]](https://travis-ci.org/bioforensics/MicroHapulator)
[![Test coverage][codecovbadge]](https://codecov.io/github/bioforensics/MicroHapulator)
![Platform support][platformbadge]
[![install with bioconda][biocondabadge]](http://bioconda.github.io/recipes/microhapulator/README.html)
[![BSD licensed][licensebadge]](https://github.com/bioforensics/MicroHapDB/blob/master/LICENSE.txt)

# MicroHapulator

Daniel Standage, 2018-2019  
https://github.com/bioforensics/microhapulator

**MicroHapulator** is a package for analyzing and interpreting microhaplotype sequence data for forensic analysis.
It provides tools for genotyping as well as deterministic and probabilistic interpretation of simple and complex samples.
It also provides tools for constructing mock genotypes and simulating MPS sequencing of simple and complex samples.
MicroHapulator relies on microhap annotations and population allele frequencies from [MicroHapDB](https://github.com/bioforensics/microhapdb) and the MPS error models included with [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq/).


## Installation

Installation with bioconda is recommended.

```
conda install -c bioconda microhapulator
```

> ***NOTE**: If you'd prefer to install with pip, see [.travis.yml](.travis.yml) and [setup.py](setup.py) for hints.*


## Interactive demo

Click the badge below to launch a quick interactive demo of MicroHapulator.

[![Binder][binderbadge]](https://mybinder.org/v2/gh/bioforensics/MicroHapulator/master?filepath=binder%2FREADME.md)


## Usage

```
usage: mhpl8r [-h] [-v] [-l F] [--tee] subcmd ...

≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
  __  __ _            _  _                _      _
 |  \/  (_)__ _ _ ___| || |__ _ _ __ _  _| |__ _| |_ ___ _ _
 | |\/| | / _| '_/ _ \ __ / _` | '_ \ || | / _` |  _/ _ \ '_|
 |_|  |_|_\__|_| \___/_||_\__,_| .__/\_,_|_\__,_|\__\___/_|
                               |_|
≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
Invoke `mhpl8r <subcmd> --help` and replace `<subcmd>` with one of the
subcommands listed below to see instructions for that operation.

Subcommands:
  subcmd             contain, contrib, diff, dist, mix, prob, seq, sim, type,
                     unite

Global arguments:
  -h, --help         show this help message and exit
  -v, --version      show program's version number and exit
  -l F, --logfile F  log file for diagnostic messages, warnings, and errors
  --tee              write diagnostic output to logfile AND terminal (stderr)

```


## Contributing

We welcome contributions from the community!
Feel free to ask questions, make suggestions, or report bugs using the [issue tracker](https://github.com/bioforensics/MicroHapulator/issues).
If you are interested in submitting patches, [docs/DEVEL.md](docs/DEVEL.md) contains a few suggestions for a development setup.
All contributors are expected to abide by the project's [Code of Conduct](docs/CONDUCT.md).


[travisbadge]: https://img.shields.io/travis/bioforensics/MicroHapulator.svg
[codecovbadge]: https://img.shields.io/codecov/c/github/bioforensics/MicroHapulator.svg
[licensebadge]: https://img.shields.io/badge/license-BSD-blue.svg
[biocondabadge]: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
[binderbadge]: https://mybinder.org/badge_logo.svg
[platformbadge]: https://img.shields.io/badge/Platforms-linux--64%2Cosx--64-orange.svg
