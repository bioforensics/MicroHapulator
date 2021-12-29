[![MicroHapulator build status][travisbadge]](https://travis-ci.org/bioforensics/MicroHapulator)
[![Test coverage][codecovbadge]](https://codecov.io/github/bioforensics/MicroHapulator)
![Platform support][platformbadge]
[![install with bioconda][biocondabadge]](http://bioconda.github.io/recipes/microhapulator/README.html)
[![BSD licensed][licensebadge]](https://github.com/bioforensics/MicroHapDB/blob/master/LICENSE.txt)

# MicroHapulator

Daniel Standage, 2018-2021
https://github.com/bioforensics/microhapulator

**MicroHapulator** is software for empirical haplotype calling, analysis, and basic forensic interpretation of microhaplotypes from Illumina sequence data.
Also included are tools for constructing mock genotypes and simulating Illumina sequencing of single- and multiple-contributor samples for testing and development purposes.
MicroHapulator can be configured for use with any published or custom panel.
[MicroHapDB](https://github.com/bioforensics/microhapdb) can be used to prepare marker definitions, population frequencies, and reference sequences—see [the configuration documentation](google.com) for more details.


## Installation

Installation with the [Conda package manager](https://docs.conda.io/en/latest/) is recommended.
Full installation instructions are available [here](google.com).

```
conda install -c bioconda microhapulator
```


## Interactive demo

Click the badge below to launch a quick interactive demo of MicroHapulator.

- [![Binder][binderbadge]](https://mybinder.org/v2/gh/bioforensics/MicroHapulator/master?filepath=binder%2Fdemo-cli.ipynb) MicroHapulator CLI (`mhpl8r --help`)
- [![Binder][binderbadge]](https://mybinder.org/v2/gh/bioforensics/MicroHapulator/master?filepath=binder%2Fdemo-api.ipynb) MicroHapulator Python API (`import microhapulator`)
- [![Binder][binderbadge]](https://mybinder.org/v2/gh/bioforensics/MicroHapulator/master?filepath=binder%2Fdemo-sim.ipynb) Mock data simulation


## Usage

```
usage: mhpl8r [-h] [-v] subcmd ...

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
  subcmd         balance, contain, contrib, diff, dist, mix, prob, seq, sim,
                 type, unite

Global arguments:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit

```


## Contributing

MicroHapulator was created by Daniel Standage, with the support of the National Bioforensic Analysis Center (NBFAC).
But we welcome contributions from the broader community!
Feel free to ask questions, make suggestions, or report bugs using the [issue tracker](https://github.com/bioforensics/MicroHapulator/issues).
If you are interested in submitting patches, [the developer documentation](google.com) contains a few suggestions for a development setup.
All contributors are expected to abide by the project's [Code of Conduct](google.com).


[travisbadge]: https://img.shields.io/travis/bioforensics/MicroHapulator.svg
[codecovbadge]: https://img.shields.io/codecov/c/github/bioforensics/MicroHapulator.svg
[licensebadge]: https://img.shields.io/badge/license-BSD-blue.svg
[biocondabadge]: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
[binderbadge]: https://mybinder.org/badge_logo.svg
[platformbadge]: https://img.shields.io/badge/Platforms-linux--64%2Cosx--64-orange.svg
