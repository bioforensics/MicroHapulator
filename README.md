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
[MicroHapDB](https://github.com/bioforensics/microhapdb) can be used to prepare marker definitions, population frequencies, and reference sequences—see [the configuration documentation](https://microhapulator.readthedocs.io/en/latest/config.html) for more details.


## Installation

Installation with the [Conda package manager](https://docs.conda.io/en/latest/) is recommended.
Full installation instructions are available [here](https://microhapulator.readthedocs.io/en/latest/install.html).

```
conda install -c bioconda microhapulator
```


## Interactive demo

Click the badge below to launch a quick interactive demo of MicroHapulator.

- [![Binder][binderbadge]](https://mybinder.org/v2/gh/bioforensics/MicroHapulator/master?filepath=binder%2Fdemo-cli.ipynb) MicroHapulator CLI (`mhpl8r --help`)
- [![Binder][binderbadge]](https://mybinder.org/v2/gh/bioforensics/MicroHapulator/master?filepath=binder%2Fdemo-api.ipynb) MicroHapulator Python API (`import microhapulator`)
- [![Binder][binderbadge]](https://mybinder.org/v2/gh/bioforensics/MicroHapulator/master?filepath=binder%2Fdemo-sim.ipynb) Mock data simulation


## Usage

The [MicroHapulator user manual](https://microhapulator.readthedocs.io/en/latest/manual.html) is a good place to start.
There is also a reference for the [command line interface](https://microhapulator.readthedocs.io/en/latest/cli.html) and the [Python API](https://microhapulator.readthedocs.io/en/latest/api.html).


## Contributing

MicroHapulator was created by Daniel Standage, with the support of the National Bioforensic Analysis Center (NBFAC).
Contributions from the broader community are welcomed!
Feel free to ask questions, make suggestions, or report bugs using the [issue tracker](https://github.com/bioforensics/MicroHapulator/issues).
If you are interested in submitting patches, [the developer documentation](https://microhapulator.readthedocs.io/en/latest/devel.html) contains a few suggestions for a development setup.
All contributors are expected to abide by the project's [Code of Conduct](https://microhapulator.readthedocs.io/en/latest/conduct.html).


[travisbadge]: https://img.shields.io/travis/bioforensics/MicroHapulator.svg
[codecovbadge]: https://img.shields.io/codecov/c/github/bioforensics/MicroHapulator.svg
[licensebadge]: https://img.shields.io/badge/license-BSD-blue.svg
[biocondabadge]: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
[binderbadge]: https://mybinder.org/badge_logo.svg
[platformbadge]: https://img.shields.io/badge/Platforms-linux--64%2Cosx--64-orange.svg
