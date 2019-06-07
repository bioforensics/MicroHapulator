# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## Unreleased

### Added
- New `contain` module for calculating the proportion of alleles from one sample present in another sample (see #41).
- New `--dry-run` mode for `sim` and `mixture` modules (see #41).


## [0.3] 2019-06-06

### Added
- New interactive demo using Jupyter notebook and mybinder.org (see #32).
- Scripts for simulating 8k individuals with demographics roughly matching those of the U.S.A. (see #13).

### Changed
- Updated `dist` and `genotype` modules to better support equality and distance comparisons between simulated and observed genotypes (see #34).
- Simulated and observed genotype data now use a single unified JSON representation, enforced by JSON schema validation (see #38).


## [0.2] 2019-05-21

### Added
- New `getrefr` module for installing human reference genome (see #30).
- Improved error handling related to BAM index files (see #31).

### Changed
- Replaced pip/PyPI installation instructions with bioconda installation instructions (see 1b61c59200).

### Fixed
- Corrected usage statement for `mhpl8r contrib` (see #31).


## [0.1.3] 2019-05-14

- mark failing test
- exclude `__init__.py` from test invocation


## [0.1.2] 2019-05-07

Updating file manifest, troubleshooting bioconda build.


## [0.1.1] 2019-05-07

Minor fix to `setup.py` for distribution purposes.


## [0.1] 2019-05-02

Initial release of MicroHapulator!

Command-line entry point: `mhpl8r`
Command-line operations:
- `contrib`: estimate number of contributors in a sample
- `dist`: compute naive Hamming distance between two sample genotypes
- `mixture`: simulate a multiple-contributor sample
- `refr`: retrieve reference genome cutouts for a specified microhap panel
- `sim`: simulate a single-contributor sample
- `type`: infer genotype from mapped reads
