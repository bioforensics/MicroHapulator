# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

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
