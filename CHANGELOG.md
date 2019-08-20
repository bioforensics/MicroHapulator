# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## Unreleased

### Added
- New `contain` module for calculating the proportion of alleles from one sample present in another sample (see #41).
- New `--dry-run` mode for `sim` and `mixture` modules (see #41).
- New `prob` module for calculating likelihood ratio tests based on the random match probability (see #43).
- New `seq` module focused entirely on sequencing samples where genotypes have already been simulated (see #45).
- New `mix` module for merging simulated genotypes into a simulated mixture sample (see #45).
- New `unite` module for "mating" two genotypes to create a simulated "offspring" genotype (see #47).
- New `diff` modfule for showing the differences between two genotypes (see #58, #60).

### Changed
- The `sim` module no longer performs simulated sequencing (now handled by new `seq` module) and instead focuses entirely on haplotype simulation (see #45).
- The `type` module now dynamically selects either an automatic threshold or a fixed threshold based on effective coverage (see #51, #61).
- Moved simulation scripts to notebook directory, reimplemented as a Snakemake workflow (see #50, #57, #59).

### Fixed
- Corrected a bug in the USA panel code ensuring that all loci have allele frequency data for all relevant populations (see #56).

### Removed
- Dropped the `mixture` module, whose functionality is now covered by the more granular `sim`, `mix`, and `seq` modules (see #45).


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
