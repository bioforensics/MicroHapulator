# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## Unreleased

### Changed
- Updated mybinder demo (see #69).
- Simulated Illumina sequencing now uses 1 thread by default, which paradoxically should lead to better performance (#71).
- Moved panel definition code moved out of the core code and into dedicated notebooks (#74).
- Replaced `MissingBAMIndexError` with BAM auto-indexing code (#78).

### Fixed
- Corrected a bug with Fastq headers in `mhpl8r seq` module (see #71).
- Corrected a bug resulting from attempting to do set operations on `None` (see #75).

## [0.4] 2019-11-05

### Added
- New `contain` module for calculating the proportion of alleles from one sample present in another sample (see #41).
- New `--dry-run` mode for `sim` and `mixture` modules (see #41).
- New `prob` module for calculating likelihood ratio tests based on the random match probability (see #43).
- New `seq` module focused entirely on sequencing samples where profiles have already been simulated (see #45).
- New `mix` module for merging simulated profiles into a simulated mixture sample profile (see #45).
- New `unite` module for "mating" two profiles to create a simulated "offspring" profile (see #47).
- New `diff` modfule for showing the differences between two profiles (see #58, #60).

### Changed
- Huge refactoring effort to accommodate for recent changes to MicroHapDB's Python API (see #66).
- The `sim` module no longer performs simulated sequencing (now handled by new `seq` module) and instead focuses entirely on haplotype simulation (see #45).
- The `type` module now dynamically selects either an automatic threshold or a fixed threshold based on effective coverage (see #51, #61).
- Moved simulation scripts to notebook directory, reimplemented as a Snakemake workflow (see #50, #57, #59).

### Fixed
- Corrected a bug in the USA panel code ensuring that all loci have allele frequency data for all relevant populations (see #56).
- Corrected a problem with the `contrib` API that prevented it from working directly on `Profile` objects (see #70).

### Removed
- Dropped the `mixture` module, whose functionality is now covered by the more granular `sim`, `mix`, and `seq` modules (see #45).
- Dropped the `LocusContext` class in favor of MicroHapDB's `TargetAmplicon` class designed for a similar purpose (see #67).


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
