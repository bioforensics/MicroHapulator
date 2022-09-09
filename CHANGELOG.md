# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## Unreleased

### Added
- Profiles compatible with probgen programs now included in pipe output (#135).
- Haplotype call plots now included in the pipe HTML report (#136).
- Added offtarget module to count reads that map to off target loci in hg38 (#143).
- Added typing rate and mapping rate information per marker to the main pipe HTML report (#146).
- Added marker detail HTML report (#146).

### Changed
- Exposed static and dynamic threshold configuration to pipe CLI (#135).
- Updated plot colors in the pipe HTML report (#139).
- Updated pipe HTML report to conditionally plot read length histograms or tables depending on uniformity (#140).


## [0.6.1] 2022-04-27

### Fixed
- Bug with Snakefile not being included in package data.


## [0.6] 2022-04-27

### Added
- New API and CLI entry points for computing and visualizing heterozygote balance (#122, #131).
- New `typing_rate` method for the TypingResult class (#127).
- New API function for plotting distribution of read lengths (#128).
- New CLI entry point for downloading GRCh38 (#130).
- New end-to-end microhap analysis pipeline and report (#129, #132).

### Changed
- Interlocus balance code updated to support generating high-resolution graphics and performing a chi-square goodness-of-fit test (#121, #131, #132).

### Fixed
- Bug with filtering/genotype calling for markers with no valid reads (#123).
- Set 0.01 as a more reasonable default frequency for rare alleles than 0.001 (#131, #132).


## [0.5] 2022-01-20

### Added
- New `--base-qual` parameter for `mhpl8r type` to set the minimum required base quality when iterating over reads in a pileup (#83).
- New `mhpl8r balance` subcommand for calculating and visualizing interlocus balance (#85).
- Users can now supply marker definitions, frequences, and reference sequences as TSV/FASTA files instead of MicroHapDB references (#93).
- Configuration file examples in `microhapulator/data/configs/` (#105).
- New `mhpl8r filter` subcommand with support for marker-specific thresholds (#113, #114).
- New `mhpl8r convert` subcommand for converting genotype calls into a format compatible with probgen tools such as LRMix Studio and EuroForMix (#115).

### Changed
- Updated mybinder demo (see #69, #110, #113).
- Simulated Illumina sequencing now uses 1 thread by default, which paradoxically leads to better performance (#71).
- Moved panel definition code moved out of the core code and into dedicated notebooks (#74).
- Replaced `MissingBAMIndexError` with BAM auto-indexing code (#78).
- Improved read names and choice of interleaved or paired output for `mhpl8r seq` (#80).
- Updated filtering of haplotype calls / typing results
    - Replaced `--threshold` argument with `--static` and `--dynamic`, disabled both by default (#82, #83).
    - Split `mhpl8r type` subcommand into `type` and `filter`, with `--static` and `--dynamic` arguments only relevant to the latter (#113).
- Changed the default pysam pileup `max_depth` parameter, overriding 8000 with 1e6 and exposing as a CLI parameter (#87, #113).
- Removed dependency on MicroHapDB for marker definitions, frequencies, and sequences (#93).
- Refactored CLI and Python API, adding new `microhapulator.api` module to serve as main entry point (#98, c98bf6c78ef4).
- Replaced the "ObservedProfile" terminology with the more appropriate "TypingResult" (#99).
- Documentation now uses Sphinx to render markdown as HTML (c98bf6c78e, #101, #102, #105, #106, #110).
- Updated JSON schema for simulated profiles and typing results (#109).

### Fixed
- Corrected a bug with Fastq headers in `mhpl8r seq` module (#71).
- Corrected a bug resulting from attempting to do set operations on `None` (#75).
- Corrected a bug with RMP implementation (#86).
- Set minimum versions for runtime dependencies (#97).
- Corrected the `--dynamic` filter to operate on total haplotype counts rather than average counts (#114).


## [0.4.1] 2019-11-06

Fix minor metadata typo


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
