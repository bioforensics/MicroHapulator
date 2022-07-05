# MicroHapulator

**MicroHapulator** is a software package designed for analysis of sequence data from NGS-based microhaplotype assays.
Its core features include an empirical haplotype calling algorithm and tools for filtering haplotypes by read count.
The MicroHapulator package also provides an end-to-end quality control and analysis workflow that is invoked with a single command:
outputs of this workflow include a comprehensive quality control report as well as profiles suitable for interpretation with probabilistic interpretation programs such as [LRmix Studio](https://github.com/smartrank/lrmixstudio) or [EuroForMix](http://www.euroformix.com/).

Beyond its core features, MicroHapulator also operations supporting basic forensic interpretation (e.g. estimating number of contributors, computing random match probabilities, performing basic likelihood ratio tests) and tools for constructing mock genotypes and simulating Illumina sequencing of simple and complex samples.

MicroHapulator can be configured for use with any published or custom panel.
[MicroHapDB](https://github.com/bioforensics/microhapdb) can be used to prepare marker definitions, population frequencies, and reference sequences—see [the configuration documentation](config.md) for more details.


## Getting Started

- If you're interested in MicroHapulator's empirical haplotype calling procedure, see [this page](typing).
- If you want to give MicroHapulator a try before installing it on your machine, <a href="https://mybinder.org/v2/gh/bioforensics/MicroHapulator/main?filepath=binder%2Fdemo.ipynb" target="_blank">click here</a> to launch an interactive demo in the cloud with simulated data.
- If you're new to forensic DNA analysis, [this primer](primer.md) is the recommended entry point for the documentation.

The table of contents below includes additional links to various manuals, instructions, and reference materials.


## Table of Contents

```{toctree}
:maxdepth: 2

install
primer
typing
manual
config
devel
api
cli
conduct
```


## Indices and tables

- {ref}`genindex`
- {ref}`modindex`
- {ref}`search`
