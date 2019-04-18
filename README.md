[![MicroHapulator build status][travisbadge]](https://travis-ci.org/bioforensics/MicroHapulator)
[![Test coverage][codecovbadge]](https://codecov.io/github/bioforensics/MicroHapulator)
[![BSD licensed][licensebadge]](https://github.com/bioforensics/MicroHapDB/blob/master/LICENSE.txt)

# MicroHapulator

Daniel Standage, 2018-2019  
https://github.com/bioforensics/microhapulator

**MicroHapulator** is a package for simulating and analyzing microhaplotype sequence data for forensic analysis.
In addition to simulating targeted sequencing of selected microhap loci, it can be used for calling genotypes and performing sample matching and mixture analysis.
MicroHapulator relies on microhap annotations and population allele frequencies from [MicroHapDB](https://github.com/bioforensics/microhapdb) and the Illumina error models included with [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq/).


## Installation

Installation with bioconda is pending.
For now, install the following installation procedure is recommended.

```
pip install git+https://github.com/bioforensics/happer.git
pip install git+https://github.com/bioforensics/MicroHapDB.git
pip install git+https://github.com/bioforensics/MicroHapulator.git
```


## Usage

<!-- replaceme:start -->
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
  subcmd             contrib, dist, refr, sim, type

Global arguments:
  -h, --help         show this help message and exit
  -v, --version      show program's version number and exit
  -l F, --logfile F  log file for diagnostic messages, warnings, and errors
  --tee              write diagnostic output to logfile AND terminal (stderr)

```
<!-- replaceme.end -->


## Demo

```
# Grab sequences of microhap panel "beta"
mhpl8r refr --out panel.fasta hg38.fasta beta

# Simulate 100000 reads from a European American individual
mhpl8r sim --panel beta --num-reads 100000 --out myreads.fastq hg38.fasta SA004250L

# Align reads to reference genome
bwa mem hg38.fasta myreads.fastq | samtools view -bS | samtools sort -o myreads.bam -
samtools index myreads.bam

# Make genotype calls
mhpl8r type --out genotype.json panel.fasta myreads.bam

# Estimate number of contributors in the sample
mhpl8r --json genotype.json
```


[travisbadge]: https://img.shields.io/travis/bioforensics/MicroHapulator.svg
[codecovbadge]: https://img.shields.io/codecov/c/github/bioforensics/MicroHapulator.svg
[licensebadge]: https://img.shields.io/badge/license-BSD-blue.svg
