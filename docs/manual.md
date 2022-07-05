# User Manual: End-to-end Workflow

## Overview

**MicroHapulator** is software for empirical haplotype calling and analysis of microhaplotype (MH) from NGS data.
It can be configured for use with any published or custom panel.
In addition to comprehensive quality control checks and basic interpretation features, MicroHapulator output can be readily formatted for compatibility with popular probabilistic genotyping programs for robust forensic interpretation.

The MicroHapulator software is run by calling the `mhpl8r` command in a UNIX shell terminal (using e.g. the Windows Subsystem for Linux or the "Terminal" program in MacOS).
Using the terminal can be intimidating for the uninitiated, so this manual is written to be acommodating of all users.
Forensic practicioners routinely perform tasks and analyses that are much more complex than MicroHapulator's interface, so we are confident that the command line will not be an obstacle to interested users.


## End-to-end Workflow Execution

Most users will want to run the entire MicroHapulator workflow end-to-end with a single command.
The `mhpl8r pipe` procedure accepts one or more samples as input and—resources permitting—will process multiple samples in parallel.
This command executes all core preprocessing and analysis workflow steps, as well as some additional quality control checks, and collates all results into a single human-readable report.

Prior to executing the workflow, the user will need to have the following inputs ready.

- MH marker reference sequences in FASTA format (see [the config docs](config.md))
- MH marker definitions (SNP offsets) in TSV format (see [the config docs](config.md))
- A pair of NGS read files in FASTQ format for each sample
    - Each FASTQ file name must include the name of the corresponding sample
    - Sample names must be unique, and one sample name cannot be contained in another sample name
    - All FASTQ files must be contained in a single directory
        - If desired, the directory may be further divided into subdirectories by e.g. experiment, run, project, case, and so on; MicroHapulator will scan all subdirectories to locate FASTQ files
        - It is OK if this directory includes FASTQ files for samples that the user does not want to analyze; any FASTQ file that does not match the user-provided sample name(s) will be ignored

> *MicroHapulator has been tested successfully on single-end Ion S5 data, but at the moment the end-to-end workflow has only been tested on paired-end reads from the Illumina MiSeq instrument.*

With all of this in place, the user can invoke the end-to-end workflow as follows.

```bash
mhpl8r pipe marker-refr.fasta marker-defn.tsv /path/to/fastq/dir/ SAMPLE1 SAMPLE2 --workdir caseXYZ
```

MicroHapulator will generate many intermediate data files throughout the workflow.
These will be stored in the *working directory* (`caseXYZ` in this example), organized by sample.
MicroHapulator will create this working directory if it doesn't already exist.
When the workflow completes, the final report will be located at `caseXYZ/report.html`.
This can be viewed in your favorite Web browser.
The typing result and genotype prediction for each sample will be located at `caseXYZ/analysis/SAMPLENAME/SAMPLENAME-type.json`.

The user has the option to configure various aspects of the workflow's behavior.
For a detailed description of these options, run `mhpl8r pipe --help` in the terminal or consult the [CLI docs](cli.md) for the `mhpl8r pipe` command.


## Stepwise Workflow Execution

While end-to-end execution is most convenient, users do have the option to execute individual tasks in a stepwise fashion, which gives maximum flexibility for workflow configuration.
This section provides a detailed description of each step of the core MicroHapulator workflow, as well as some operations that are supplementary to the primary workflow.

The core workflow consists of steps related to data preprocessing, haplotype calling, and quality control.

- The preprocessing steps merge each read pair into a single sequence and then align the merged read to MH reference sequences.
- The haplotype calling step examines each read alignment one at a time to genotype and phase the targeted SNPs and tally observed haplotypes. A filtering step is also performed to ignore low-abundance haplotypes likely resulting from sequencing error.
- Various quality control steps can be performed to assess factors that may impact or compromise accurate haplotype calling or forensic interpretation, such as interlocus balance and heterozygote balance.

MicroHapulator also provides tools for basic interpretation which are not included in the core workflow.

As with the end-to-end workflow, it is important to configure MicroHapulator correctly.
Be sure to read [the config docs](config.md) and prepare the three configuration files before beginning the workflow.

### Preprocessing: read merging

Paired-end NGS reads should be merged before they are aligned to the reference sequences.
The [FLASH program](https://ccb.jhu.edu/software/FLASH/) is recommended for this task.
(Run `flash --help` for instructions on configuring FLASH.)

```
flash EVD1-reads-R1.fastq.gz EVD1-reads-R2.fastq.gz --output-prefix EVD1 --allow-outies --threads 8
```

This step can be skipped for single-end reads.

### Preprocessing: read alignment

Merged reads must then be aligned to marker reference sequences before haplotypes can be called.
In this manual we use the [bwa mem](http://bio-bwa.sourceforge.net/bwa.shtml) algorithm, but other algorithms such as [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) would also be appropriate to use here.

First, the marker reference sequences must be indexed to allow NGS read alignment.
(Note, this only needs to be done once for each reference sequence file.)

```
bwa index marker-refr.fasta
```

The we align, sort, and index the merged reads.

```
bwa mem marker-refr.fasta EVD1.extendedFrags.fastq | samtools view -bS | samtools sort -o EVD1-reads.bam -
samtools index EVD1-reads.bam
```

### Haplotype calling

With NGS reads merged and aligned to the reference, we can proceed with haplotype calling.

```
mhpl8r type marker-defn.tsv EVD1-reads.bam --out EVD1-result-raw.json
```

The initial typing result will likely include numerous erroneous haplotypes due to sequencing errors.
It is typical to apply a combination of static and dynamic abundance (in this case, read count) thresholds to distinguish between true and false alleles.

```
mhpl8r filter EVD1-result-raw.json --static 10 --dynamic 0.02 --out EVD1-result.json
```

### Quality control

Interlocus balance looks at whether there are any markers with a disproportionately high or low number of aligned reads.
This variation could be due to any combination of factors, such as primer kinetics, off-target amplification, or stochastic effects in sequencing.
With the application of appropriate thresholds, interlocus imbalance shouldn't cause serious problems with forensic interpretation of a sample except in cases of extreme imbalance or the presence of DNA contributor(s) at very low levels in a sample.

Heterozygote balance looks only at markers with two alleles, and compares the relative abundance of the major and minor alleles.
Large differences in abundance between major and minor alleles can be a source of allelic drop-out, and should be accounted for in interpretation.

The following commands can be used to plot histograms of inter<u>loc</u>us and <u>het</u>erozygote balance.

```
mhpl8r locbalance EVD1-result.json --figure EVD1-loc-balance.png --title EVD1
mhpl8r hetbalance EVD1-result.json --figure EVD1-het-balance.png --title EVD1
```



### Analysis and interpretation

MicroHapulator provides some tools for analysis and basic interpretation of typing results.

```bash
# Estimate minimum number of DNA contributors in a sample
mhpl8r contrib --json EVD1-result.json

# Compute the random match probability (RMP) of a single-source sample
mhpl8r prob frequencies.tsv EVD1-result.json

# Show the differences between two profiles
mhpl8r diff EVD1-result.json EVD2-result.json

# Perform an RMP-based likelihood ratio (LR) test
mhpl8r prob frequencies.tsv EVD1-result.json REF1-result.json

# Perform a mixture containment query
mhpl8r contain EVD1-result.json REF2-result.json
```

## Appendix: Simulated Data

In addition to tools for haplotype calling, analysis, and interpretation, MicroHapulator also implements several tools for simulating genetic profiles and Illumina sequencing of those profiles.
These were used extensively in the early stages of MicroHapulator's development, and may be useful for testing more broadly.

```
# Simulate a mock profile—no haplotype counts, just genotype
mhpl8r sim frequencies.tsv --out mock1.json

# Combine mock profiles into a "mixture" profile
mhpl8r mix mock1.json mock2.json mock3.json --out mixture1.json

# Simulate reproduction to create "related" profiles
mhpl8r unite mom.json dad.json --out kid.json

# Simulate Illumina MiSeq sequencing of a mock profile
mhpl8r seq mypanel-def.tsv panel.fasta mock1.json --out mock1-reads.fastq --num-reads 250000
```
