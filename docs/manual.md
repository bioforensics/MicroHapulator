# User Manual

## Overview

Microhaplotypes are versatile and effective genetic markers for human forensics, matching or exceeding the power of conventional STR markers for identification and mixture deconvolution without stutter artifacts that complicate analysis and forensic interpretation.
In the literature, a microhaplotype has been defined as a set of SNPs [^f] residing in close proximity in the genome, such that a single NGS read is capable of spanning the entire set.
This allows each read to not only genotype individual SNPs, but also to empirically phase the SNP alleles and determine the haplotype.
Existing microhap-based methods for human forensics have for the most part adopted SNP genotyping and statistical haplotype phasing procedures designed for genome-wide studies, procedures that do not leverage the unique nature of microhaplotypes.
MicroHapulator implements the first empirical haplotype caller designed specifically for microhaplotypes.

The MicroHapulator software is run by calling the `mhpl8r` command in a shell terminal.
Using the terminal can be intimidating for the uninitiated, so this manual is written to be acommodating of all users.
Forensic practicioners routinely perform tasks and analyses that are much more complex than MicroHapulator's interface, so we are confident that the command line interface will not be an obstacle to interested users.


## Processing workflow

The MicroHapulator data processing workflow is broken up roughly into three phases.

1. preprocessing
2. haplotype calling
3. analysis and interpretation

In the preprocessing phase, read pairs are merged and then aligned to microhaplotype reference sequences.
In the haplotype calling phase, read alignments are processed one at a time to genotype and phase the SNPs of interest and to tally observed haplotypes.
In the analysis and interpretation phase, various characteristics of the haplotype calls are examined to assess one or more forensic hypotheses.

Before executing the workflow, it is important to make sure to configure MicroHapulator correctly.
Be sure to read [the config docs](config.md) and prepare the three configuration files before beginning the workflow.

### Preprocessing

NGS reads from a sample must be aligned to marker reference sequences before haplotypes can be called.
In this manual we use the [bwa mem](http://bio-bwa.sourceforge.net/bwa.shtml) algorithm, but other algorithms such as [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) would also be appropriate to use here.

First, the marker reference sequences must be indexed to allow NGS read alignment.
(Note, this only needs to be done once for each reference sequence file.)

```
bwa index panel.fasta
```

Next, paired-end reads must be merged before they are aligned to the reference sequences.
The [FLASH program](https://ccb.jhu.edu/software/FLASH/) is recommended for this task.
(Run `flash --help` for instructions on configuring FLASH.)

```
flash EVD1-reads-R1.fastq.gz EVD1-reads-R2.fastq.gz --output-prefix EVD1 --allow-outies --threads 8
```

And finally with read pairs merged, we can align, sort, and index the NGS reads.

```
bwa mem panel.fasta EVD1.extendedFrags.fastq | samtools view -bS | samtools sort -o EVD1-reads.bam -
samtools index EVD1-reads.bam
```

### Haplotype calling

With NGS reads merged and aligned to the reference, we can proceed with haplotype calling.

```
mhpl8r type mypanel-def.tsv EVD1-reads.bam --out EVD1-result-raw.json
```

The initial typing result will likely include numerous erroneous haplotypes due to sequencing errors.
It is typical to apply a combination of static and dynamic abundance (in this case, read count) thresholds to distinguish between true and false alleles.

```
mhpl8r filter EVD1-result-raw.json --static 10 --dynamic 0.02 --out EVD1-result.json
```

### Analysis and interpretation

MicroHapulator provides some tools for analysis and basic interpretation of typing results.

```bash
# Examine interlocus balance
mhpl8r balance EVD1-result.json

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

## Appendix: simulated data

In addition to tools for haplotype calling, analysis, and interpretation, MicroHapulator also implements several tools for simulating genetic profiles and Illumina sequencing of those profiles.

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


[^f]: Less commonly, short insertion/deletion polymorphisms (INDELs) are used to define a microhaplotype marker, but these markers are incompatible with MicroHapulator.
