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

Under construction.

```
flash evid-repl1-R1.fastq evid-repl1-R2.fastq --output-prefix evid-repl1 --allow-outies --threads 8
bwa index panel.fasta
bwa mem panel.fasta evid-repl1.extendedFrags.fastq | samtools view -bS | samtools sort -o evid-repl1.bam -
samtools index evid-repl1.bam
```

### Haplotype calling

Under construction.

```
mhpl8r type mypanel-def.tsv evid-repl1.bam --out evid-repl1-tallies.json
mhpl8r type mypanel-def.tsv evid-repl1.bam --out evid-repl1-typing-result.json --effcov 0.25 --static 10 --dynamic 0.8
```

### Analysis and interpretation

Under construction.

```
mhpl8r balance evid-repl1-tallies.json
mhpl8r contrib --json evid-repl1-typing-result.json
mhpl8r prob frequencies.tsv evid-repl1-typing-result.json
mhpl8r diff evid-repl1-typing-result.json evid-repl2-typing-result.json
mhpl8r dist evid-repl1-typing-result.json evid-repl2-typing-result.json
mhpl8r prob frequencies.tsv evid-repl1-typing-result.json reference-typing-result.json
```

## Appendix: simulated data

In addition to tools for haplotype calling, analysis, and interpretation, MicroHapulator also implements several tools for simulating genetic profiles and Illumina sequencing of those profiles.

Under construction.

```
mhpl8r sim frequencies.tsv --out mock1.json
mhpl8r mix mock1.json mock2.json mock3.json --out mixture1.json
mhpl8r unite mom.json dad.json --out kid.json
mhpl8r seq mypanel-def.tsv panel.fasta mock1.json --out mock1-reads.fastq --num-reads 250000
```


[^f]: Less commonly, short insertion/deletion polymorphisms (INDELs) are used to define a microhaplotype marker, but these markers are incompatible with MicroHapulator.
