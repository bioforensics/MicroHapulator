# Configuration

MicroHapulator's core operations depend on three sets of metadata describing the panel used for microhaplotype sequencing.
More information for each is provided below.
Not every piece of metadata is used for each operation, but it is recommended that all three be prepared at the same time, prior to any analysis or interpretation.

- microhaplotype marker reference sequences
- microhaplotype marker definitions
- population haplotype frequencies

This document will demonstrate how to prepare these files using the [MicroHapDB](https://github.com/bioforensics/MicroHapDB) database, which contains a comprehensive collection of published microhaplotypes.
A trivially small panel composed of three arbitrary markers (identifiers: mh01KK-205, mh03USC-3qC, mh18CP-005) will be used as an example to show how the data is formatted.
When preparing a full panel composed of dozens of markers, it is recommended that marker identifiers be placed in a plain text file, one identifier per line, and that the configuration data be written to directly to files rather than printed to the terminal screen.
This will also be demonstrated.


## Reference sequences

A reference sequence must be provided for each marker, and aggregated into a single FASTA file.
At a minimum, each reference sequence must contain the SNPs associated with the microhaplotype, but additional flanking sequence can also be included.
The following command demonstrates how to use MicroHapDB to retrieve these reference sequences.

```
$ microhapdb marker --format=fasta --delta=25 --min-length=200 mh01KK-205 mh03USC-3qC mh18CP-005
>mh01KK-205 PermID=MHDBM-1f7eaca2 GRCh38:chr1:18396197-18396351 variants=25,46,134,179 Xref=SI664550A
CACCAGTTCTCATGAATCTGAGGAATTCTTCCTCCTAGCTACTTCCTTCCTTTTCCCTCATTACATCCCTGCCAAGGACA
AATTCTGCCATTTGCATGGCAGGACTCCTCCAAAAAGGGGCTTCCTCCCTTTCCGTTAGTAAAGGAAGAGGTTACCTGAG
ACTTGACTTAACCTCCTTGGGAGGGAACATGCTTTCACTGTTGCG
>mh03USC-3qC PermID=MHDBM-eacabfd9 GRCh38:chr3:196653025-196653121 variants=52,61,71,111,148
TAGCATTGAAATGATGCCTTGTAATTTACTAAATCTGCAACTATGCAGCCTTATTTCATGGCGGGCAGTGGTGGTGATCC
CAGGTTTCAGGGGCGGGGAAGGGTGCTGGGGGGATCCTGAGGTCAGGAACCCGTACACCTCTGCTTCTGCCCTCTCTTCC
CTGTGCCGGCCACAAGGCAATGACTCCTGTGTGGGTGCAGA
>mh18CP-005 PermID=MHDBM-a85754d3 GRCh38:chr18:8892864-8892907 variants=78,107,110,121 Xref=SI664898P
GAGATTCTGTCTCAAAAAATAAAAAATTAAAAAAAATTTTTTTAAACCCAAAATATTACTGCAGATGTCCTTATACGCAG
TGGTGTTAGTTTTAGAAACTGATTCTACGGGTATGCTTGCTCGTGTGTAAAATTATTCATATACAAATTATTTATGACAG
TATTGTTTCTAGTAGTAAAATATCGGAAATATTCTAAATG
```


## Marker definitions

MicroHapulator needs to know the location of every SNP of interest in the corresponding reference sequence.
The list of SNP positions for a microhaplotype is its *marker definition*, and this information is provided in a tab-separated tabular plain text (TSV) file.
The **Marker** column contains the identifier (name, label, or designator) of a microhaplotype in the panel, and the **Offset** column contains the distance of one SNP from the beginning of the reference sequence.
For example, if a SNP of interest is the very first nucleotide in the reference, it has a distance of 0 from the beginning of the sequence and thus its offset is `0`.
If a SNP is the 10th nucleotide, its offset is `9`.
The **Chrom** and **OffsetHg38** columns indicate the position of each SNP in the GRCh38 reference human genome assembly.
**As of version MicroHapulator version 0.8, these two columns are now required.**

The following command shows how to use MicroHapDB (version 0.8 or greater) to prepare a marker definition file.
Note that it is identical to the previous command, except that the `--format=fasta` setting was changed to `--format=offsets`.

```
$ microhapdb marker --format=offsets --delta=25 --min-length=200 mh01KK-205 mh03USC-3qC mh18CP-005
Marker	Offset	Chrom	OffsetHg38
mh01KK-205	25	chr1	18396197
mh01KK-205	46	chr1	18396218
mh01KK-205	134	chr1	18396306
mh01KK-205	179	chr1	18396351
mh03USC-3qC	52	chr3	196653025
mh03USC-3qC	61	chr3	196653034
mh03USC-3qC	71	chr3	196653044
mh03USC-3qC	111	chr3	196653084
mh03USC-3qC	148	chr3	196653121
mh18CP-005	78	chr18	8892864
mh18CP-005	107	chr18	8892893
mh18CP-005	110	chr18	8892896
mh18CP-005	121	chr18	8892907
```


## Population frequencies

Performing forensic interpretation or simulating mock profiles depends on reliable estimates of population microhaplotype frequencies.
These must also be provided to MicroHapulator as a tab-separated tabular plain text (TSV) file.
The **Marker** column contains the name/label/designator of a microhaplotype in the panel, the **Haplotype** column contains a comma-separated list of SNP alleles, and the **Frequency** column contains the relative prevalance of that haplotype in the population of interest.

MicroHapDB contains population frequency estimates from 26 global populations in the [1000 Genomes Project](https://www.internationalgenome.org/) for most of its markers.
MicroHapDB (version 0.7 or greater) can format this frequency data for use with MicroHapulator.
Using the correct population identifier (running `microhapdb population` beforehand if needed), haplotype frequencies can be retrieved and formatted as follows.
(The "PUR" population, "Puerto Ricans from Puerto Rico", is used for this example.)

```
$ microhapdb frequency --format mhpl8r --population PUR --marker mh01KK-205 mh03USC-3qC mh18CP-005
Marker	Haplotype	Frequency
mh01KK-205	C,C,A,G	0.149
mh01KK-205	T,C,A,G	0.361
mh01KK-205	T,T,A,A	0.183
mh01KK-205	T,T,A,G	0.13
mh01KK-205	T,T,G,G	0.178
mh03USC-3qC	A,C,C,A,G	0.043
mh03USC-3qC	A,C,C,G,G	0.12
mh03USC-3qC	A,C,C,G,T	0.034
mh03USC-3qC	A,C,T,G,G	0.269
mh03USC-3qC	A,C,T,G,T	0.212
mh03USC-3qC	G,C,C,A,G	0.005
mh03USC-3qC	G,C,C,G,G	0.014
mh03USC-3qC	G,C,T,G,T	0.005
mh03USC-3qC	G,T,C,A,G	0.226
mh03USC-3qC	G,T,C,G,G	0.067
mh03USC-3qC	G,T,C,G,T	0.005
mh18CP-005	A,C,A,T	0.38
mh18CP-005	A,C,G,C	0.255
mh18CP-005	A,T,A,C	0.255
mh18CP-005	A,T,G,C	0.043
mh18CP-005	G,C,A,C	0.01
mh18CP-005	G,T,A,C	0.058
```


## Summary

As mentioned at the beginning of this document, you're much better off writing the config data directly to files rather than printing it to your screen.
If you have the marker identifiers for your panel (one identifier per line) in a plain text file named, say, `mypanel.txt`, you would create your config files like so.
(Replace "PUR" with the appropriate population.)

```
$ microhapdb marker --format=fasta --delta=25 --min-length=200 --panel=mypanel.txt > mypanel-refr.fasta
$ microhapdb marker --format=offsets --delta=25 --min-length=200 --panel=mypanel.txt > mypanel-defn.tsv
$ microhapdb frequency --format=mhpl8r --population=PUR --marker --panel=mypanel.txt > mypanel-freq-pr.tsv
```

These commands will create three files: `mypanel-refr.fasta` with the reference sequences, `mypanel-defn.tsv` with the marker definitions, and `mypanel-freq-pr.tsv` with the haplotype frequencies.


## Example configuration files

The [`microhapulator/data/configs/`](https://github.com/bioforensics/MicroHapulator/tree/master/microhapulator/data/configs/) directory in the MicroHapulator source code distribution contains example configuration files for a published panel.


## What if my data isn't in MicroHapDB?

If you have marker and/or frequency data that you would like to submit to MicroHapDB, that is always welcome!
See [this page](https://github.com/bioforensics/MicroHapDB#adding-markers-to-microhapdb) for details.

In any case, the examples above show how the reference sequences, marker definitions, and haplotype frequencies should be formatted.
So if your data is not included in MicroHapDB, you should still be able to configure MicroHapulator correctly, it will just take some extra time to prepare the files manually.
