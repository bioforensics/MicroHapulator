# Configuration

MicroHapulator's core operations depend on three sets of metadata describing the panel used for microhaplotype sequencing.
More information for each is provided below.
Not every piece of metadata is used for each operation, but it is recommended that all three be prepared at the same time, prior to any analysis or interpretation.

- microhaplotype marker reference sequences
- microhaplotype marker definitions
- population haplotype frequencies


## Reference sequences

A reference sequence must be provided for each marker, and aggregated into a single FASTA file.
At a minimum, each reference sequence must contain the SNPs associated with the microhaplotype, but additional flanking sequence can also be included.

The [MicroHapDB](https://github.com/bioforensics/MicroHapDB) database contains a comprehensive collection of published microhaplotypes and can be used to retrieve reference sequences for published or custom panels.
The example below shows how to use MicroHapDB retrieve reference sequences for a (trivially small) three-marker panel.

```
$ microhapdb marker --format=fasta --delta=25 --min-length=200 mh02USC-2pA mh08USC-8qA mh17USC-17qA
>mh02USC-2pA PermID=MHDBM-1734fe04 GRCh38:chr2:10810991-10811069 variants=61,105,112,139
TGTGAGCAGATGGCTCTGCCATCTGTGGTTTTCCATTCCTCCATGGCAGTTAAGGATGCTATAAATTCACATCAGTAACA
TCTCTGTCCCACATCTGCTAGGAGGAGTCTACAGTCCCAGAGCAGTCATCTAGAAACCATGAGGTTTCCCAGAAAGGCTG
CCATTTATTAAAATGAGTTAATTAATAGACTCACGGACTCT
>mh08USC-8qA PermID=MHDBM-d80aad76 GRCh38:chr8:62032395-62032463 variants=66,96,134
AGCCACTTGCAATGGTTATCTTGAAAATTGGAGGAATGTCCTCAGGATCCATTGAGTAGCAACAACATGAATGATAATTT
TTGCATCCCAAGTTCAGAAGAAGCAACCGACATTGAGTTTTTCTTGAATCCTTCTAAGAGAGACCTTCAGATAAATATTC
CAGACTGGCTTATCTCCTCTCTATCTCATCTATAATTTTGT
>mh17USC-17qA PermID=MHDBM-cd7a9041 GRCh38:chr17:27762200-27762287 variants=56,60,94,143
CCTGCTCTGCTTCCCTCAGTGACCCTGACCTTTCCTTCTTGGCTTCAGCACAATGAAAATCATCTAATTATTTCTATGCC
TATCTCTCCACTGCCGCCTCTAGATGCTCCCTGAGGGCAGGGCTGATCTCTGTTGTATTTCCATACCCACCCCAGTGCCA
GGCACCCACTGAATAGCTGAGTGAGTGAATGAATAATCAA
```

For larger / full-sized panels, you may want to place the marker names one per line in a text file (e.g. `mypanel.txt`) and write the reference sequences directly to a FASTA file (e.g. `refrseqs.fasta`) like so.

```
$ microhapdb marker --format=fasta --delta=25 --min-length=200 --panel=mypanel.txt > refrseqs.fasta
$ # Count the number of entries in the FASTA file as a sanity check
$ grep -c '^>' refrseqs.fasta
26
```


## Marker definitions

MicroHapulator needs to know the location of every SNP of interest in the corresponding reference sequence.
The list of SNP positions for a microhaplotype is its *marker definition*, and this information is provided in a tab-separated tabular plain text (TSV) file.
The **Marker** column contains the name/label/designator of a microhaplotype in the panel, and the **Offset** column contains the distance of one SNP from the beginning of the reference sequence.
For example, if a SNP of interest is the very first nucleotide in the reference, it has a distance of 0 from the beginning of the sequence and thus its **Offset** is `0`.
If a SNP is the 10th nucleotide, its offset is `9`.

[MicroHapDB](https://github.com/bioforensics/MicroHapDB) (version 0.7 or greater) can also be used to prepare marker definition files.
Simply change the `--format=fasta` setting to `--format=offsets`.

```
$ microhapdb marker --format=offsets --delta=25 --min-length=200 mh02USC-2pA mh08USC-8qA mh17USC-17qA
Marker	Offset
mh02USC-2pA	61
mh02USC-2pA	105
mh02USC-2pA	112
mh02USC-2pA	139
mh08USC-8qA	66
mh08USC-8qA	96
mh08USC-8qA	134
mh17USC-17qA	56
mh17USC-17qA	60
mh17USC-17qA	94
mh17USC-17qA	143
```

As before, you probably want to write the result directly to a TSV file.

```
$ microhapdb marker --format=offsets --delta=25 --min-length=200 --panel=mypanel.txt > panel-def.tsv
```


## Population frequencies

Simulating mock profiles and performing forensic interpretation depends on reliable estimates of population microhaplotype frequencies.
These are also provied as a tab-separated tabular plain text (TSV) file.
The **Marker** column contains the name/label/designator of a microhaplotype in the panel, the **Haplotype** column contains a comma-separated list of SNP alleles, and the **Frequency** column contains the relative prevalance of that haplotype in the population of interest.

[MicroHapDB](https://github.com/bioforensics/MicroHapDB) contains population frequency estimates for a comprehensive set of published microhaplotype markers, including 26 global populations from the [1000 Genomes Project](https://www.internationalgenome.org/).
MicroHapDB (version 0.7 or greater) can format this frequency data for use with MicroHapulator.
Using the correct population identifier (running `microhapdb population` beforehand if needed), haplotype frequencies can be retrieved and formatted as follows.

```
$ microhapdb frequency --format mhpl8r --population PUR --marker mh02USC-2pA mh08USC-8qA mh17USC-17qA
Marker	Haplotype	Frequency
mh02USC-2pA	A,A,G,A	0.038
mh02USC-2pA	T,A,A,T	0.428
mh02USC-2pA	T,A,G,A	0.144
mh02USC-2pA	T,A,G,T	0.125
mh02USC-2pA	T,C,G,A	0.255
mh02USC-2pA	T,C,G,T	0.01
mh08USC-8qA	A,A,C	0.154
mh08USC-8qA	A,A,T	0.308
mh08USC-8qA	A,G,T	0.351
mh08USC-8qA	C,A,C	0.188
mh17USC-17qA	A,C,C,T	0.26
mh17USC-17qA	A,C,T,T	0.245
mh17USC-17qA	G,A,T,T	0.197
mh17USC-17qA	G,C,C,C	0.298
```

And of course, it's always a good idea to write the results directly to a file.

```
$ microhapdb frequency --format=mhpl8r --population=PUR --panel=mypanel.txt > freqs.tsv
```

## Footnote

If the marker and/or frequency data data for your panel is not available in MicroHapDB, you should be able prepare the files manually, using the examples above as a reference for formatting.
In any case, it should go without saying, but just to be clear: a microhaplotype must use a single name/label/designator across *all* configuration files.
Using one name for the microhaplotype in the marker definition file and a different name for the microhaplotype in the population frequency file or reference sequence file will lead to issues.
