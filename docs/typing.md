# Calling Haplotypes with MicroHapulator

MicroHapulator makes haplotype calls by aligning NGS reads to marker reference sequences and examining each aligned read.
If a given read spans all SNPs at the microhaplotype (MH) marker, the read reveals not only the allele of each individual SNP, but the *phase* of the alleles.
We call this aggregate allele grouping the *MH allele*.

Consider the following mock example.
The first line shows the reference sequence for the `mh01USC-1pD` MH marker.
Each subsequent line represents an NGS read aligned to the marker sequence.
The `*` symbols denote the locations of the SNPs that define the marker, and the `.` symbols denote locations where the NGS read matches the reference.

> *MicroHapulator considers only those SNPs that are explicitly indicated in the MH marker definition (see [the config docs](config.md) for details).
> Currently, it does not report or incorporate any additional rare SNPs that may be present at the locus.*

```
                       *                     *           *
AAATAGCTGGGCTAATAATGAACTGAAGCAAAGTCAACTGAAATGTCCTGGGCAGCTCCAGAAACTCCAGAATGGGGAGGA
.......................C.....................C.......
  .....................C.....................C...........A...
       .....A..........C.....................T...........C........
        ...............C.....................T...........C.........
          .............C.....................C...........A...........
            ...........C.....................C...........A.............
                 ......C.....................C...........A.................
                 ......C.....................C...........A.................
                    ...C.....................G...........A....................
                     ..C.....................T...........C.....................
```

The aligned reads allow us to determine the number of times each MH allele occurs.
We call this tally of MH allele counts a *typing result*.
The typing result for this example is as follows:

- the `C,C,A` haplotype is observed 5 times
- the `C,T,C` haplotype is observed 3 times
- the `C,G,A` haplotype is observed 1 time.

If we then filter out the `C,G,A` haplotype as erroneous (see below), we can infer a diploid `C,C,A` / `C,T,C` genotype for this marker.
We call this th *genotype call* or *genotype prediction*.

> It is helpful to point about a few observations about this example.
> - The first aligned read does not span all SNPs at the marker, so it is discarded.
> - The third aligned read shows an `A` at the 13th position of the marker reference sequence. Whether this reflects true genetic variation or is a technical artifact resulting from sequencing error, it is ignored by MicroHapulator because it is not one of the three SNPs of interest.
> - The `C,G,A` haplotype is only observed once and is likely a false haplotype resulting from sequencing error at the second SNP of interest in the haplotype. When dozens or hundreds (or thousands!) of reads are successfully sequenced and aligned to the marker reference, it is usually simple to distinguish signal (true haplotypes) from noise (false haplotypes resulting from sequencing error). When the depth of sequencing coverage is low, as it is in this example, it can be more difficult to distinguish signal from noise. Determining appropriate per-marker thresholds (detection thresholds and analytical thresholds) for filtering will typically require a non-trivial amount of testing with the laboratory's NGS sequencing instrument(s).

While this mock example is helpful in building intuition about haplotype calling, manual visual examination is not feasible for performing this task on dozens of markers and (potentially) millions of NGS reads.
The MicroHapulator software provides tools that automate haplotype calling and genotype prediction, as well as assist with basic interpretation of the forensic typing result.
