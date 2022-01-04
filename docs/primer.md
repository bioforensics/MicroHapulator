# Primer on Forensic DNA Typing

> *While MicroHapulator is designed only for analysis and interpretation of NGS-based microhaplotype assays, the following primer provides a basic introduction to the broader topic of forensic DNA analysis and basic related concepts.*

Forensic DNA typing has a variety of applications in anthropology, medicine, criminal justice, and investigations of missing persons.
The origin of a DNA sample collected for a study or investigation may be known (e.g. if it is collected from the victim of a crime) or unknown (e.g. if it is collected from body fluid present at a crime scene).
In either case, the principle behind DNA typing methods is straightforward: analyzing DNA to build a genetic profile that can then be compared to other profiles to confirm or refute hypothesized matches.

## Genetic profiles

A "genetic profile" here refers to the alleles present in a DNA sample as measured at specific well-studied genetic *markers*.
For instance, consider a hypothetical marker for which four alleles (A, B, C, and D) are known to be present in a population.
Also note that the *relative frequencies* of these alleles throughout the entire population is known: say, A=45%, B=11%, C=20%, and D=24%.
Then, if we analyze a DNA sample and determine that it has allele C at this marker, we can say that the probability of observing this DNA sample by chance is 0.2 (20%).
Forensic scientists call this the *random match probability*, or *RMP*.

The RMP associated with a single genetic marker is inconclusive: an allele present in 20% of the population cannot distinguish individuals on its own.
The power of genetic profiles is the confidence that comes from analyzing multiple markers simultaneously.
Continuing our hypothetical scenario, consider another marker (with allele frequences of X=60%, Y=8%, and Z=32%), and an analysis revealing that our DNA sample has allele Y at this marker.
If we look only at this second marker, our RMP of 0.08 is again inconclusive.
But if these two markers are genetically independent (the chance of inheriting one allele is unrelated to the chance of inheriting the other allele, as is the case for markers residing on different chromosomes) then the profile's overall RMP is the *product* of the per-marker RMPs: 0.2 × 0.08 = 0.016, a marginally more conclusive result.

Forensic DNA typing typically involves profiles composed of more than a dozen independent genetic markers, and RMPs computed for full profiles can in the best case reach 1 × 10<sup>-10</sup> or smaller.
Probabilities of this magnitude are more than suitable to make confident claims about the origin of matching DNA samples.


## Different types of genetic markers

### STRs

The genetic markers most commonly used for forensic DNA typing are *short tandem repeat polymorphisms*, or *STRs*.
These markers are accordian-like DNA sequences containing a target region with a short repeated sequence 3-5 nucleotides in length.
The number of consecutive repeat units can vary between individuals, and thus the allele reported for a STR marker typically indicates the number of repeat units present.
Consider the following STR marker, characterized by tandem `GATA` repeats.
Three alleles are shown: 7 (containing 7 `GATA` repeat units), 10, and 11.
(The `-` symbols are not part of the sequence and were only added for display purposes.)

```
ATACAGACAGAAGACAGGTG---GATA-GATA-GATA-GATA-GATA-GATA-GATA-----------------------TCATTGAAAGACAAAACAGA
ATACAGACAGAAGACAGGTG---GATA-GATA-GATA-GATA-GATA-GATA-GATA-GATA-GATA-GATA--------TCATTGAAAGACAAAACAGA
ATACAGACAGAAGACAGGTG---GATA-GATA-GATA-GATA-GATA-GATA-GATA-GATA-GATA-GATA-GATA---TCATTGAAAGACAAAACAGA
```

<small>(Example borrowed shamelessly from [this slide deck](https://www.nist.gov/system/files/documents/mml/bmd/genetics/03-Vallone-NGS-talk-STRs-final.pdf).)</small>

### SNPs

Another common type of genetic marker used in forensic DNA typing is the *single nucleotide polymorphism*, or *SNP*.
As its name suggests, a SNP is a single nucleotide that is present in the population with two or more alleles.
For example, the two DNA sequences below differ only at a single position, the top sequence containing a `C` and the bottom sequence containing a `G`.

```
GATTACAAAGA
     *
GATTAGAAAGA
```

Most known SNPs in the human genome are *biallelic*, meaning they assume one of two alleles.
But it's not uncommon for a SNP to be *multiallelic*, capable of assuming three or in some cases all four possible alleles.


### MHs

*Microhaplotypes*, or *MHs*, are a more recently proposed type of genetic marker.
MHs are composed of a number of SNPs residing in close proximity in the genome, such that they can be observed simultaneously by a single NGS sequence read (i.e., spanning roughly 300 nucleotides or fewer).
Each individual SNP in a MH is typically biallelic, but when considered together the SNPs determine a short (you could even say a "micro") haplotype that can be present in numerous allelic combinations in a population.
The example below shows a MH composed of three SNPs: the `*` symbols indicate the location of the SNPs, while the `.` symbols indicate DNA that is shared among all eight observed haplotypes.

```
                       *                     *           *
AAATAGCTGGGCTAATAATGAACTGAAGCAAAGTCAACTGAAATGTCCTGGGCAGCTCCAGAAACTCCAGAATGGGGAGGA
.......................C.....................C...........A.......................
.......................C.....................C...........C.......................
.......................C.....................T...........A.......................
.......................C.....................T...........C.......................
.......................T.....................C...........A.......................
.......................T.....................C...........C.......................
.......................T.....................T...........A.......................
.......................T.....................T...........C.......................
```
