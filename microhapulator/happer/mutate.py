# -------------------------------------------------------------------------------------------------
# Copyright (c) 2018, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from . import seqio
from .allele import parse_alleles
from .mutablestring import MutableString
import re


class PloidyMismatchError(ValueError):
    pass


def populate_haplotype_index(allelestream):
    """Parse a BED file to populate an index of haplotypes.

    Returns a tuple containing the ploidy and the haplotype index.

    The index is a dictionary mapping each sequence ID to a list of haplotypes
    for the corresponding sequence. Each list then contains nested lists of
    alleles representing each haplotype. The nested lists are sorted by genomic
    position.

    haplotypes = {
        'chr1': [
            [AlleleObj1a, AlleleObj1b, ..., AlleleObj1N],  # haplotype 1
            [AlleleObj2a, AlleleObj2b, ..., AlleleObj2N],  # haplotype 2
        ],
        'chr2': [
            [AlleleObj3a, AlleleObj3b, ..., AlleleObj3N],  # haplotype 1
            [AlleleObj4a, AlleleObj3b, ..., AlleleObj3N],  # haplotype 2
        ],
        ...
    }
    """
    n = None
    haplotypes = dict()
    for genotype in allelestream:
        seqid = genotype[0].seqid
        n = len(genotype)
        if seqid not in haplotypes:
            haplotypes[seqid] = list()
            while len(haplotypes[seqid]) < n:
                haplotypes[seqid].append(list())
        else:
            testn = len(haplotypes[seqid])
            if n != testn:
                message = "ploidy confusion: {:d} vs {:d}".format(n, testn)
                raise PloidyMismatchError(message)
        for index, allele in zip(haplotypes[seqid], genotype):
            index.append(allele)
    for seqid in haplotypes:
        for hap in haplotypes[seqid]:
            hap.sort()
    if n is None:
        raise ValueError("no input provided to populate haplotype index")
    return n, haplotypes


def mutate(seqstream, alleles):
    ploidy, haplotypes = populate_haplotype_index(parse_alleles(alleles))

    for defline, sequence in seqstream:
        haploseqs = list()
        while len(haploseqs) < ploidy:
            haploseqs.append(MutableString(sequence))
        seqid = defline.strip().split()[0]
        if seqid in haplotypes:
            seqhaps = haplotypes[seqid]
            for sequence, allelelist in zip(haploseqs, seqhaps):
                for allele in reversed(allelelist):
                    sequence[allele.start : allele.end] = allele.seq
        for n, hs in enumerate(haploseqs, 1):
            newdefline = re.sub(r"^(\S+)", r"\1:hap{}".format(n), defline)
            yield newdefline, str(hs)


def main(args):
    seq = open(args.seqfile, "r")
    als = open(args.bed, "r")
    out = open(args.out, "w")
    seqstream = seqio.parse_fasta(seq)
    for defline, haploseq in mutate(seqstream, als):
        print(">", defline, sep="", file=out)
        seqio.format(haploseq, out)
    out.close()
