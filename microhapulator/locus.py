#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


from numpy.random import choice
import microhapdb


def default_panel():
    loci = microhapdb.loci.query('Source == "ALFRED"').\
        sort_values('AvgAe', ascending=False).\
        drop_duplicates('Chrom')
    return list(loci.ID)


def validate_loci(panel):
    valid_loci = microhapdb.standardize_ids(panel)
    if len(valid_loci) < len(panel):
        message = 'panel includes duplicate and/or invalid locus IDs'
        microhapdb.plog('[MicroHapulator::loci] WARNING', message)
    return valid_loci


def sample_panel(popids, loci):
    for haplotype, popid in enumerate(popids):
        for locusid in loci:
            f = microhapdb.frequencies
            allelefreqs = f[(f.Population == popid) & (f.Locus == locusid)]
            if len(allelefreqs) == 0:
                message = 'no allele frequencies available'
                message += ' for population "{pop}"'.format(pop=popid)
                message += ' at locus "{loc}"'.format(loc=locusid)
                message += '; in "relaxed" mode, drawing an allele uniformly'
                microhapdb.plog('[MicroHapulator::loci] WARNING:', message)
                alleles = list(f[f.Locus == locusid].Allele.unique())
                sampled_allele = choice(alleles)
            else:
                alleles = list(allelefreqs.Allele)
                freqs = list(allelefreqs.Frequency)
                freqs = [x / sum(freqs) for x in freqs]
                sampled_allele = choice(alleles, p=freqs)
            yield haplotype, locusid, sampled_allele
