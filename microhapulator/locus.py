#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from math import ceil
import microhapdb
import microhapulator
from numpy.random import choice


class LocusContext(object):
    """A class for resolving the context around a microhaplotype locus.

    >>> row = microhapdb.id_xref('mh01KK-172').iloc[0]
    >>> context = LocusContext(row)
    >>> len(context)
    350
    >>> context.offset
    1551391
    >>> context.global_to_local(1551522)
    131
    >>> context.local_to_global(287)
    1551678
    """
    def __init__(self, rowdata, mindelta=30, minlen=350):
        self._data = rowdata
        start = rowdata['Start'] - mindelta
        end = rowdata['End'] + mindelta
        initlen = end - start
        if initlen < minlen:
            lendiff = minlen - initlen
            bpextend = ceil(lendiff / 2)
            start -= bpextend
            end += bpextend
        self.start = start
        self.end = end

    def __len__(self):
        return self.end - self.start

    @property
    def chrom(self):
        return str(self._data['Chrom'])

    @property
    def offset(self):
        return self.start

    def global_to_local(self, coord):
        if coord < self.start or coord > self.end:
            return None
        return coord - self.start

    def local_to_global(self, coord):
        if coord >= len(self):
            return None
        return coord + self.start

    def sequence(self, seqindex, prechr=False):
        seqid = 'chr' + self.chrom if prechr else self.chrom
        return str(seqindex[seqid][self.start:self.end].seq).upper()

    def defline(self):
        coords = [self.global_to_local(x) for x in microhapdb.allele_positions(self._data['ID'])]
        return '{loc} GRCh38:{chrom}:{s}-{e} variants={var}'.format(
            loc=self._data['ID'], chrom=self.chrom, s=self.start,
            e=self.end, var=':'.join(map(str, coords))
        )


def default_panel():
    loci = microhapdb.loci.query('Source == "ALFRED"').\
        sort_values('AvgAe', ascending=False).\
        drop_duplicates('Chrom')
    return list(loci.ID)


def validate_loci(panel):
    valid_loci = microhapdb.standardize_ids(panel)
    if len(valid_loci) < len(panel):
        message = 'panel includes duplicate and/or invalid locus IDs'
        microhapulator.plog('[MicroHapulator::loci] WARNING', message)
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
                microhapulator.plog('[MicroHapulator::loci] WARNING:', message)
                alleles = list(f[f.Locus == locusid].Allele.unique())
                sampled_allele = choice(alleles)
            else:
                alleles = list(allelefreqs.Allele)
                freqs = list(allelefreqs.Frequency)
                freqs = [x / sum(freqs) for x in freqs]
                sampled_allele = choice(alleles, p=freqs)
            yield haplotype, locusid, sampled_allele
