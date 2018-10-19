#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from math import ceil
import microhapulator


class LocusContext(object):
    """A class for resolving the context around a microhaplotype locus.

    >>> row = microhapulator.bogus_loci[0]
    >>> context = LocusContext(row)
    >>> len(context)
    301
    >>> context.offset
    983
    >>> context.global_to_local(1111)
    128
    >>> context.local_to_global(101)
    1084
    """
    def __init__(self, rowdata, mindelta=30, minlen=300):
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
        return seqindex[seqid][self.start:self.end].seq

    def defline(self):
        return '{loc} {chrom}:{s}-{e}'.format(
            loc=self._data['Name'], chrom=self.chrom, s=self.start,
            e=self.end
        )
