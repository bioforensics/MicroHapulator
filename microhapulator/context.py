#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from math import ceil
import microhapdb
import microhapulator


class LocusContext(object):
    """A class for resolving the context around a microhaplotype locus.

    >>> row = microhapdb.id_xref('mh01KK-172').iloc[0]
    >>> context = microhapulator.LocusContext(row)
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
        return '{loc} GRCh38:{chrom}:{s}-{e}'.format(
            loc=self._data['ID'], chrom=self.chrom, s=self.start,
            e=self.end
        )
