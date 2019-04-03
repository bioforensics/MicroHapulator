#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from collections import defaultdict
from io import StringIO
import microhapdb
import microhapulator


class Genotype(object):
    """Genotype represented by phased alleles at a number of specified loci.

    Alleles are stored in a dictionary, with microhap locus ID/name as the key
    and a list as the value. Each list contains 2 items, the haplotype/phase 0
    allele and the hap/phase 1 allele.

    >>> gt = microhapulator.Genotype()
    >>> gt.add(0, 'mh21KK-315', 'G,C,T')
    >>> gt.add(1, 'mh21KK-315', 'A,T,C')
    >>> gt.add(0, 'mh21KK-316', 'A,C,G,T')
    >>> gt.add(1, 'mh21KK-316', 'A,T,G,C')
    >>> print(gt)
    mh21KK-315	77	78	G|A
    mh21KK-315	182	183	C|T
    mh21KK-315	222	223	T|C
    mh21KK-316	83	84	A|A
    mh21KK-316	107	108	C|T
    mh21KK-316	154	155	G|G
    mh21KK-316	217	218	T|C
    """
    def __init__(self):
        self._data = defaultdict(lambda: [None] * 2)
        self._contexts = dict()

    def add(self, hapid, locusid, allele):
        assert hapid in (0, 1)
        self._data[locusid][hapid] = allele
        if locusid not in self._contexts:
            locus = microhapdb.id_xref(locusid).iloc[0]
            context = microhapulator.LocusContext(locus)
            self._contexts[locusid] = context

    def seqstream(self, seqindex, prechr=False):
        for locusid, context in sorted(self._contexts.items()):
            yield context.defline(), context.sequence(seqindex, prechr=prechr)

    @property
    def bedstream(self):
        for locusid in sorted(self._data):
            context = self._contexts[locusid]
            alleles_0 = self._data[locusid][0].split(',')
            alleles_1 = self._data[locusid][1].split(',')
            coords = microhapdb.allele_positions(locusid)
            for a0, a1, coord in zip(alleles_0, alleles_1, coords):
                localcoord = context.global_to_local(coord)
                allelestr = a0 + '|' + a1
                yield '\t'.join(
                    (locusid, str(localcoord), str(localcoord + 1), allelestr)
                )

    def __str__(self):
        out = StringIO()
        for line in self.bedstream:
            print(line, file=out)
        return out.getvalue()
