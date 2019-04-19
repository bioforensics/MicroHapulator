#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from collections import defaultdict
from io import StringIO
import json
import microhapdb
import microhapulator


class SimulatedGenotype(object):
    """Genotype represented by phased alleles at a number of specified loci.

    Alleles are stored in a dictionary, with microhap locus ID/name as the key
    and a list as the value. Each list contains 2 items, the haplotype/phase 0
    allele and the hap/phase 1 allele.

    >>> gt = SimulatedGenotype()
    >>> gt.add(0, 'mh21KK-315', 'G,C,T')
    >>> gt.add(1, 'mh21KK-315', 'A,T,C')
    >>> gt.add(0, 'mh21KK-316', 'A,C,G,T')
    >>> gt.add(1, 'mh21KK-316', 'A,T,G,C')
    >>> print(gt)
    mh21KK-315 102     103     G|A
    mh21KK-315 207     208     C|T
    mh21KK-315 247     248     T|C
    mh21KK-316 108     109     A|A
    mh21KK-316 132     133     C|T
    mh21KK-316 179     180     G|G
    mh21KK-316 242     243     T|C
    """
    def __init__(self, frombed=None):
        self._data = defaultdict(lambda: [None] * 2)
        self._contexts = dict()
        if frombed:
            self._populate(frombed)

    def _populate(self, bedstream):
        locus_alleles = defaultdict(lambda: [list(), list()])
        for line in bedstream:
            line = line.strip()
            if line == '':
                continue
            locusid, start, end, allelestr = line.split('\t')
            alleles = allelestr.split('|')
            assert len(alleles) == 2
            for i, a in enumerate(alleles):
                locus_alleles[locusid][i].append(a)
        for locusid, allele_list in locus_alleles.items():
            for i, allele in enumerate(allele_list):
                self.add(i, locusid, ','.join(allele))

    def add(self, hapid, locusid, allele):
        assert hapid in (0, 1)
        self._data[locusid][hapid] = allele
        if locusid not in self._contexts:
            locus = microhapdb.id_xref(locusid).iloc[0]
            context = microhapulator.panel.LocusContext(locus)
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

    def __eq__(self, other):
        if type(other) != type(self) and type(other) != ObservedGenotype:
            return False
        if type(other) == ObservedGenotype:
            if self._data.keys() != other.data.keys():
                return False
            for locusid in self._data:
                if set(self._data[locusid]) != set(other.data[locusid]['genotype']):
                    return False
            return True
        else:
            return self._data == other._data


class ObservedGenotype(object):
    def __init__(self, filename=None):
        self.data = dict()
        if filename:
            with microhapulator.open(filename, 'r') as fh:
                self.data = json.load(fh)

    def record_coverage(self, locusid, cov_by_pos, ndiscarded=0):
        self.data[locusid] = {
            'mean_coverage': round(sum(cov_by_pos) / len(cov_by_pos), 1),
            'min_coverage': min(cov_by_pos),
            'max_coverage': max(cov_by_pos),
            'num_discarded_reads': ndiscarded,
            'allele_counts': dict(),
        }

    def record_allele(self, locusid, allele, count):
        self.data[locusid]['allele_counts'][allele] = count

    def infer(self):
        for locusid, locusdata in self.data.items():
            allelecounts = locusdata['allele_counts']
            avgcount = sum(allelecounts.values()) / len(allelecounts.values())
            gt = set()
            for allele, count in allelecounts.items():
                if count * 4 < avgcount:
                    continue
                gt.add(allele)
            self.data[locusid]['genotype'] = sorted(gt)

    def dump(self, file=None):
        if file is None:
            return json.dumps(self.data, indent=4, sort_keys=True)
        else:
            return json.dump(self.data, file, indent=4, sort_keys=True)

    def alleles(self):
        a = dict()
        for locusid in sorted(self.data):
            a[locusid] = ':'.join(self.data[locusid]['genotype'])
        return a
