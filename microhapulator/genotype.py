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
import jsonschema
import microhapdb
import microhapulator


# with microhapulator.open(microhapulator.package_file('genotype-schema.json'), 'r') as fh:
#     schema = json.load(fh)


class Genotype(object):
    def __init__(self, fromfile=None):
        if fromfile:
            if isinstance(fromfile, str):
                with microhapulator.open(fromfile, 'r') as fh:
                    self.data = json.load(fh)
            else:
                self.data = json.load(fromfile)
#            jsonschema.validate(instance=self.data, schema=schema)
        else:
            self.data = self.initialize()

    @property
    def gttype(self):
        return 'base'

    def seqstream(self, seqindex, prechr=False):
        for locusid in sorted(self.loci()):
            canonid = microhapdb.id_xref(locusid).iloc[0]
            context = microhapulator.panel.LocusContext(canonid)
            yield context.defline(), context.sequence(seqindex, prechr=prechr)

    def initialize(self):
        data = {
            'version': microhapulator.__version__,
            'type': self.gttype,
            'ploidy': 0,
            'loci': dict(),
        }
        return data

    def haplotypes(self):
        hapids = set()
        for locusid, locusdata in self.data['loci'].items():
            for alleledata in locusdata['genotype']:
                if 'haplotype' in alleledata:
                    hapids.add(alleledata['haplotype'])
        assert sorted(hapids) == sorted(range(len(hapids)))
        return hapids

    def loci(self):
        return set(list(self.data['loci']))

    def alleles(self, locusid, haplotype=None):
        if locusid not in self.data['loci']:
            return None
        if haplotype is not None:
            return set(
                [a['allele'] for a in self.data['loci'][locusid]['genotype']
                 if 'haplotype' in a and a['haplotype'] == haplotype]
            )
        return set([a['allele'] for a in self.data['loci'][locusid]['genotype']])

    def dump(self, filename):
        with microhapulator.open(filename, 'w') as fh:
            json.dump(self.data, fh, indent=4, sort_keys=True)

    def __eq__(self, other):
        if not isinstance(other, Genotype):
            return False
        if self.loci() != other.loci():
            return False
        for locusid in self.loci():
            if self.alleles(locusid) != other.alleles(locusid):
                return False
        return True

    def __str__(self):
        return json.dumps(self.data, indent=4, sort_keys=True)


class SimulatedGenotype(Genotype):
    """Genotype represented by phased alleles at a number of specified loci.

    Alleles are stored in a dictionary, with microhap locus ID/name as the key
    and a list as the value. Typically each list contains 2 items, the
    haplotype/phase 0 allele and the hap/phase 1 allele. However, if the
    genotype represents a mixture there may be more than 2 haplotypes, as
    indicated by the `ploidy` parameter.

    >>> gt = SimulatedGenotype()
    >>> gt.add(0, 'mh21KK-315', 'G,C,T')
    >>> gt.add(1, 'mh21KK-315', 'A,T,C')
    >>> gt.add(0, 'mh21KK-316', 'A,C,G,T')
    >>> gt.add(1, 'mh21KK-316', 'A,T,G,C')
    >>> print(gt.bedstr)
    mh21KK-315 102     103     G|A
    mh21KK-315 207     208     C|T
    mh21KK-315 247     248     T|C
    mh21KK-316 108     109     A|A
    mh21KK-316 132     133     C|T
    mh21KK-316 179     180     G|G
    mh21KK-316 242     243     T|C
    """
    def populate_from_bed(bedfile):
        with microhapulator.open(bedfile, 'r') as fh:
            line = next(fh)
            ploidy = line.count('|') + 1
            fh.seek(0)
            locus_alleles = defaultdict(lambda: [list() for _ in range(ploidy)])
            for line in fh:
                line = line.strip()
                if line == '':
                    continue
                locusid, start, end, allelestr = line.split('\t')
                alleles = allelestr.split('|')
                for i, a in enumerate(alleles):
                    locus_alleles[locusid][i].append(a)
            genotype = SimulatedGenotype(ploidy=ploidy)
            for locusid, allele_list in locus_alleles.items():
                for i, allele in enumerate(allele_list):
                    genotype.add(i, locusid, ','.join(allele))
            return genotype

    def merge(genotypes):
        ploidy = 2 * len(genotypes)
        gt = SimulatedGenotype(ploidy=ploidy)
        offset = 0
        for genotype in genotypes:
            for locusid, locusdata in sorted(genotype.data['loci'].items()):
                for allele in locusdata['genotype']:
                    gt.add(offset + allele['haplotype'], locusid, allele['allele'])
            offset += 2
        return gt

    def __init__(self, fromfile=None, ploidy=0):
        super(SimulatedGenotype, self).__init__(fromfile=fromfile)
        if ploidy:
            self.data['ploidy'] = ploidy

    def add(self, hapid, locusid, allele):
        if self.data['ploidy'] > 0:
            assert hapid in range(self.data['ploidy'])
        if locusid not in self.data['loci']:
            self.data['loci'][locusid] = {'genotype': list()}
        self.data['loci'][locusid]['genotype'].append({'allele': allele, 'haplotype': hapid})

    @property
    def gttype(self):
        return 'SimulatedGenotype'

    @property
    def bedstream(self):
        hapids = self.haplotypes()
        for locusid in sorted(self.loci()):
            locusdata = microhapdb.id_xref(locusid).iloc[0]
            context = microhapulator.panel.LocusContext(locusdata)
            coords = microhapdb.allele_positions(locusdata['ID'])
            coords = list(map(context.global_to_local, coords))
            locusvars = [list() for _ in range(len(coords))]
            for haplotype in sorted(hapids):
                allele = self.alleles(locusid, haplotype=haplotype).pop()
                for var, varlist in zip(allele.split(','), locusvars):
                    varlist.append(var)
            for coord, var in zip(coords, locusvars):
                allelestr = '|'.join(var)
                yield '\t'.join((locusid, str(coord), str(coord + 1), allelestr))

    @property
    def bedstr(self):
        out = StringIO()
        for line in self.bedstream:
            print(line, file=out)
        return out.getvalue()


class ObservedGenotype(Genotype):
    def __init__(self, fromfile=None):
        super(ObservedGenotype, self).__init__(fromfile=fromfile)

    def record_coverage(self, locusid, cov_by_pos, ndiscarded=0):
        self.data['loci'][locusid] = {
            'mean_coverage': round(sum(cov_by_pos) / len(cov_by_pos), 1),
            'min_coverage': min(cov_by_pos),
            'max_coverage': max(cov_by_pos),
            'num_discarded_reads': ndiscarded,
            'allele_counts': dict(),
        }

    def record_allele(self, locusid, allele, count):
        self.data['loci'][locusid]['allele_counts'][allele] = count

    def infer(self, threshold=None):
        for locusid, locusdata in self.data['loci'].items():
            allelecounts = locusdata['allele_counts']
            avgcount = sum(allelecounts.values()) / len(allelecounts.values())
            gt = set()
            for allele, count in allelecounts.items():
                if threshold is not None:
                    if count < threshold:
                        continue
                else:
                    if count * 4 < avgcount:
                        continue
                gt.add(allele)
            self.data['loci'][locusid]['genotype'] = [
                {'allele': a, 'haplotype': None} for a in sorted(gt)
            ]
