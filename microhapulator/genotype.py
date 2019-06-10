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


def load_schema():
    with microhapulator.open(microhapulator.package_file('genotype-schema.json'), 'r') as fh:
        return json.load(fh)


schema = None


class Genotype(object):
    def __init__(self, fromfile=None):
        global schema
        if fromfile:
            if isinstance(fromfile, str):
                with microhapulator.open(fromfile, 'r') as fh:
                    self.data = json.load(fh)
            else:
                self.data = json.load(fromfile)
            if schema is None:
                schema = load_schema()
            jsonschema.validate(instance=self.data, schema=schema)
        else:
            self.data = self.initialize()

    @property
    def gttype(self):
        return 'base'

    def initialize(self):
        data = {
            'version': microhapulator.__version__,
            'type': self.gttype,
            'ploidy': None,
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

    def rand_match_prob(self, popid):
        """Compute the random match probability of this genotype.

        Given a set of population allele frequencies, the random match
        probability is the product of the allele frequencies of each allele
        observed in the genotype.
        """
        prob = 1.0
        for locus in sorted(self.loci()):
            alleles = self.alleles(locus)
            alleles = [*alleles] * 2 if len(alleles) == 1 else alleles
            for allele in sorted(alleles):
                result = microhapdb.frequencies[
                    (microhapdb.frequencies.Population == popid) &
                    (microhapdb.frequencies.Locus == locus) &
                    (microhapdb.frequencies.Allele == allele)
                ]
                frequency = None
                if len(result) > 0:
                    assert len(result) == 1
                    frequency = list(result.Frequency)[0]
                if frequency is None or frequency == 0.0:
                    if frequency is None:
                        message = 'No population allele frequency data for '
                    else:
                        message = 'Frequency=0.0 for '
                    message += 'allele "{:s}" at locus "{:s}" '.format(allele, locus)
                    message += 'for population "{:s}"; '.format(popid)
                    message += 'using RMP=0.001 for this allele'
                    microhapulator.plog('[MicroHapulator::genotype] WARNING:', message)
                    prob *= 0.001
                else:
                    prob *= frequency
        return prob

    def rmp_lr_test(self, other, popid, erate=0.001):
        """Compute a likelihood ratio test for the random match probability.

        The likelihood ratio test compares the probability that the two samples
        come from the same source to the probability that the samples come from
        two unrelated sources. The first probability (numerator) should be 1.0,
        with any missing or incongruent alleles treated as genotyping error.
        The second probability (denominator) is the random match probability.
        The test only makes sense when the two genotypes being compared are
        identical or nearly identical.
        """
        assert self.data['ploidy'] == 2 and other.data['ploidy'] == 2
        mismatches = 0
        for locus in self.loci():
            selfalleles = self.alleles(locus)
            otheralleles = other.alleles(locus)
            if selfalleles == otheralleles:
                pass
            elif len(selfalleles & otheralleles) == 1:
                mismatches += 1
            else:
                mismatches += 2
        numerator = erate ** mismatches
        denominator = self.rand_match_prob(popid)
        return numerator / denominator

    def dump(self, outfile):
        if isinstance(outfile, str):
            with microhapulator.open(outfile, 'w') as fh:
                json.dump(self.data, fh, indent=4, sort_keys=True)
        else:
            json.dump(self.data, outfile, indent=4, sort_keys=True)

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
        if self.data['ploidy'] is not None and self.data['ploidy'] > 0:
            assert hapid in range(self.data['ploidy'])
        if locusid not in self.data['loci']:
            self.data['loci'][locusid] = {'genotype': list()}
        self.data['loci'][locusid]['genotype'].append({'allele': allele, 'haplotype': hapid})

    def seqstream(self, seqindex, prechr=False):
        for locusid in sorted(self.loci()):
            canonid = microhapdb.id_xref(locusid).iloc[0]
            context = microhapulator.panel.LocusContext(canonid)
            yield context.defline(), context.sequence(seqindex, prechr=prechr)

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
