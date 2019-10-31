#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from collections import defaultdict
from happer.mutate import mutate
from io import StringIO
import json
import jsonschema
import microhapdb
import microhapulator
from numpy.random import choice


def load_schema():
    with microhapulator.open(microhapulator.package_file('profile-schema.json'), 'r') as fh:
        return json.load(fh)


schema = None


class Profile(object):
    def unite(mom, dad):
        """Simulate the creation of a new profile from a mother and father."""
        gt = SimulatedProfile(ploidy=2)
        allloci = mom.loci() | dad.loci()
        commonloci = mom.loci() & dad.loci()
        if len(commonloci) == 0:
            raise ValueError('mom and dad profiles have no loci in common')
        notshared = allloci - commonloci
        if len(notshared) > 0:
            message = 'loci not common to mom and dad profiles are excluded: '
            message += ', '.join(notshared)
            microhapulator.plog('[MicroHapulator::profile]', message)
        for parent, hapid in zip((mom, dad), (0, 1)):
            for locus in sorted(commonloci):
                haploallele = choice(sorted(parent.alleles(locus)))
                gt.add(hapid, locus, haploallele)
        return gt

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

    @property
    def ploidy(self):
        return self.data['ploidy']

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

    def markers(self):
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
        """Compute the random match probability of this profile.

        Given a set of population allele frequencies, the random match
        probability is the product of the allele frequencies of each allele
        observed in the profile.
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
                    microhapulator.plog('[MicroHapulator::profile] WARNING:', message)
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
        The test only makes sense when the two profiles being compared are
        identical or nearly identical.
        """
        assert self.data['ploidy'] in (2, None)
        assert other.data['ploidy'] in (2, None)
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
        if not isinstance(other, Profile):
            return False
        if self.loci() != other.loci():
            return False
        for locusid in self.loci():
            if self.alleles(locusid) != other.alleles(locusid):
                return False
        return True

    def __str__(self):
        return json.dumps(self.data, indent=4, sort_keys=True)

    @property
    def bedstream(self):
        hapids = self.haplotypes()
        for locusid in sorted(self.loci()):
            result = microhapdb.markers[microhapdb.markers.Name == locusid]
            if len(result) == 0:
                raise ValueError('unknown marker identifier "{:s}"'.format(locusid))
            locusdata = result.iloc[0]
            context = microhapulator.panel.LocusContext(locusdata)
            coords = list(map(int, locusdata.Offsets.split(',')))
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
    def seqstream(self):
        for locusid in sorted(self.loci()):
            canonid = microhapdb.markers[microhapdb.markers.Name == locusid].iloc[0]
            context = microhapulator.panel.LocusContext(canonid)
            yield context.defline, context.sequence

    @property
    def bedstr(self):
        out = StringIO()
        for line in self.bedstream:
            print(line, file=out)
        return out.getvalue()

    @property
    def haploseqs(self):
        '''Apply genotype to reference and construct full haplotype sequences.'''
        mutator = mutate(self.seqstream, self.bedstream)
        for defline, sequence in mutator:
            yield defline, sequence

    def unmix(self):
        assert self.ploidy % 2 == 0
        ncontrib = int(self.ploidy / 2)
        profiles = [SimulatedProfile() for _ in range(ncontrib)]
        for locus in self.loci():
            for contrib in range(ncontrib):
                for hap in range(2):
                    haplotype = (2 * contrib) + hap
                    allele = self.alleles(locus, haplotype=haplotype).pop()
                    profiles[contrib].add(hap, locus, allele)
        return profiles


class SimulatedProfile(Profile):
    """Profile represented by phased alleles at a number of specified loci.

    Alleles are stored in a dictionary, with microhap locus ID/name as the key
    and a list as the value. Typically each list contains 2 items, the
    haplotype/phase 0 allele and the hap/phase 1 allele. However, if the
    profile represents a mixture there may be more than 2 haplotypes, as
    indicated by the `ploidy` parameter.

    >>> gt = SimulatedProfile()
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
            profile = SimulatedProfile(ploidy=ploidy)
            for locusid, allele_list in locus_alleles.items():
                for i, allele in enumerate(allele_list):
                    profile.add(i, locusid, ','.join(allele))
            return profile

    def merge(profiles):
        ploidy = 2 * len(profiles)
        gt = SimulatedProfile(ploidy=ploidy)
        offset = 0
        for profile in profiles:
            for locusid, locusdata in sorted(profile.data['loci'].items()):
                for allele in locusdata['genotype']:
                    gt.add(offset + allele['haplotype'], locusid, allele['allele'])
            offset += 2
        return gt

    def __init__(self, fromfile=None, ploidy=0):
        super(SimulatedProfile, self).__init__(fromfile=fromfile)
        if ploidy:
            self.data['ploidy'] = ploidy

    def add(self, hapid, locusid, allele):
        if self.data['ploidy'] is not None and self.data['ploidy'] > 0:
            assert hapid in range(self.data['ploidy'])
        if locusid not in self.data['loci']:
            self.data['loci'][locusid] = {'genotype': list()}
        self.data['loci'][locusid]['genotype'].append({'allele': allele, 'haplotype': hapid})

    @property
    def gttype(self):
        return 'SimulatedProfile'


class ObservedProfile(Profile):
    def __init__(self, fromfile=None):
        super(ObservedProfile, self).__init__(fromfile=fromfile)

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

    def infer(self, threshold=10):
        for locusid, locusdata in self.data['loci'].items():
            allelecounts = locusdata['allele_counts']
            eff_cov = 1.0 - (locusdata['num_discarded_reads'] / locusdata['max_coverage'])
            avgcount = 0.0
            if len(allelecounts.values()) > 0:
                avgcount = sum(allelecounts.values()) / len(allelecounts.values())
            gt = set()
            for allele, count in allelecounts.items():
                if eff_cov < 0.25:
                    # Low effective coverage --> use static cutoff
                    if count < threshold:
                        continue
                else:
                    # High effective coverage --> compute cutoff dynamically
                    if count * 4 < avgcount:
                        continue
                gt.add(allele)
            self.data['loci'][locusid]['genotype'] = [
                {'allele': a, 'haplotype': None} for a in sorted(gt)
            ]

    @property
    def gttype(self):
        return 'ObservedProfile'
