#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from itertools import combinations
from math import ceil
import microhapdb
import microhapulator
from numpy.random import choice
import pandas


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


def panel_loci(panellist):
    if panellist is None or panellist == ['alpha']:
        loci = panel_alpha()
    elif panellist == ['beta']:
        loci = panel_beta()
    elif panellist == ['usa']:
        loci = panel_usa()
    elif panellist == ['allpops']:
        loci = panel_allpops()
    else:
        loci = panellist
    return validate_loci(loci)


def panel_alpha():
    '''Initial minimum-effort panel selection.

    There are many factors to consider when designing a panel. This method
    ignores most of those considerations and simply grabs the microhap locus
    from each autosome with the largest effective number of alleles (Ae)
    averaged across all populations.
    '''
    loci = microhapdb.loci.query('Source == "ALFRED"').\
        sort_values('AvgAe', ascending=False).\
        drop_duplicates('Chrom')
    return list(loci.ID)


def panel_beta():
    '''First attempt at optimizing panel design.

    Slightly less unsophisticated approach to panel selection than panel_alpha,
    but still ignoring some important considerations. This method focuses on a
    few simple filters and simple operations.
    - discard any microhap not present in ALFRED
    - discard any microhap with an average Ae of less than 2.0
    - discard any microhap that spans more than 250 bp
    - for each chromosome, grab the 3 microhaps with the highest combined
      average Ae such that no 2 microhaps occur within 25 Mb; if this criterion
      is too strict, reduce the distance and then the number of desired
      microhaps from 3 to 2 until a compatible set is selected
    - combine microhaps from all chromosomes and select the top 50 by AvgAe
    '''
    loci = microhapdb.loci.copy()
    loci['Length'] = loci['End'] - loci['Start']
    locusids = set()
    for chromid in loci.Chrom.unique():
        chromloci = loci[
            (loci.Source == 'ALFRED') &
            (loci.Chrom == chromid) &
            (loci.AvgAe > 2.0) &
            (loci.Length <= 250)
        ]

        def trycombos(n=3, dist=25e6):
            opt_ae, opt_loci = None, None
            for testlocusids in combinations(chromloci.ID, n):
                testloci = loci[loci.ID.isin(testlocusids)]
                for coord1, coord2 in combinations(testloci.Start, 2):
                    if abs(coord1 - coord2) < dist:
                        break
                else:
                    ae = sum(testloci.AvgAe) / len(testloci.AvgAe)
                    if opt_ae is None or ae > opt_ae:
                        opt_ae = ae
                        opt_loci = testlocusids
            return opt_loci
        params = (
            (3, 25e6), (3, 20e6), (2, 25e6), (2, 20e6),
            (3, 15e6), (2, 15e6),
            (3, 10e6), (2, 10e6), (2, 7.5e6),
        )
        for n, dist in params:
            testloci = trycombos(n=n, dist=dist)
            if testloci is not None:
                break
        assert testloci is not None, chromid
        locusids.update(testloci)
    panel = loci[loci.ID.isin(locusids)].sort_values('AvgAe').head(50)
    return list(panel.ID)


def panel_usa():
    '''Loci with frequency data for populations in a sample USA demographic.

    Panel containing the top 100 loci (ranked by average Ae) for which
    population allele frequency data is available for all 19 ALFRED populations
    in a rough simulation of USA demographics.
    '''
    demo_file = microhapulator.package_file('usa-demographics.tsv')
    usa_demo = pandas.read_csv(demo_file, sep='\t')

    prelim_panel = set()
    pops = set(usa_demo.Population)
    for locusid in microhapdb.loci[microhapdb.loci.Source == "ALFRED"].ID.unique():
        testpops = microhapdb.frequencies[
            microhapdb.frequencies.Locus == locusid
        ].Population.unique()
        if pops - set(testpops) == set():
            prelim_panel.add(locusid)

    final_panel = microhapdb.loci[
        microhapdb.loci.ID.isin(prelim_panel)
    ].sort_values('AvgAe', ascending=False).head(100)
    return list(final_panel.ID)


def panel_allpops():
    '''Loci with frequency data for all ALFRED populations.

    Panel containing only loci for which population allele frequency data is
    available for all 96 ALFRED populations.
    '''
    panel = set()
    for locusid in microhapdb.loci[microhapdb.loci.Source == "ALFRED"].ID.unique():
        pops = microhapdb.frequencies[microhapdb.frequencies.Locus == locusid].Population.unique()
        if len(pops) == 96:
            panel.add(locusid)
    return sorted(panel)


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
                microhapulator.plog('[MicroHapulator::panel] WARNING:', message)
                alleles = list(f[f.Locus == locusid].Allele.unique())
                sampled_allele = choice(alleles)
            else:
                alleles = list(allelefreqs.Allele)
                freqs = list(allelefreqs.Frequency)
                freqs = [x / sum(freqs) for x in freqs]
                sampled_allele = choice(alleles, p=freqs)
            yield haplotype, locusid, sampled_allele


def validate_populations(popids):
    if len(popids) not in (1, 2):
        message = 'please provide only 1 or 2 population IDs'
        raise ValueError(message)
    popids = sorted(set(popids))
    haplopops = microhapdb.standardize_ids(popids)
    if len(haplopops) < len(popids):
        raise ValueError('invalid or duplicated population ID(s)')
    if len(haplopops) == 1:
        haplopops = haplopops * 2
    return haplopops


def check_loci_for_population(loci, popid):
    loci = validate_loci(loci)
    if len(loci) == 0:
        return list()
    freq = microhapdb.frequencies
    allelefreqs = freq[(freq.Population == popid) & (freq.Locus.isin(loci))]
    idsfound = list(allelefreqs.Locus.unique())
    nodata = set(loci) - set(idsfound)
    if len(nodata) > 0:
        message = 'no allele frequency data available'
        message += ' for population "{pop:s}"'.format(pop=popid)
        message += ' at the following microhap loci: {loc:s}'.format(loc=','.join(nodata))
        microhapulator.plog('[MicroHapulator::panel] WARNING:', message)
    return sorted(set(loci) & set(idsfound))


def exclude_loci_missing_data(loci, popids):
    loci_to_keep = None
    for popid in popids:
        temp_loci = check_loci_for_population(loci, popid)
        if loci_to_keep is None:
            loci_to_keep = set(temp_loci)
        else:
            loci_to_keep = set(loci_to_keep) & set(temp_loci)
    return sorted(loci_to_keep)
