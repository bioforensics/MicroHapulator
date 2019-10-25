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

    >>> row = microhapdb.markers[microhapdb.markers.Name == 'mh01KK-172'].iloc[0]
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
        start = int(rowdata.Offsets.split(',')[0]) - mindelta
        end = int(rowdata.Offsets.split(',')[-1]) + mindelta + 1
        initlen = end - start
        if initlen < minlen:
            lendiff = minlen - initlen
            bpextend = ceil(lendiff / 2)
            start -= bpextend
            end += bpextend
        self.start = start
        self.end = end
        ampdata = microhapdb.sequences[microhapdb.sequences.Marker == rowdata.Name].iloc[0]
        ampstart = self.start - ampdata.LeftFlank
        ampend = self.end - ampdata.LeftFlank
        self._sequence = ampdata.Sequence[ampstart:ampend]

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

    @property
    def sequence(self):
        return self._sequence

    @property
    def defline(self):
        gcoords = list(map(int, self._data.Offsets.split(',')))
        coords = [self.global_to_local(x) for x in gcoords]
        return '{loc} GRCh38:{chrom}:{s}-{e} variants={var}'.format(
            loc=self._data.Name, chrom=self.chrom, s=self.start,
            e=self.end, var=':'.join(map(str, coords))
        )


def panel_markers(panellist):
    if panellist is None or panellist == ['alpha']:
        markers = panel_alpha()
    elif panellist == ['beta']:
        markers = panel_beta()
    elif panellist == ['usa']:
        markers = panel_usa()
    elif panellist == ['allpops']:
        markers = panel_allpops()
    else:
        markers = panellist
    return validate_markers(markers)


def panel_alpha():
    '''Initial minimum-effort panel selection.

    There are many factors to consider when designing a panel. This method ignores most of those
    considerations and simply grabs the microhap marker from each autosome with the largest
    effective number of alleles (Ae) averaged across all populations. Only populations from ALFRED
    were considered as these were the most comprehensive published resource at the time.
    '''
    markers = microhapdb.markers.query('Source == "ALFRED"').\
        sort_values('AvgAe', ascending=False).\
        drop_duplicates('Chrom')
    return list(markers.Name)


def panel_beta():
    '''First attempt at optimizing panel design.

    Slightly less unsophisticated approach to panel selection than alpha,
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
    markers = microhapdb.markers.copy()
    markers['Start'] = markers.Offsets.apply(lambda x: min(map(int, x.split(','))))
    markers['End'] = markers.Offsets.apply(lambda x: max(map(int, x.split(','))) + 1)
    markers['Length'] = markers.apply(lambda x: x.End - x.Start, axis=1)
    markerids = set()
    for chromid in markers.Chrom.unique():
        chromloci = markers[
            (markers.Source == 'ALFRED')
            & (markers.Chrom == chromid)
            & (markers.AvgAe > 2.0)
            & (markers.Length <= 250)
        ]

        def trycombos(n=3, dist=25e6):
            opt_ae, opt_markers = None, None
            for testmarkerids in combinations(chromloci.Name, n):
                testmarkers = markers[markers.Name.isin(testmarkerids)]
                for coord1, coord2 in combinations(testmarkers.Start, 2):
                    if abs(coord1 - coord2) < dist:
                        break
                else:
                    ae = sum(testmarkers.AvgAe) / len(testmarkers.AvgAe)
                    if opt_ae is None or ae > opt_ae:
                        opt_ae = ae
                        opt_markers = testmarkerids
            return opt_markers
        params = (
            (3, 25e6), (3, 20e6), (2, 25e6), (2, 20e6),
            (3, 15e6), (2, 15e6),
            (3, 10e6), (2, 10e6), (2, 7.5e6),
        )
        for n, dist in params:
            testmarkers = trycombos(n=n, dist=dist)
            if testmarkers is not None:
                markerids.update(testmarkers)
                break
    panel = markers[markers.Name.isin(markerids)].sort_values('AvgAe').head(50)
    return list(panel.Name)


def panel_usa():
    '''Loci with frequency data for populations in a sample USA demographic.

    Panel containing the top 50 loci (ranked by average Ae) composed of 3 or
    more variants and for which population allele frequency data is available
    for all 19 ALFRED sub-populations in a mock population roughly reflecting
    demographics of the United States.
    '''
    demo_file = microhapulator.package_file('usa-demographics.tsv')
    usa_demo = pandas.read_csv(demo_file, sep='\t')
    pops = set(usa_demo.Population)
    markerids = exclude_markers_missing_freq_data(
        microhapdb.markers.Name.unique(), pops, suppress=True
    )
    markers = microhapdb.markers.copy()
    markers['NumVariants'] = markers.Offsets.apply(lambda x: x.count(',') + 1)
    panel = markers[
        (markers.Name.isin(markerids)) & (markers.NumVariants > 2)
    ].sort_values('AvgAe', ascending=False).head(50)
    return list(panel.Name)


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


def validate_markers(panel):
    valid_markers = microhapdb.marker.standardize_ids(panel)
    if len(valid_markers) < len(panel):
        message = 'panel includes duplicate and/or invalid marker IDs'
        microhapulator.plog('[MicroHapulator::panel] WARNING', message)
    return valid_markers


def sample_panel(popids, markers):
    for haplotype, popid in enumerate(popids):
        for markerid in markers:
            f = microhapdb.frequencies
            allelefreqs = f[(f.Population == popid) & (f.Marker == markerid)]
            if len(allelefreqs) == 0:
                message = 'no allele frequencies available'
                message += ' for population "{pop}"'.format(pop=popid)
                message += ' at marker "{loc}"'.format(loc=markerid)
                message += '; in "relaxed" mode, drawing an allele uniformly'
                microhapulator.plog('[MicroHapulator::panel] WARNING:', message)
                alleles = list(f[f.Marker == markerid].Allele.unique())
                sampled_allele = choice(alleles)
            else:
                alleles = list(allelefreqs.Allele)
                freqs = list(allelefreqs.Frequency)
                freqs = [x / sum(freqs) for x in freqs]
                sampled_allele = choice(alleles, p=freqs)
            yield haplotype, markerid, sampled_allele


def validate_populations(popids):
    if len(popids) not in (1, 2):
        message = 'please provide only 1 or 2 population IDs'
        raise ValueError(message)
    popids = sorted(set(popids))
    haplopops = list(microhapdb.population.standardize_ids(popids))
    if len(haplopops) < len(popids):
        raise ValueError('invalid or duplicated population ID(s)')
    if len(haplopops) == 1:
        haplopops = haplopops * 2
    return haplopops


def check_markers_for_population(markers, popid, suppress=False):
    markers = validate_markers(markers)
    if len(markers) == 0:
        return list()
    freq = microhapdb.frequencies
    allelefreqs = freq[(freq.Population == popid) & (freq.Marker.isin(markers))]
    idsfound = list(allelefreqs.Marker.unique())
    nodata = set(markers) - set(idsfound)
    if len(nodata) > 0 and not suppress:
        message = 'no allele frequency data available'
        message += ' for population "{pop:s}"'.format(pop=popid)
        message += ' at the following microhap markers: {loc:s}'.format(loc=','.join(nodata))
        microhapulator.plog('[MicroHapulator::panel] WARNING:', message)
    return sorted(set(markers) & set(idsfound))


def exclude_markers_missing_freq_data(markers, popids, suppress=False):
    markers_to_keep = set(markers)
    for popid in popids:
        temp_markers = check_markers_for_population(markers, popid, suppress=suppress)
        markers_to_keep = markers_to_keep & set(temp_markers)
    return sorted(markers_to_keep)
