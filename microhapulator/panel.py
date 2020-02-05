#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from itertools import combinations
import microhapdb
import microhapulator
from numpy.random import choice
import pandas


# def panel_markers(panellist):
#     if panellist is None or panellist == ['alpha']:
#         markers = panel_alpha()
#     elif panellist == ['beta']:
#         markers = panel_beta()
#     elif panellist == ['usa']:
#         markers = panel_usa()
#     elif panellist == ['allpops']:
#         markers = panel_allpops()
#     else:
#         markers = panellist
#     return validate_markers(markers)


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
