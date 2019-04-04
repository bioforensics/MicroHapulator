#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


import microhapdb
from microhapulator.locus import validate_loci
from sys import stderr


def validate_populations(popids):
    if len(popids) not in (1, 2):
        message = 'please provide only 1 or 2 population IDs'
        raise ValueError(message)
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
    f = microhapdb.frequencies
    allelefreqs = f[(f.Population == popid) & (f.Locus.isin(loci))]
    idsfound = list(allelefreqs.Locus.unique())
    nodata = set(loci) - set(idsfound)
    if len(nodata) > 0:
        message = 'no allele frequency data available'
        message += ' for population "{pop:s}"'.format(pop=popid)
        message += ' at the following microhap loci: {loc:s}'.format(loc=','.join(nodata))
        print('[MicroHapulator::loci] WARNING:', message, file=stderr)
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
