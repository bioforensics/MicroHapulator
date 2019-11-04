#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from happer.mutate import mutate
import microhapulator
from microhapulator.profile import SimulatedProfile
from microhapulator.panel import panel_markers, exclude_markers_missing_freq_data
from microhapulator.panel import validate_populations, sample_panel
import os
import numpy.random


def resolve_panel(panellist):
    panel = list()
    for value in panellist:
        if os.path.isfile(value):
            with open(value, 'r') as fh:
                for line in fh:
                    panel.append(line.strip())
        else:
            panel.append(value)
    return panel


def sim(popids, panel, seed=None, relaxed=False):
    '''Simulate a diploid genotype profile for the specified panel.

    For each locus in the panel simulate a diploid microhap profile. For each
    haplotype, alleles are selected according to population allele frequencies
    corresponding to the populations indicated in `popids`. If there is no
    population allele frequency data for a particular locus it is excluded
    (except in `relaxed` mode, in which case an allele is sampled from a
    uniform distribution).
    '''
    haplopops = validate_populations(popids)
    markers = panel_markers(panel)
    if not relaxed:
        markers = exclude_markers_missing_freq_data(markers, haplopops)
    if len(markers) == 0:
        raise ValueError('invalid panel: {}'.format(panel))
    profile = SimulatedProfile(ploidy=2)
    if seed is None:
        seed = numpy.random.randint(2**32-1)
    numpy.random.seed(seed)
    for haplotype, locus, allele in sample_panel(haplopops, markers):
        profile.add(haplotype, locus, allele)
    profile.data['metadata'] = {
        'MaternalHaploPop': haplopops[0],
        'PaternalHaploPop': haplopops[1],
        'HaploSeed': seed,
    }
    message = 'simulated microhaplotype variation at {loc:d} markers'.format(loc=len(markers))
    microhapulator.plog('[MicroHapulator::sim]', message)
    return profile


def main(args):
    panel = resolve_panel(args.panel)
    profile = sim(args.popid, panel, seed=args.seed, relaxed=args.relaxed)
    with microhapulator.open(args.out, 'w') as fh:
        profile.dump(fh)
        message = 'profile JSON written to {:s}'.format(fh.name)
        microhapulator.plog('[MicroHapulator::sim]', message)
    if args.haplo_seq:
        with microhapulator.open(args.haplo_seq, 'w') as fh:
            for defline, sequence in profile.haploseqs:
                print('>', defline, '\n', sequence, sep='', file=fh)
            message = 'haplotype sequences written to {:s}'.format(fh.name)
            microhapulator.plog('[MicroHapulator::sim]', message)
