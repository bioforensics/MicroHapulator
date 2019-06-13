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
from microhapulator.genotype import SimulatedGenotype
from microhapulator.panel import panel_loci, exclude_loci_missing_data
from microhapulator.panel import validate_populations, sample_panel
import numpy.random
from pyfaidx import Fasta as Fastaidx


def sim(popids, panel, seed=None, relaxed=False):
    '''Simulate a diploid genotype for the specified panel.

    For each locus in the panel simulate a diploid microhap genotype. For each
    haplotype, alleles are selected according to population allele frequencies
    corresponding to the populations indicated in `popids`. If there is no
    population allele frequency data for a particular locus it is excluded
    (except in `relaxed` mode, in which case an allele is sampled from a
    uniform distribution).
    '''
    haplopops = validate_populations(popids)
    loci = panel_loci(panel)
    if not relaxed:
        loci = exclude_loci_missing_data(loci, haplopops)
    if loci in (None, list()):
        raise ValueError('invalid panel: {}'.format(panel))
    genotype = SimulatedGenotype(ploidy=2)
    if seed:
        numpy.random.seed(seed)
    for haplotype, locus, allele in sample_panel(haplopops, loci):
        genotype.add(haplotype, locus, allele)
    message = 'simulated microhaplotype variation at {loc:d} loci'.format(loc=len(loci))
    microhapulator.plog('[MicroHapulator::sim]', message)
    return genotype


def main(args):
    genotype = sim(args.popid, args.panel, seed=args.seed, relaxed=args.relaxed)
    with microhapulator.open(args.out, 'w') as fh:
        genotype.dump(fh)
        message = 'genotype JSON written to {:s}'.format(fh.name)
        microhapulator.plog('[MicroHapulator::sim]', message)
    if args.haplo_seq:
        with microhapulator.open(args.haplo_seq, 'w') as fh:
            for defline, sequence in genotype.haploseqs:
                print('>', defline, '\n', sequence, sep='', file=fh)
            message = 'haplotype sequences written to {:s}'.format(fh.name)
            microhapulator.plog('[MicroHapulator::sim]', message)
