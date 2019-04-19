#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from io import StringIO
import microhapdb
import microhapulator
from microhapulator.genotype import SimulatedGenotype, ObservedGenotype
from microhapulator.tests import data_file
import numpy
import pytest


def test_sim_genotype_roundtrip():
    seed = numpy.random.randint(1, 2**32 - 1)
    print('DEBUG seed:', seed)
    numpy.random.seed(seed)
    output = StringIO()
    panel = ['MHDBL000011', 'MHDBL000091', 'MHDBL000171', 'MHDBL000159', 'MHDBL000069']
    populations = ['MHDBP000004', 'MHDBP000022']
    genotype = SimulatedGenotype()
    simulator = microhapulator.panel.sample_panel(populations, panel)
    for haplotype, locusid, allele in simulator:
        genotype.add(haplotype, locusid, allele)
    print(genotype, file=output)
    bedlines = output.getvalue().split('\n')
    testgenotype = SimulatedGenotype(frombed=bedlines)
    assert testgenotype == genotype


@pytest.mark.parametrize('simgtfile,obsgtfile,same', [
    ('gttest.bed.gz', 'gttest.json', True),
    ('gttest-mismatch.bed.gz', 'gttest.json', False),
])
def test_sim_obs_genotype_compare(simgtfile, obsgtfile, same):
    with microhapulator.open(data_file(simgtfile), 'r') as fh:
        simgt = SimulatedGenotype(frombed=fh)
    obsgt = ObservedGenotype(filename=data_file(obsgtfile))
    assert (simgt == obsgt) is same
    assert (obsgt == simgt) is same
