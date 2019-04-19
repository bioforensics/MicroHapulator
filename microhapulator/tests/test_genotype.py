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


def test_sim_obs_genotype_equal():
    with microhapulator.open(data_file('gttest.bed.gz'), 'r') as fh:
        simgt = SimulatedGenotype(frombed=fh)
    obsgt = ObservedGenotype(filename=data_file('gttest.json'))
    assert simgt == obsgt
    assert obsgt == simgt


def test_sim_obs_genotype_not_equal():
    with microhapulator.open(data_file('gttest-mismatch1.bed.gz'), 'r') as fh:
        simgt = SimulatedGenotype(frombed=fh)
    assert simgt is not None
    assert simgt != 42
    assert simgt != 3.14159
    assert simgt != 'A,C,C,T'

    obsgt = ObservedGenotype(filename=data_file('gttest.json'))
    assert simgt != obsgt
    assert obsgt != simgt

    with microhapulator.open(data_file('gttest-mismatch2.bed.gz'), 'r') as fh:
        simgt = SimulatedGenotype(frombed=fh)
