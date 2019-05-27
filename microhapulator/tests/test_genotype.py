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


def test_alleles():
    with microhapulator.open(data_file('gttest.bed.gz'), 'r') as fh:
        simgt = SimulatedGenotype(frombed=fh)
    obsgt = ObservedGenotype(filename=data_file('gttest.json'))
    assert simgt.alleles('BoGuSlOcUs') is None
    assert obsgt.alleles('BoGuSlOcUs') is None


def test_sim_obs_genotype_equal():
    with microhapulator.open(data_file('gttest.bed.gz'), 'r') as fh:
        simgt = SimulatedGenotype(frombed=fh)
    obsgt = ObservedGenotype(filename=data_file('gttest.json'))
    assert simgt == obsgt
    assert obsgt == simgt


def test_sim_obs_genotype_not_equal():
    with microhapulator.open(data_file('gttest-mismatch1.bed.gz'), 'r') as fh:
        simgt1 = SimulatedGenotype(frombed=fh)
    assert simgt1 is not None
    assert simgt1 != 42
    assert simgt1 != 3.14159
    assert simgt1 != 'A,C,C,T'

    obsgt1 = ObservedGenotype(filename=data_file('gttest.json'))
    assert simgt1 != obsgt1
    assert obsgt1 != simgt1
    assert obsgt1 != 1985
    assert obsgt1 != 98.6

    with microhapulator.open(data_file('gttest-mismatch2.bed.gz'), 'r') as fh:
        simgt2 = SimulatedGenotype(frombed=fh)
    assert simgt1 != simgt2
    assert simgt2 != obsgt1
    assert obsgt1 != simgt2

    obsgt2 = ObservedGenotype(filename=data_file('gttest-altered.json'))
    assert obsgt1 != obsgt2


def test_merge_sim_genotypes():
    gt1 = SimulatedGenotype()
    gt1.add(0, 'MHDBL000033', 'C,G,G')
    gt1.add(1, 'MHDBL000033', 'C,G,G')
    gt1.add(0, 'MHDBL000182', 'A,C')
    gt1.add(1, 'MHDBL000182', 'A,T')
    gt2 = SimulatedGenotype()
    gt2.add(0, 'MHDBL000033', 'C,T,A')
    gt2.add(1, 'MHDBL000033', 'C,T,G')
    gt2.add(0, 'MHDBL000182', 'A,T')
    gt2.add(1, 'MHDBL000182', 'A,T')
    gt3 = SimulatedGenotype()
    gt3.add(0, 'MHDBL000033', 'C,G,G')
    gt3.add(1, 'MHDBL000033', 'T,G,G')
    gt3.add(0, 'MHDBL000182', 'G,C')
    gt3.add(1, 'MHDBL000182', 'G,T')

    genotype = microhapulator.genotype.merge_simulated_genotypes([gt1, gt2, gt3])
    assert str(genotype) == (
        'MHDBL000033\t162\t163\tC|C|C|C|C|T\n'
        'MHDBL000033\t163\t164\tG|G|T|T|G|G\n'
        'MHDBL000033\t187\t188\tG|G|A|G|G|G\n'
        'MHDBL000182\t121\t122\tA|A|A|A|G|G\n'
        'MHDBL000182\t228\t229\tC|T|T|T|C|T\n'
    )
