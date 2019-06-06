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
from tempfile import NamedTemporaryFile


def test_sim_genotype_roundtrip():
    seed = numpy.random.randint(1, 2**32 - 1)
    print('DEBUG seed:', seed)
    numpy.random.seed(seed)
    panel = ['MHDBL000011', 'MHDBL000091', 'MHDBL000171', 'MHDBL000159', 'MHDBL000069']
    populations = ['MHDBP000004', 'MHDBP000022']
    simulator = microhapulator.panel.sample_panel(populations, panel)

    genotype = SimulatedGenotype(ploidy=2)
    for haplotype, locusid, allele in simulator:
        genotype.add(haplotype, locusid, allele)
    with NamedTemporaryFile() as outfile:
        genotype.dump(outfile.name)
        testgenotype = SimulatedGenotype(fromfile=outfile.name)
        assert testgenotype == genotype

    output = StringIO()
    genotype.dump(output)
    assert output.getvalue() == str(genotype)


def test_alleles():
    simgt = SimulatedGenotype.populate_from_bed(data_file('gttest.bed.gz'))
    obsgt = ObservedGenotype(fromfile=data_file('gttest.json'))
    assert simgt.alleles('BoGuSlOcUs') is None
    assert obsgt.alleles('BoGuSlOcUs') is None
    assert simgt.alleles('MHDBL000135') == set(['G,C,T', 'G,T,C'])
    assert obsgt.alleles('MHDBL000135') == set(['G,C,T', 'G,T,C'])
    assert simgt.alleles('MHDBL000135', haplotype=0) == set(['G,C,T'])
    assert simgt.alleles('MHDBL000135', haplotype=1) == set(['G,T,C'])
    assert obsgt.alleles('MHDBL000135', haplotype=0) == set()


def test_haplotypes():
    simgt = SimulatedGenotype.populate_from_bed(data_file('gttest-mismatch1.bed.gz'))
    assert simgt.haplotypes() == set([0, 1])
    obsgt = ObservedGenotype(data_file('pashtun-sim/test-output.json'))
    assert obsgt.haplotypes() == set()


def test_sim_obs_genotype_equality():
    simgt = SimulatedGenotype.populate_from_bed(data_file('gttest.bed.gz'))
    obsgt = ObservedGenotype(fromfile=data_file('gttest.json'))
    assert simgt == obsgt
    assert obsgt == simgt


def test_sim_obs_genotype_not_equal():
    simgt1 = SimulatedGenotype.populate_from_bed(data_file('gttest-mismatch1.bed.gz'))
    assert simgt1 is not None
    assert simgt1 != 42
    assert simgt1 != 3.14159
    assert simgt1 != 'A,C,C,T'

    obsgt1 = ObservedGenotype(fromfile=data_file('gttest.json'))
    assert simgt1 != obsgt1
    assert obsgt1 != simgt1
    assert obsgt1 != 1985
    assert obsgt1 != 98.6

    simgt2 = SimulatedGenotype.populate_from_bed(data_file('gttest-mismatch2.bed.gz'))
    assert simgt1 != simgt2
    assert simgt2 != obsgt1
    assert obsgt1 != simgt2

    obsgt2 = ObservedGenotype(fromfile=data_file('gttest-altered.json'))
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

    genotype = microhapulator.genotype.SimulatedGenotype.merge([gt1, gt2, gt3])
    print(genotype.bedstr)
    assert genotype.bedstr == (
        'MHDBL000033\t162\t163\tC|C|C|C|C|T\n'
        'MHDBL000033\t163\t164\tG|G|T|T|G|G\n'
        'MHDBL000033\t187\t188\tG|G|A|G|G|G\n'
        'MHDBL000182\t121\t122\tA|A|A|A|G|G\n'
        'MHDBL000182\t228\t229\tC|T|T|T|C|T\n'
    )
