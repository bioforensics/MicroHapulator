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
from microhapulator.profile import SimulatedProfile, ObservedProfile
from microhapulator.tests import data_file
import numpy
import pytest
from tempfile import NamedTemporaryFile


def test_profile_genotype_roundtrip():
    seed = numpy.random.randint(1, 2**32 - 1)
    print('DEBUG seed:', seed)
    numpy.random.seed(seed)
    panel = ['mh01KK-072', 'mh17KK-014', 'mh05CP-006', 'mh04CP-007', 'mh14KK-101']
    populations = ['SA004047P', 'SA004250L']
    simulator = microhapulator.panel.sample_panel(populations, panel)

    genotype = SimulatedProfile(ploidy=2)
    for haplotype, locusid, allele in simulator:
        genotype.add(haplotype, locusid, allele)
    with NamedTemporaryFile() as outfile:
        genotype.dump(outfile.name)
        testgenotype = SimulatedProfile(fromfile=outfile.name)
        assert testgenotype == genotype

    output = StringIO()
    genotype.dump(output)
    assert output.getvalue() == str(genotype)


def test_alleles():
    simprof = SimulatedProfile.populate_from_bed(data_file('gttest.bed.gz'))
    obsprof = ObservedProfile(fromfile=data_file('gttest.json'))
    assert simprof.alleles('BoGuSlOcUs') is None
    assert obsprof.alleles('BoGuSlOcUs') is None
    assert simprof.alleles('MHDBL000135') == set(['G,C,T', 'G,T,C'])
    assert obsprof.alleles('MHDBL000135') == set(['G,C,T', 'G,T,C'])
    assert simprof.alleles('MHDBL000135', haplotype=0) == set(['G,C,T'])
    assert simprof.alleles('MHDBL000135', haplotype=1) == set(['G,T,C'])
    assert obsprof.alleles('MHDBL000135', haplotype=0) == set()


def test_haplotypes():
    simgt = SimulatedProfile.populate_from_bed(data_file('gttest-mismatch1.bed.gz'))
    assert simgt.haplotypes() == set([0, 1])
    obsgt = ObservedProfile(data_file('pashtun-sim/test-output.json'))
    assert obsgt.haplotypes() == set()


def test_sim_obs_profile_equality():
    simgt = SimulatedProfile.populate_from_bed(data_file('gttest.bed.gz'))
    obsgt = ObservedProfile(fromfile=data_file('gttest.json'))
    assert simgt == obsgt
    assert obsgt == simgt


def test_sim_obs_profile_not_equal():
    simgt1 = SimulatedProfile.populate_from_bed(data_file('gttest-mismatch1.bed.gz'))
    assert simgt1 is not None
    assert simgt1 != 42
    assert simgt1 != 3.14159
    assert simgt1 != 'A,C,C,T'

    obsgt1 = ObservedProfile(fromfile=data_file('gttest.json'))
    assert simgt1 != obsgt1
    assert obsgt1 != simgt1
    assert obsgt1 != 1985
    assert obsgt1 != 98.6

    simgt2 = SimulatedProfile.populate_from_bed(data_file('gttest-mismatch2.bed.gz'))
    assert simgt1 != simgt2
    assert simgt2 != obsgt1
    assert obsgt1 != simgt2

    obsgt2 = ObservedProfile(fromfile=data_file('gttest-altered.json'))
    assert obsgt1 != obsgt2


def test_merge_sim_genotypes():
    prof1 = SimulatedProfile()
    prof1.add(0, 'mh11CP-004', 'C,G,G')
    prof1.add(1, 'mh11CP-004', 'C,G,G')
    prof1.add(0, 'mh05KK-123', 'A,C')
    prof1.add(1, 'mh05KK-123', 'A,T')
    prof2 = SimulatedProfile()
    prof2.add(0, 'mh11CP-004', 'C,T,A')
    prof2.add(1, 'mh11CP-004', 'C,T,G')
    prof2.add(0, 'mh05KK-123', 'A,T')
    prof2.add(1, 'mh05KK-123', 'A,T')
    prof3 = SimulatedProfile()
    prof3.add(0, 'mh11CP-004', 'C,G,G')
    prof3.add(1, 'mh11CP-004', 'T,G,G')
    prof3.add(0, 'mh05KK-123', 'G,C')
    prof3.add(1, 'mh05KK-123', 'G,T')

    profile = microhapulator.profile.SimulatedProfile.merge([prof1, prof2, prof3])
    print(profile.bedstr)
    assert profile.bedstr == (
        'mh05KK-123\t121\t122\tA|A|A|A|G|G\n'
        'mh05KK-123\t228\t229\tC|T|T|T|C|T\n'
        'mh11CP-004\t162\t163\tC|C|C|C|C|T\n'
        'mh11CP-004\t163\t164\tG|G|T|T|G|G\n'
        'mh11CP-004\t187\t188\tG|G|A|G|G|G\n'
    )