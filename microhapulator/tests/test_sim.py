#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import filecmp
import microhapulator
from microhapulator.genotype import SimulatedGenotype
from microhapulator.tests import data_file
import pytest
import shutil
import tempfile


def test_meaning_of_life():
    genotype = microhapulator.sim.sim(['MHDBP000022', 'MHDBP000022'], ['beta'], seed=42)
    testgenotype = SimulatedGenotype(fromfile=data_file('meaning-of-life.json.gz'))
    assert genotype == testgenotype


@pytest.mark.parametrize('relaxmode,testfile', [
    (False, 'red-strict-gt.json'),
    (True, 'red-relaxed-gt.json'),
])
def test_sim_relaxed(relaxmode, testfile):
    genotype = microhapulator.sim.sim(
        ['MHDBP000003'], ['MHDBL000013', 'MHDBL000212', 'MHDBL000197'],
        seed=54321, relaxed=relaxmode
    )
    testgenotype = SimulatedGenotype(fromfile=data_file(testfile))
    assert genotype == testgenotype


def test_main():
    with tempfile.NamedTemporaryFile(suffix='.gt.json') as outfile:
        arglist = [
            'sim', '--out', outfile.name, '--seed', '1985', 'MHDBP000022', 'MHDBP000022',
            'usa'
        ]
        args = microhapulator.cli.get_parser().parse_args(arglist)
        microhapulator.sim.main(args)
        gt = SimulatedGenotype(fromfile=outfile.name)
        testgt = SimulatedGenotype(fromfile=data_file('bitusa-gt.json'))
        import subprocess
        subprocess.check_call(['cp', outfile.name, 'DUDE.json'])
        assert gt == testgt


def test_main_haplo_seq():
    tempdir = tempfile.mkdtemp()
    try:
        arglist = [
            'sim', '--seed', '293847', '--out', tempdir + '/genotype.json',
            '--haplo-seq', tempdir + '/haplo.fasta', 'MHDBP000004', 'MHDBP000004',
            'MHDBL000197', 'MHDBL000066'
        ]
        args = microhapulator.cli.get_parser().parse_args(arglist)
        microhapulator.sim.main(args)
        gt = SimulatedGenotype(fromfile=tempdir + '/genotype.json')
        testgt = SimulatedGenotype(fromfile=data_file('orange-sim-gt.json'))
        assert gt == testgt
        assert filecmp.cmp(tempdir + '/haplo.fasta', data_file('orange-haplo.fasta'))
    finally:
        shutil.rmtree(tempdir)


def test_no_seed():
    genotype = microhapulator.sim.sim(['MHDBP000004'], ['MHDBL000197', 'MHDBL000066'])
    assert len(genotype.data['loci']) == 2
    assert sorted(genotype.data['loci']) == ['MHDBL000066', 'MHDBL000197']


def test_bad_panel():
    with pytest.raises(ValueError, match=r'invalid panel') as ve:
        microhapulator.sim.sim(['MHDBP000004'], ['DUUUUDE', 'SWEEEET'])
