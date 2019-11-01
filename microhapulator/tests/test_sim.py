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
from microhapulator.profile import SimulatedProfile
from microhapulator.tests import data_file
import pytest
import shutil
import tempfile


def test_meaning_of_life():
    profile = microhapulator.sim.sim(['SA004250L', 'SA004250L'], ['beta'], seed=42)
    testprofile = SimulatedProfile(fromfile=data_file('meaning-of-life.json.gz'))
    assert profile == testprofile


@pytest.mark.parametrize('relaxmode,testfile', [
    (False, 'red-strict-profile.json'),
    (True, 'red-relaxed-profile.json'),
])
def test_sim_relaxed(relaxmode, testfile):
    profile = microhapulator.sim.sim(
        ['SA000101C'], ['mh01KK-117', 'mh09KK-157', 'mh07CP-004'],
        seed=54321, relaxed=relaxmode
    )
    testprofile = SimulatedProfile(fromfile=data_file(testfile))
    assert profile == testprofile


def test_main():
    with tempfile.NamedTemporaryFile(suffix='-profile.json') as outfile:
        arglist = [
            'sim', '--out', outfile.name, '--seed', '1985', 'SA004250L', 'SA004250L',
            'usa'
        ]
        args = microhapulator.cli.get_parser().parse_args(arglist)
        microhapulator.sim.main(args)
        p = SimulatedProfile(fromfile=outfile.name)
        testp = SimulatedProfile(fromfile=data_file('bitusa-profile.json'))
        assert p == testp


def test_main_haplo_seq():
    tempdir = tempfile.mkdtemp()
    try:
        arglist = [
            'sim', '--seed', '293847', '--out', tempdir + '/profile.json',
            '--haplo-seq', tempdir + '/haplo.fasta', 'SA004047P', 'SA004047P',
            'mh07CP-004', 'mh14CP-003'
        ]
        args = microhapulator.cli.get_parser().parse_args(arglist)
        microhapulator.sim.main(args)
        p = SimulatedProfile(fromfile=tempdir + '/profile.json')
        testp = SimulatedProfile(fromfile=data_file('orange-sim-profile.json'))
        assert p == testp
        assert filecmp.cmp(tempdir + '/haplo.fasta', data_file('orange-haplo.fasta'))
    finally:
        shutil.rmtree(tempdir)


def test_no_seed():
    genotype = microhapulator.sim.sim(['SA004047P'], ['mh07CP-004', 'mh14CP-003'])
    assert len(genotype.data['markers']) == 2
    assert sorted(genotype.data['markers']) == ['mh07CP-004', 'mh14CP-003']


def test_bad_panel():
    with pytest.raises(ValueError, match=r'invalid panel') as ve:
        microhapulator.sim.sim(['SA004047P'], ['DUUUUDE', 'SWEEEET'])
