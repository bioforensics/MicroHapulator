#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import filecmp
import microhapdb
import microhapulator
from microhapulator.profile import Profile
from microhapulator.seq import calc_n_reads_from_proportions
from microhapulator.tests import data_file
import numpy.random
import pytest
from tempfile import NamedTemporaryFile


@pytest.mark.parametrize('n,totalreads,prop,result', [
    (3, 100, None, [33, 33, 33]),
    (3, 100, [0.5, 0.4, 0.1], [50, 40, 10]),
    (4, 4000000, [20, 30, 40, 50], [571428, 857142, 1142857, 1428571]),
])
def test_proportions(n, totalreads, prop, result):
    assert calc_n_reads_from_proportions(n, totalreads, prop) == result


def test_proportions_failure_modes():
    message = r'mismatch between contributor number and proportions'
    with pytest.raises(ValueError, match=message) as ve:
        calc_n_reads_from_proportions(3, 1000, [0.6, 0.4])


def test_even_mixture():
    seed = numpy.random.randint(1, 2**32 - 1)
    print('Seed:', seed)
    numpy.random.seed(seed)
    popids = microhapdb.populations[microhapdb.populations.Source == 'ALFRED'].ID.unique()
    profiles = list()
    for _ in range(numpy.random.randint(2, 6)):
        pops = [numpy.random.choice(popids), numpy.random.choice(popids)]
        panel = microhapulator.panel.panel_allpops()[:5]
        p = microhapulator.sim.sim(pops, panel)
        profiles.append(p)
    sequencer = microhapulator.seq.seq(profiles, threads=2, totalreads=500)
    for n, read in enumerate(sequencer):
        pass
    assert n == pytest.approx(500, abs=25)


def test_complex_genotype(capsys):
    profile = Profile(fromfile=data_file('mixture-genotype.json'))
    sequencer = microhapulator.seq.seq(list(profile.unmix()), threads=2, totalreads=200)
    for n, read in enumerate(sequencer):
        pass
    terminal = capsys.readouterr()
    assert terminal.err.count('Individual seed=') == 3


def test_uneven_mixture(capsys):
    panel = ['mh01KK-001', 'mh01KK-205', 'mh01KK-117', 'mh10KK-163']
    pops = ['SA004248S', 'SA004239S', 'SA001530J']
    profiles = [microhapulator.sim.sim([popid], panel) for popid in pops]
    sequencer = microhapulator.seq.seq(
        profiles, threads=2, totalreads=500, proportions=[0.5, 0.3, 0.2]
    )
    for read in sequencer:
        pass
    terminal = capsys.readouterr()
    assert 'numreads=250' in terminal.err
    assert 'numreads=150' in terminal.err
    assert 'numreads=100' in terminal.err


def test_mixture_failure_modes():
    panel = ['mh01KK-001', 'mh01KK-205', 'mh01KK-117', 'mh10KK-163']
    pops = ['SA004248S', 'SA004239S', 'SA001530J']
    profiles = [microhapulator.sim.sim([popid], panel) for popid in pops]

    message = r'number of profiles must match number of seeds'
    with pytest.raises(ValueError, match=message) as ve:
        for read in microhapulator.seq.seq(profiles, seeds=[42, 1776]):
            pass

    message = r'mismatch between contributor number and proportions'
    with pytest.raises(ValueError, match=message) as ve:
        for read in microhapulator.seq.seq(profiles, proportions=[0.5, 0.3, 0.1, 0.1]):
            pass

    message = r'specified proportions result in 0 reads for 1 or more individuals'
    with pytest.raises(ValueError, match=message) as ve:
        for read in microhapulator.seq.seq(profiles, totalreads=500, proportions=[1, 100, 10000]):
            pass


def test_main():
    with NamedTemporaryFile(suffix='.fastq') as outfile:
        arglist = [
            'seq', '--out', outfile.name, '--seeds', '123454321', '--num-reads', '500',
            '--signature', 'srd6Sei', data_file('orange-sim-profile.json')
        ]
        args = microhapulator.cli.get_parser().parse_args(arglist)
        microhapulator.seq.main(args)
        assert filecmp.cmp(outfile.name, data_file('orange-reads.fastq'))


@pytest.mark.parametrize('relaxmode,gtfile,signature,testfile', [
    (False, 'red-strict-profile.json', 'ihkSW9I', 'red-reads-strict.fastq'),
    (True, 'red-relaxed-profile.json', 'kSW9IlM', 'red-reads-relaxed.fastq'),
])
def test_main_relaxed(relaxmode, gtfile, signature, testfile):
    with NamedTemporaryFile(suffix='.fastq') as outfile:
        arglist = [
            'seq', '--out', outfile.name, '--seeds', '24680', '--num-reads', '100',
            '--signature', signature, data_file(gtfile)
        ]
        args = microhapulator.cli.get_parser().parse_args(arglist)
        args.relaxed = relaxmode
        microhapulator.seq.main(args)
        assert filecmp.cmp(outfile.name, data_file(testfile))


def test_main_no_seed():
    with NamedTemporaryFile(suffix='.fastq') as outfile:
        arglist = [
            'seq', '--out', outfile.name, '--num-reads', '200', '--threads', '1',
            data_file('orange-sim-profile.json')
        ]
        args = microhapulator.cli.get_parser().parse_args(arglist)
        microhapulator.seq.main(args)
        with open(outfile.name, 'r') as fh:
            filelines = fh.read().strip().split('\n')
            assert len(filelines) == 800  # 200 reads * 4 lines per read = 800 lines


def test_main_mixture(capsys):
    arglist = [
        'seq', '--seeds', '42', '1776', '--threads', '2', '--proportions', '0.8', '0.2',
        '--num-reads', '500', data_file('yellow-mix-gt.json')
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.seq.main(args)
    terminal = capsys.readouterr()
    outlines = terminal.out.strip().split('\n')
    nrecords = len(outlines) / 4
    assert nrecords == pytest.approx(500, abs=25)
    assert outlines[-3] == (
        'AGTATGTTTTAAGACTCTGAAAATTTTTGAACTCACTCCCAGAAAGTTTTACCACCTCTTCTTCTGTGT'
        'GGCCACCAGGGGGACGTAGTGTGGCCGAGACTCCAGGAGTGCCCGTGAGCACCCGAGGCGCTGAGGAGG'
        'GCTGGGTTGCAGTCTCCTGTGGTTGTACCAGCATTAAAAATCGCTGTATGTGTGTGTGTGTGTGTGTGT'
        'GTGTGCTGAGCCTAAATTTTCTTTGAGCCGCCAATACCTATTATCATGAATCCCTGCCTTGACGCTGAG'
        'GGTAGAAATTGAATTGGATATATGA'
    )
