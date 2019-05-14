#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import microhapdb
import microhapulator
from microhapulator.mixture import calc_n_reads_from_proportions
from microhapulator.tests import data_file
import numpy
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
    with pytest.raises(ValueError) as ve:
        calc_n_reads_from_proportions(3, 1000, [0.6, 0.4])
    assert 'mismatch between contributor number and proportions' in str(ve)


def test_even_mixture(capsys):
    seed = numpy.random.randint(1, 2**32 - 1)
    print('Seed:', seed)
    numpy.random.seed(seed)
    n = numpy.random.randint(2, 6)
    popids = microhapdb.populations[microhapdb.populations.Source == 'ALFRED'].ID.unique()
    indiv_pops = [(numpy.random.choice(popids), numpy.random.choice(popids)) for _ in range(n)]
    panel = microhapulator.panel.panel_allpops()[:5]
    simulator = microhapulator.mixture.mixture(
        indiv_pops, panel, 'hg38.fasta', totalreads=500, seqthreads=2,
    )
    for m, read in enumerate(simulator):
        pass
    assert m == pytest.approx(500, abs=25)


@pytest.mark.known_failing
def test_uneven_mixture(capfd):
    simulator = microhapulator.mixture.mixture(
        [['MHDBP000021'], ['MHDBP000009'], ['MHDBP000081']],
        ['MHDBL000002', 'MHDBL000003', 'MHDBL000007', 'MHDBL000013', 'MHDBL000017'],
        'hg38.fasta', seqthreads=2, totalreads=500, proportions=[0.5, 0.3, 0.2]
    )
    for read in simulator:
        pass
    terminal = capfd.readouterr()
    assert 'numreads=250' in terminal.err
    assert 'numreads=150' in terminal.err
    assert 'numreads=100' in terminal.err


def test_mixture_failure_modes():
    indivs = [['MHDBP000021'], ['MHDBP000009'], ['MHDBP000081']]
    panel = ['MHDBL000002', 'MHDBL000003', 'MHDBL000007', 'MHDBL000013', 'MHDBL000017']

    with pytest.raises(ValueError) as ve:
        simulator = microhapulator.mixture.mixture(
            indivs, panel, 'hg38.fasta', totalreads=500, hapseeds=[42, 1776]
        )
        list(simulator)
    assert 'individuals must match number of "--hap-seeds" and "--seq-seeds"' in str(ve)

    with pytest.raises(ValueError) as ve:
        simulator = microhapulator.mixture.mixture(
            indivs, panel, 'hg38.fasta', totalreads=500, proportions=[0.5, 0.3, 0.1, 0.1]
        )
        list(simulator)
    assert 'mismatch between contributor number and proportions' in str(ve)

    with pytest.raises(ValueError) as ve:
        simulator = microhapulator.mixture.mixture(
            indivs, panel, 'hg38.fasta', totalreads=500, proportions=[1, 100, 10000]
        )
        list(simulator)
    assert 'specified proportions result in 0 reads for 1 or more individuals' in str(ve)


def test_mixture_main(capsys):
    arglist = [
        'mixture', '--indiv', 'MHDBP000052', '--indiv', 'MHDBP000052',
        '--panel', 'MHDBL000002', 'MHDBL000003', 'MHDBL000007', 'MHDBL000013', 'MHDBL000017',
        '--proportions', '0.8', '0.2', '--num-reads', '500', '--seq-threads', '2',
        '--hap-seeds', '3579', '2468', '--seq-seeds', '42', '1776', '--', 'hg38.fasta'
    ]
    args = microhapulator.cli.parse_args(arglist)
    microhapulator.mixture.main(args)
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
