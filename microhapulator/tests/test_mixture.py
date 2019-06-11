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


@pytest.mark.known_failing
def test_uneven_mixture(capfd):
    simulator = microhapulator.mixture.mixture(
        [['MHDBP000021'], ['MHDBP000009'], ['MHDBP000081']],
        ['MHDBL000002', 'MHDBL000003', 'MHDBL000007', 'MHDBL000013', 'MHDBL000017'],
        seqthreads=2, totalreads=500, proportions=[0.5, 0.3, 0.2]
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
            indivs, panel, totalreads=500, hapseeds=[42, 1776]
        )
        list(simulator)
    assert 'individuals must match number of "--hap-seeds" and "--seq-seeds"' in str(ve)

    with pytest.raises(ValueError) as ve:
        simulator = microhapulator.mixture.mixture(
            indivs, panel, totalreads=500, proportions=[0.5, 0.3, 0.1, 0.1]
        )
        list(simulator)
    assert 'mismatch between contributor number and proportions' in str(ve)

    with pytest.raises(ValueError) as ve:
        simulator = microhapulator.mixture.mixture(
            indivs, panel, totalreads=500, proportions=[1, 100, 10000]
        )
        list(simulator)
    assert 'specified proportions result in 0 reads for 1 or more individuals' in str(ve)


def test_mixture_main(capsys):
    arglist = [
        'mixture', '--indiv', 'MHDBP000052', '--indiv', 'MHDBP000052',
        '--panel', 'MHDBL000002', 'MHDBL000003', 'MHDBL000007', 'MHDBL000013', 'MHDBL000017',
        '--proportions', '0.8', '0.2', '--num-reads', '500', '--seq-threads', '2',
        '--hap-seeds', '3579', '2468', '--seq-seeds', '42', '1776',
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
