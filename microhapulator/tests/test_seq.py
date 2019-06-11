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
from microhapulator.seq import calc_n_reads_from_proportions
from microhapulator.tests import data_file
import numpy.random
import pytest
import shutil
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


def test_even_mixture():
    seed = numpy.random.randint(1, 2**32 - 1)
    print('Seed:', seed)
    numpy.random.seed(seed)
    popids = microhapdb.populations[microhapdb.populations.Source == 'ALFRED'].ID.unique()
    genotypes = list()
    for _ in range(numpy.random.randint(2, 6)):
        pops = [numpy.random.choice(popids), numpy.random.choice(popids)]
        panel = microhapulator.panel.panel_allpops()[:5]
        gt = microhapulator.sim.sim(pops, panel)
        genotypes.append(gt)
    sequencer = microhapulator.seq.seq(genotypes, threads=2, totalreads=500)
    for n, read in enumerate(sequencer):
        pass
    assert n == pytest.approx(500, abs=25)


def test_complex_genotype(capsys):
    genotype = microhapulator.genotype.Genotype(fromfile=data_file('mixture-genotype.json'))
    sequencer = microhapulator.seq.seq(list(genotype.unmix()), threads=2, totalreads=200)
    for n, read in enumerate(sequencer):
        pass
    terminal = capsys.readouterr()
    assert terminal.err.count('Individual seed=') == 3


def test_main():
    with NamedTemporaryFile(suffix='.fastq') as outfile:
        arglist = [
            'seq', '--out', outfile.name, '--seeds', '123454321', '--num-reads', '500',
            '--signature', 'srd6Sei', data_file('orange-sim-gt.json')
        ]
        args = microhapulator.cli.parse_args(arglist)
        microhapulator.seq.main(args)
        assert filecmp.cmp(outfile.name, data_file('orange-reads.fastq'))


@pytest.mark.parametrize('relaxmode,gtfile,signature,testfile', [
    (False, 'red-strict-gt.json', 'ihkSW9I', 'red-reads-strict.fastq'),
    (True, 'red-relaxed-gt.json', 'kSW9IlM', 'red-reads-relaxed.fastq'),
])
def test_main_relaxed(relaxmode, gtfile, signature, testfile):
    with NamedTemporaryFile(suffix='.fastq') as outfile:
        arglist = [
            'seq', '--out', outfile.name, '--seeds', '24680', '--num-reads', '100',
            '--signature', signature, data_file(gtfile)
        ]
        args = microhapulator.cli.parse_args(arglist)
        args.relaxed = relaxmode
        microhapulator.seq.main(args)
        assert filecmp.cmp(outfile.name, data_file(testfile))


def test_main_no_seed():
    with NamedTemporaryFile(suffix='.fastq') as outfile:
        arglist = [
            'seq', '--out', outfile.name, '--num-reads', '200', '--threads', '1',
            data_file('orange-sim-gt.json')
        ]
        args = microhapulator.cli.parse_args(arglist)
        microhapulator.seq.main(args)
        with open(outfile.name, 'r') as fh:
            filelines = fh.read().strip().split('\n')
            assert len(filelines) == 800  # 200 reads * 4 lines per read = 800 lines
