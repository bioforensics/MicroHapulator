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
from microhapulator.tests import data_file
import pytest
import shutil
import tempfile


def test_main():
    tempdir = tempfile.mkdtemp()
    try:
        arglist = [
            'sim', '--panel', 'MHDBL000197', 'MHDBL000066', '--out', tempdir + '/reads.fastq',
            '--num-reads', '500', '--haploseq', tempdir + '/haplo.fasta',
            '--genotype', tempdir + '/genotype.bed', '--hap-seed', '293847',
            '--seq-seed', '123454321', 'MHDBP000004'
        ]
        args = microhapulator.cli.parse_args(arglist)
        microhapulator.sim.main(args)

        gt = microhapulator.genotype.SimulatedGenotype(fromfile=tempdir + '/genotype.bed')
        testgt = microhapulator.genotype.SimulatedGenotype(fromfile=data_file('alpha.sim.json'))
        assert gt == testgt
        assert filecmp.cmp(tempdir + '/haplo.fasta', data_file('alpha.fasta'))
        assert filecmp.cmp(tempdir + '/reads.fastq', data_file('alpha.fastq'))
    finally:
        shutil.rmtree(tempdir)


def test_main_no_haploseq():
    tempdir = tempfile.mkdtemp()
    try:
        arglist = [
            'sim', '--panel', 'MHDBL000197', 'MHDBL000066', '--out', tempdir + '/reads.fastq',
            '--num-reads', '250', '--hap-seed', '1234', '--seq-seed', '5678',
            '--seq-threads', '2', 'MHDBP000004',
        ]
        args = microhapulator.cli.parse_args(arglist)
        microhapulator.sim.main(args)
        assert filecmp.cmp(tempdir + '/reads.fastq', data_file('beta.fastq'))
    finally:
        shutil.rmtree(tempdir)


@pytest.mark.parametrize('relaxed,testfile', [
    (False, 'gamma-strict.fastq'),
    (True, 'gamma-relaxed.fastq'),
])
def test_main_relaxed(relaxed, testfile):
    tempdir = tempfile.mkdtemp()
    try:
        arglist = [
            'sim', '--panel', 'MHDBL000013', 'MHDBL000212', 'MHDBL000197', '--num-reads', '100',
            '--hap-seed', '54321', '--seq-seed', '24680', '--out', tempdir + '/reads.fastq',
            'MHDBP000003',
        ]
        args = microhapulator.cli.parse_args(arglist)
        args.relaxed = relaxed
        microhapulator.sim.main(args)
        assert filecmp.cmp(tempdir + '/reads.fastq', data_file(testfile))
    finally:
        shutil.rmtree(tempdir)


def test_main_no_seeds():
    tempdir = tempfile.mkdtemp()
    try:
        arglist = [
            'sim', '--panel', 'MHDBL000197', 'MHDBL000066', '--out', tempdir + '/reads.fastq',
            '--num-reads', '200', '--seq-threads', '1', 'MHDBP000004',
        ]
        args = microhapulator.cli.parse_args(arglist)
        microhapulator.sim.main(args)
        with open(tempdir + '/reads.fastq', 'r') as fh:
            filelines = fh.read().strip().split('\n')
            assert len(filelines) == 800  # 200 reads * 4 lines per read = 800 lines
    finally:
        shutil.rmtree(tempdir)


def test_main_bad_panel():
    with tempfile.NamedTemporaryFile(suffix='fq.gz') as outfile:
        arglist = [
            'sim', '--panel', 'DUUUUDE', 'SWEEEET', '--num-reads', '10',
            '-o', outfile.name, 'MHDBP000004'
        ]
        args = microhapulator.cli.parse_args(arglist)
        with pytest.raises(ValueError) as ve:
            microhapulator.sim.main(args)
        assert 'invalid panel' in str(ve)
