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
from microhapulator import data_file
import numpy
import pytest
import shutil
import tempfile


def test_validate_populations():
    from microhapulator.population import validate_populations as valpop
    assert valpop(['MHDBP000012']) == ['MHDBP000012', 'MHDBP000012']
    assert valpop(['SA000020C', 'SA000012D']) == ['MHDBP000023', 'MHDBP000063']
    assert valpop(['SA000025H', 'MHDBP000063']) == ['MHDBP000063', 'MHDBP000075']


def test_validate_populations_bad_ids():
    from microhapulator.population import validate_populations as valpop
    with pytest.raises(ValueError) as ve:
        _ = valpop(['BogU$pOPiD'])
    assert 'invalid or duplicated population ID(s)' in str(ve)
    with pytest.raises(ValueError) as ve:
        _ = valpop(['NotARealID', 'MHDBP000077'])
    assert 'invalid or duplicated population ID(s)' in str(ve)
    with pytest.raises(ValueError) as ve:
        _ = valpop(['MHDBP000012', 'MHDBP000012'])
    assert 'invalid or duplicated population ID(s)' in str(ve)
    with pytest.raises(ValueError) as ve:
        _ = valpop(['SA000020C', 'MHDBP000023'])
    assert 'invalid or duplicated population ID(s)' in str(ve)


def test_validate_populations_cardinality():
    from microhapulator.population import validate_populations as valpop
    with pytest.raises(ValueError) as ve:
        _ = valpop([])
    assert 'please provide only 1 or 2 population IDs' in str(ve)
    with pytest.raises(ValueError) as ve:
        _ = valpop(['MHDBP000012', 'MHDBP000023', 'MHDBP000063'])
    assert 'please provide only 1 or 2 population IDs' in str(ve)


def test_validate_loci():
    from microhapulator.locus import validate_loci as valloc
    assert valloc(['BogusId']) == []
    assert valloc(['MHDBL000114']) == ['MHDBL000114']
    assert valloc(['MHDBL000079', 'MHDBL000146', 'MHDBL000192']) == ['MHDBL000079', 'MHDBL000146', 'MHDBL000192']  # noqa
    assert valloc(['mh09KK-034', 'mh07PK-38311', 'SI664723C', 'MHDBL000049']) == ['MHDBL000049', 'MHDBL000130', 'MHDBL000198', 'MHDBL000208']  # noqa
    assert valloc(['mh09KK-034', 'mh07PK-38311', 'NotARealId', 'SI664723C', 'MHDBL000049']) == ['MHDBL000049', 'MHDBL000130', 'MHDBL000198', 'MHDBL000208']  # noqa


def test_check_loci_for_population():
    from microhapulator.population import check_loci_for_population as check
    assert check(list(), 'MHDBP000003') == list()
    assert check(['BogusLocus'], 'MHDBP000003') == list()
    assert check(['MHDBL000197'], 'MHDBP000003') == list()
    assert check(['MHDBL000197'], 'MHDBP000004') == ['MHDBL000197']


def test_exclude_loci_missing_data():
    from microhapulator.population import exclude_loci_missing_data as exclude
    assert exclude(['MHDBL000197'], ['MHDBP000003', 'MHDBP000003']) == list()
    assert exclude(['MHDBL000197'], ['MHDBP000003', 'MHDBP000004']) == list()
    assert exclude(['MHDBL000197'], ['MHDBP000004', 'MHDBP000004']) == ['MHDBL000197']
    assert exclude(['MHDBL000066'], ['MHDBP000003', 'MHDBP000004']) == list()


def test_sample_panel():
    pops = ['MHDBP000004', 'MHDBP000004']
    panel = ['MHDBL000197', 'MHDBL000066']
    numpy.random.seed(112358)
    sampler = microhapulator.locus.sample_panel(pops, panel)
    assert list(sampler) == [
        (0, 'MHDBL000197', 'A,A,T,A,T'),
        (0, 'MHDBL000066', 'G,G,G'),
        (1, 'MHDBL000197', 'T,T,A,T,C'),
        (1, 'MHDBL000066', 'G,G,G')
    ]


def test_main():
    tempdir = tempfile.mkdtemp()
    try:
        cli = microhapulator.cli.get_parser()
        arglist = [
            '--panel', 'MHDBL000197', 'MHDBL000066', '--out', tempdir + '/reads.fastq',
            '--num-reads', '500', '--haploseq', tempdir + '/haplo.fasta',
            '--genotype', tempdir + '/genotype.bed', '--hap-seed', '293847',
            '--seq-seed', '123454321', 'hg38.fasta', 'MHDBP000004'
        ]
        args = cli.parse_args(arglist)
        microhapulator.cli.main(args)

        with open(tempdir + '/haplo.fasta', 'r') as fh:
            print('DEBUG\n', fh.read())

        assert filecmp.cmp(tempdir + '/genotype.bed', data_file('alpha.bed'), shallow=False)
        assert filecmp.cmp(tempdir + '/haplo.fasta', data_file('alpha.fasta'), shallow=False)
        assert filecmp.cmp(tempdir + '/reads.fastq', data_file('alpha.fastq'), shallow=False)
    finally:
        shutil.rmtree(tempdir)
