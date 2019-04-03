#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import microhapulator
import numpy
import pytest


def test_validate_populations():
    from microhapulator.cli import validate_populations as valpop
    assert valpop(['MHDBP000012']) == ['MHDBP000012', 'MHDBP000012']
    assert valpop(['MHDBP000012', 'MHDBP000012']) == ['MHDBP000012', 'MHDBP000012']
    assert valpop(['SA000020C', 'SA000012D']) == ['MHDBP000023', 'MHDBP000063']
    assert valpop(['SA000025H', 'MHDBP000063']) == ['MHDBP000075', 'MHDBP000063']


def test_validate_populations_bad_ids():
    from microhapulator.cli import validate_populations as valpop
    with pytest.raises(ValueError) as ve:
        _ = valpop(['BogU$pOPiD'])
    assert 'invalid population ID' in str(ve)
    with pytest.raises(ValueError) as ve:
        _ = valpop(['NotARealID', 'MHDBP000077'])
    assert 'invalid population ID' in str(ve)


def test_validate_populations_cardinality():
    from microhapulator.cli import validate_populations as valpop
    with pytest.raises(ValueError) as ve:
        _ = valpop([])
    assert 'please provide only 1 or 2 population IDs' in str(ve)
    with pytest.raises(ValueError) as ve:
        _ = valpop(['MHDBP000012', 'MHDBP000023', 'MHDBP000063'])
    assert 'please provide only 1 or 2 population IDs' in str(ve)


def test_validate_loci():
    from microhapulator.cli import validate_loci as valloc
    assert valloc(['MHDBP000003', 'MHDBP000003'], panel=['MHDBL000197']) == list()
    assert valloc(['MHDBP000003', 'MHDBP000004'], panel=['MHDBL000197']) == list()
    assert valloc(['MHDBP000004', 'MHDBP000004'], panel=['MHDBL000197']) == ['MHDBL000197']
    assert valloc(['MHDBP000003', 'MHDBP000004'], panel=['MHDBL000066']) == list()


def test_sample_panel():
    pops = ['MHDBP000004', 'MHDBP000004']
    panel = ['MHDBL000197', 'MHDBL000066']
    numpy.random.seed(112358)
    sampler = microhapulator.cli.sample_panel(pops, panel)
    assert list(sampler) == [
        (0, 'MHDBL000197', 'A,A,T,A,T'),
        (0, 'MHDBL000066', 'G,G,G'),
        (1, 'MHDBL000197', 'T,T,A,T,C'),
        (1, 'MHDBL000066', 'G,G,G')
    ]
