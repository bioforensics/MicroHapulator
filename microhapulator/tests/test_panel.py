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
from microhapulator.panel import LocusContext, panel_alpha
import numpy
import pytest


def test_panel_alpha():
    defaultpanel = panel_alpha()
    assert len(defaultpanel) == 22
    assert 'MHDBL000017' in defaultpanel
    assert 'MHDBL000018' not in defaultpanel


@pytest.mark.parametrize('panel,size,included,excluded', [
    (None, 22, 'MHDBL000017', 'MHDBL000018'),
    (['alpha'], 22, 'MHDBL000017', 'MHDBL000018'),
    (['beta'], 50, 'MHDBL000013', 'MHDBL000014'),
    (['mh01KK-117', 'mh01KK-002'], 2, 'MHDBL000014', 'MHDBL000015'),
])
def test_panel_loci(panel, size, included, excluded):
    loci = microhapulator.panel.panel_loci(panel)
    assert len(loci) == size
    assert included in loci
    assert excluded not in loci


def test_validate_populations():
    from microhapulator.panel import validate_populations as valpop
    assert valpop(['MHDBP000012']) == ['MHDBP000012', 'MHDBP000012']
    assert valpop(['SA000020C', 'SA000012D']) == ['MHDBP000023', 'MHDBP000063']
    assert valpop(['SA000025H', 'MHDBP000063']) == ['MHDBP000063', 'MHDBP000075']


def test_validate_populations_bad_ids():
    from microhapulator.panel import validate_populations as valpop
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
    from microhapulator.panel import validate_populations as valpop
    with pytest.raises(ValueError) as ve:
        _ = valpop([])
    assert 'please provide only 1 or 2 population IDs' in str(ve)
    with pytest.raises(ValueError) as ve:
        _ = valpop(['MHDBP000012', 'MHDBP000023', 'MHDBP000063'])
    assert 'please provide only 1 or 2 population IDs' in str(ve)


def test_validate_loci():
    from microhapulator.panel import validate_loci as valloc
    assert valloc(['BogusId']) == []
    assert valloc(['MHDBL000114']) == ['MHDBL000114']
    assert valloc(['MHDBL000079', 'MHDBL000146', 'MHDBL000192']) == ['MHDBL000079', 'MHDBL000146', 'MHDBL000192']  # noqa
    assert valloc(['mh09KK-034', 'mh07PK-38311', 'SI664723C', 'MHDBL000049']) == ['MHDBL000049', 'MHDBL000130', 'MHDBL000198', 'MHDBL000208']  # noqa
    assert valloc(['mh09KK-034', 'mh07PK-38311', 'NotARealId', 'SI664723C', 'MHDBL000049']) == ['MHDBL000049', 'MHDBL000130', 'MHDBL000198', 'MHDBL000208']  # noqa


def test_check_loci_for_population():
    from microhapulator.panel import check_loci_for_population as check
    assert check(list(), 'MHDBP000003') == list()
    assert check(['BogusLocus'], 'MHDBP000003') == list()
    assert check(['MHDBL000197'], 'MHDBP000003') == list()
    assert check(['MHDBL000197'], 'MHDBP000004') == ['MHDBL000197']


def test_exclude_loci_missing_data():
    from microhapulator.panel import exclude_loci_missing_data as exclude
    assert exclude(['MHDBL000197'], ['MHDBP000003', 'MHDBP000003']) == list()
    assert exclude(['MHDBL000197'], ['MHDBP000003', 'MHDBP000004']) == list()
    assert exclude(['MHDBL000197'], ['MHDBP000004', 'MHDBP000004']) == ['MHDBL000197']
    assert exclude(['MHDBL000066'], ['MHDBP000003', 'MHDBP000004']) == list()


def test_sample_panel():
    pops = ['MHDBP000004', 'MHDBP000004']
    panel = ['MHDBL000197', 'MHDBL000066']
    numpy.random.seed(112358)
    sampler = microhapulator.panel.sample_panel(pops, panel)
    assert list(sampler) == [
        (0, 'MHDBL000197', 'A,A,T,A,T'),
        (0, 'MHDBL000066', 'G,G,G'),
        (1, 'MHDBL000197', 'T,T,A,T,C'),
        (1, 'MHDBL000066', 'G,G,G')
    ]


def test_sample_panel_relaxed(capsys):
    pops = ['MHDBP000003', 'MHDBP000004']
    panel = ['MHDBL000172', 'MHDBL000105']
    numpy.random.seed(1776)
    microhapulator.teelog = True
    sampler = microhapulator.panel.sample_panel(pops, panel)
    assert list(sampler) == [
        (0, 'MHDBL000172', 'C,A,G,A'),
        (0, 'MHDBL000105', 'A,T,A,C'),
        (1, 'MHDBL000172', 'C,A,A,G'),
        (1, 'MHDBL000105', 'A,C,A,T')
    ]
    microhapulator.teelog = False
    out, err = capsys.readouterr()
    assert 'no allele frequencies available' in err
    assert 'for population "MHDBP000003" at locus "MHDBL000105"' in err
    assert 'in "relaxed" mode, drawing an allele uniformly' in err


def test_context():
    i, locus = next(microhapdb.loci[microhapdb.loci.ID == 'MHDBL000172'].iterrows())
    c = LocusContext(locus)
    assert c.chrom == 'chr5'
    assert c.global_to_local(1000) is None
    assert c.local_to_global(1000) is None

    c = LocusContext(locus, minlen=100)
    assert len(c) == 197
    c = LocusContext(locus, minlen=100, mindelta=10)
    assert len(c) == 157
