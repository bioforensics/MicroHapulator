# -------------------------------------------------------------------------------------------------
# Copyright (c) 2018, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from microhapulator.happer.allele import Allele, InvalidGenomicCoordinateError, parse_alleles
from microhapulator.happer.tests import data_file
import pytest
from random import shuffle


@pytest.fixture
def allele1():
    return Allele("contig1776", 2468, 2469, "A")


@pytest.fixture
def allele2():
    return Allele("contig1776", 2468, 2469, "CAT", refrseq="G")


@pytest.fixture
def allele3():
    return Allele("contig1812", 13579, 13580, "C", refrseq="T")


@pytest.fixture
def allele4():
    return Allele("contig1812", 13579, 13580, "G", refrseq="T")


@pytest.fixture
def allele5():
    return Allele("contig1776", 123456788, 123456789, "T")


@pytest.fixture
def allele6():
    return Allele("chr17", 1944, 1945, "GATTACA")


def test_basic(allele1):
    assert allele1.seqid == "contig1776"
    assert allele1.start == 2468
    assert allele1.end == 2469
    assert allele1.seq == "A"
    assert allele1.refr is None
    assert len(allele1) == 1
    assert allele1.refrlength == 1


def test_bad_coord():
    with pytest.raises(InvalidGenomicCoordinateError) as e:
        Allele("scaffold1492", 3000, 2999, "C")
    assert "allele end cannot be before allele start" in str(e)


def test_indel(allele2):
    assert allele2.refrlength == 1
    assert len(allele2) == 3
    assert allele2.refr == "G"
    assert allele2 < Allele("contig1776", 2468, 2500, "ATG")


def test_transform(allele3, allele4):
    allele3.transform(100)
    assert allele3.start == 13679
    assert allele3.end == 13680

    allele4.transform(-1000)
    assert allele4.start == 12579
    assert allele4.end == 12580

    with pytest.raises(InvalidGenomicCoordinateError) as e:
        allele3.transform(-20000)
    assert "invalid allele transformation" in str(e)


def test_compare(allele1, allele2, allele3, allele4, allele5, allele6):
    assert allele1 == Allele("contig1776", 2468, 2469, "A", refrseq="G")
    assert allele3 < allele4
    for _ in range(5):
        alleles = [allele2, allele4, allele3, allele1, allele6, allele5]
        shuffle(alleles)
        sa = sorted(alleles)
        print([a.slug for a in sa])
        assert sa == [allele6, allele1, allele2, allele5, allele3, allele4]


def test_parse_alleles():
    infile = data_file("ind1.bed")
    with open(infile, "r") as instream:
        genotypes = list(parse_alleles(instream))
    assert len(genotypes) == 4
    assert [len(g) for g in genotypes] == [2, 2, 2, 2]
    assert genotypes[0][0].slug == "chr1:55-56"
    assert genotypes[3][1].slug == "chr1:88-89"
