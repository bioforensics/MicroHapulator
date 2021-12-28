# -------------------------------------------------------------------------------------------------
# Copyright (c) 2019, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

import microhapulator
from microhapulator.profile import SimulatedProfile, ObservedProfile
from microhapulator.tests import data_file
import numpy
import pandas as pd
import pytest


def test_profile_roundtrip(tmp_path):
    seed = numpy.random.randint(1, 2 ** 32 - 1)
    freqs = pd.read_csv(data_file("freq/asw5-freq.tsv"), sep="\t")
    profile = microhapulator.op.sim(freqs, seed=seed)
    profile.dump(tmp_path / "profile.json")
    test = SimulatedProfile(fromfile=tmp_path / "profile.json")
    assert profile == test
    assert str(profile) == str(test)


def test_alleles():
    simprof = SimulatedProfile.populate_from_bed(data_file("gttest.bed.gz"))
    obsprof = ObservedProfile(fromfile=data_file("prof/gttest.json"))
    assert simprof.alleles("BoGuSlOcUs") == set()
    assert obsprof.alleles("BoGuSlOcUs") == set()
    assert simprof.alleles("MHDBL000135") == set(["G,C,T", "G,T,C"])
    assert obsprof.alleles("MHDBL000135") == set(["G,C,T", "G,T,C"])
    assert simprof.alleles("MHDBL000135", haplotype=0) == set(["G,C,T"])
    assert simprof.alleles("MHDBL000135", haplotype=1) == set(["G,T,C"])
    assert obsprof.alleles("MHDBL000135", haplotype=0) == set()


def test_haplotypes():
    simprof = SimulatedProfile.populate_from_bed(data_file("gttest-mismatch1.bed.gz"))
    assert simprof.haplotypes() == set([0, 1])
    obsprof = ObservedProfile(data_file("pashtun-sim/test-output.json"))
    assert obsprof.haplotypes() == set()


def test_sim_obs_profile_equality():
    simprof = SimulatedProfile.populate_from_bed(data_file("gttest.bed.gz"))
    obsprof = ObservedProfile(fromfile=data_file("prof/gttest.json"))
    assert simprof == obsprof
    assert obsprof == simprof


def test_sim_obs_profile_not_equal():
    simprof1 = SimulatedProfile.populate_from_bed(data_file("gttest-mismatch1.bed.gz"))
    assert simprof1 is not None
    assert simprof1 != 42
    assert simprof1 != 3.14159
    assert simprof1 != "A,C,C,T"

    obsprof1 = ObservedProfile(fromfile=data_file("prof/gttest.json"))
    assert simprof1 != obsprof1
    assert obsprof1 != simprof1
    assert obsprof1 != 1985
    assert obsprof1 != 98.6

    simprof2 = SimulatedProfile.populate_from_bed(data_file("gttest-mismatch2.bed.gz"))
    assert simprof1 != simprof2
    assert simprof2 != obsprof1
    assert obsprof1 != simprof2

    obsprof2 = ObservedProfile(fromfile=data_file("prof/gttest-altered.json"))
    assert obsprof1 != obsprof2


def test_merge_sim_genotypes():
    prof1 = SimulatedProfile()
    prof1.add(0, "mh11CP-004", "C,G,G")
    prof1.add(1, "mh11CP-004", "C,G,G")
    prof1.add(0, "mh05KK-123", "A,C")
    prof1.add(1, "mh05KK-123", "A,T")
    prof2 = SimulatedProfile()
    prof2.add(0, "mh11CP-004", "C,T,A")
    prof2.add(1, "mh11CP-004", "C,T,G")
    prof2.add(0, "mh05KK-123", "A,T")
    prof2.add(1, "mh05KK-123", "A,T")
    prof3 = SimulatedProfile()
    prof3.add(0, "mh11CP-004", "C,G,G")
    prof3.add(1, "mh11CP-004", "T,G,G")
    prof3.add(0, "mh05KK-123", "G,C")
    prof3.add(1, "mh05KK-123", "G,T")
    profile = SimulatedProfile.merge([prof1, prof2, prof3])
    markers = pd.read_csv(data_file("def/loc2-offsets.tsv"), sep="\t")
    output = profile.bedstr(markers)
    print(output)
    assert output == (
        "mh05KK-123\t121\t122\tA|A|A|A|G|G\n"
        "mh05KK-123\t228\t229\tC|T|T|T|C|T\n"
        "mh11CP-004\t162\t163\tC|C|C|C|C|T\n"
        "mh11CP-004\t163\t164\tG|G|T|T|G|G\n"
        "mh11CP-004\t187\t188\tG|G|A|G|G|G\n"
    )


def test_bed_error():
    p = SimulatedProfile()
    p.add(0, "BOGUS", "A,C,C")
    p.add(1, "BOGUS", "A,C,C")
    markers = pd.read_csv(data_file("def/loc2-offsets.tsv"), sep="\t")
    with pytest.raises(ValueError, match=r"unknown marker identifier 'BOGUS'"):
        print(p.bedstr(markers))
