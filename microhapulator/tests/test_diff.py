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
from microhapulator.profile import SimulatedProfile
from microhapulator.tests import data_file
from tempfile import NamedTemporaryFile


def test_diff_basic():
    gt1 = SimulatedProfile(fromfile=data_file("prof/diff-comp-1.json"))
    gt2 = SimulatedProfile(fromfile=data_file("prof/diff-comp-2.json"))
    diff = list(microhapulator.op.diff(gt1, gt2))
    assert diff == [
        ("MHDBL000140", {"C,C,A,A"}, {"C,C,T,A"}),
        ("MHDBL000163", {"A,A,G,A,T"}, {"C,G,A,A,T"}),
    ]


def test_diff_large():
    gt1 = SimulatedProfile(fromfile=data_file("prof/diff-comp-1.json"))
    gt2 = SimulatedProfile(fromfile=data_file("prof/diff-comp-3.json"))
    diff = list(microhapulator.op.diff(gt1, gt2))
    loci = [d[0] for d in diff]
    print(diff[9], diff[17], diff[21])
    assert loci == [
        "MHDBL000002",
        "MHDBL000003",
        "MHDBL000007",
        "MHDBL000013",
        "MHDBL000017",
        "MHDBL000018",
        "MHDBL000030",
        "MHDBL000036",
        "MHDBL000038",
        "MHDBL000047",
        "MHDBL000058",
        "MHDBL000061",
        "MHDBL000076",
        "MHDBL000079",
        "MHDBL000082",
        "MHDBL000085",
        "MHDBL000088",
        "MHDBL000101",
        "MHDBL000106",
        "MHDBL000108",
        "MHDBL000111",
        "MHDBL000112",
        "MHDBL000122",
        "MHDBL000124",
        "MHDBL000128",
        "MHDBL000129",
        "MHDBL000135",
        "MHDBL000136",
        "MHDBL000138",
        "MHDBL000140",
        "MHDBL000144",
        "MHDBL000152",
        "MHDBL000154",
        "MHDBL000163",
        "MHDBL000181",
        "MHDBL000183",
        "MHDBL000194",
        "MHDBL000210",
        "MHDBL000211",
        "MHDBL000212",
    ]
    assert diff[9] == ("MHDBL000047", set(), {"T,T"})
    assert diff[17] == ("MHDBL000101", {"C,C,C,T"}, {"T,C,C,C"})
    assert diff[21] == ("MHDBL000112", {"G,G,A,C"}, set())


def test_diff_cli():
    f1 = data_file("prof/diff-comp-1.json")
    f2 = data_file("prof/diff-comp-3.json")
    with NamedTemporaryFile(suffix=".json") as outfile:
        arglist = ["diff", "-o", outfile.name, f1, f2]
        args = microhapulator.cli.get_parser().parse_args(arglist)
        microhapulator.cli.diff.main(args)
        with microhapulator.open(outfile.name, "r") as fh:
            output = fh.read().strip()
        with microhapulator.open(data_file("diff-comp-1-3.txt"), "r") as fh:
            testoutput = fh.read().strip()
        assert output == testoutput


def test_diff2():
    gt1 = SimulatedProfile(fromfile=data_file("prof/euramer-sim-gt.json"))
    gt2 = SimulatedProfile(fromfile=data_file("prof/euramer-inf-gt.json"))
    diff = list(microhapulator.op.diff(gt1, gt2))
    assert diff == [("MHDBL000018", set(), {"T,G,C,T,A"})]


def test_diff_nonmatching_alleles():
    p1 = SimulatedProfile(fromfile=data_file("prof/red-strict-profile.json"))
    p2 = SimulatedProfile(fromfile=data_file("prof/red-relaxed-profile.json"))
    diff = list(microhapulator.op.diff(p1, p2))
    print(diff)
    assert diff == [
        ("mh07CP-004", set(), {"T,T,T,A,T", "A,A,T,A,T"}),
        ("mh09KK-157", set(), {"G,C,C,A,T"}),
    ]
