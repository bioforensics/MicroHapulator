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
from microhapulator.profile import SimulatedProfile
from microhapulator.tests import data_file
import pytest
import shutil
import tempfile


def test_meaning_of_life():
    panel = [
        "mh01KK-205",
        "mh01CP-016",
        "mh01KK-117",
        "mh10CP-003",
        "mh10KK-170",
        "mh10KK-101",
        "mh11KK-180",
        "mh11KK-037",
        "mh11KK-191",
        "mh12CP-008",
        "mh12KK-202",
        "mh12KK-046",
        "mh13KK-213",
        "mh13KK-223",
        "mh14CP-003",
        "mh14KK-048",
        "mh15KK-069",
        "mh15KK-104",
        "mh16KK-255",
        "mh17CP-001",
        "mh17KK-052",
        "mh17KK-055",
        "mh18CP-005",
        "mh18KK-293",
        "mh19KK-299",
        "mh19KK-057",
        "mh02KK-138",
        "mh02KK-136",
        "mh20KK-307",
        "mh20KK-058",
        "mh21KK-315",
        "mh22KK-060",
        "mh22KK-061",
        "mh03CP-005",
        "mh03KK-007",
        "mh03KK-150",
        "mh04KK-030",
        "mh04KK-013",
        "mh04KK-017",
        "mh05KK-020",
        "mh06CP-007",
        "mh06KK-008",
        "mh07CP-004",
        "mh07KK-030",
        "mh07KK-031",
        "mh08KK-039",
        "mh08KK-032",
        "mh09KK-033",
        "mh09KK-153",
        "mh09KK-157",
    ]
    profile = microhapulator.sim.sim(["SA004250L", "SA004250L"], panel, seed=42)
    testprofile = SimulatedProfile(fromfile=data_file("meaning-of-life.json.gz"))
    assert profile == testprofile


@pytest.mark.parametrize(
    "relaxmode,testfile",
    [
        (False, "red-strict-profile.json"),
        (True, "red-relaxed-profile.json"),
    ],
)
def test_sim_relaxed(relaxmode, testfile):
    profile = microhapulator.sim.sim(
        ["SA000101C"], ["mh01KK-117", "mh09KK-157", "mh07CP-004"], seed=54321, relaxed=relaxmode
    )
    testprofile = SimulatedProfile(fromfile=data_file(testfile))
    assert profile == testprofile


def test_main():
    with tempfile.NamedTemporaryFile(suffix="-profile.json") as outfile:
        arglist = [
            "sim",
            "--out",
            outfile.name,
            "--seed",
            "1985",
            "SA004250L",
            "SA004250L",
            "mh13KK-218",
            "mh05KK-170",
            "mh21KK-320",
            "mh10KK-163",
            "mh10KK-169",
            "mh02KK-134",
            "mh16KK-049",
            "mh06KK-008",
            "mh21KK-324",
            "mh11KK-180",
            "mh13KK-217",
            "mh21KK-315",
            "mh04KK-030",
            "mh01KK-117",
            "mh19KK-299",
            "mh01KK-205",
            "mh13KK-223",
            "mh02KK-136",
            "mh04KK-013",
            "mh13KK-213",
            "mh16KK-255",
            "mh20KK-307",
            "mh09KK-157",
            "mh13KK-225",
            "mh22KK-061",
            "mh18KK-293",
            "mh03KK-150",
            "mh01KK-001",
            "mh21KK-316",
            "mh11KK-191",
            "mh01NK-001",
            "mh12KK-202",
            "mh09KK-153",
            "mh17KK-272",
            "mh16KK-302",
            "mh09KK-152",
            "mh10KK-170",
            "mh09KK-033",
            "mh20KK-058",
            "mh11KK-187",
            "mh04KK-017",
            "mh01KK-172",
            "mh06KK-030",
            "mh18KK-285",
            "mh01KK-106",
            "mh06KK-025",
            "mh02KK-004",
            "mh09KK-020",
            "mh07KK-030",
            "mh02KK-213",
        ]
        args = microhapulator.cli.get_parser().parse_args(arglist)
        microhapulator.sim.main(args)
        p = SimulatedProfile(fromfile=outfile.name)
        testp = SimulatedProfile(fromfile=data_file("bitusa-profile.json"))
        assert p == testp


def test_main_haplo_seq(tmp_path):
    profile = str(tmp_path / "profile.json")
    hapseq = str(tmp_path / "haplo.fasta")
    arglist = [
        "sim",
        "--seed",
        "293847",
        "--out",
        profile,
        "--haplo-seq",
        hapseq,
        "SA004047P",
        "SA004047P",
        "mh07CP-004",
        "mh14CP-003",
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.sim.main(args)
    p = SimulatedProfile(fromfile=profile)
    testp = SimulatedProfile(fromfile=data_file("orange-sim-profile.json"))
    assert p == testp
    assert filecmp.cmp(hapseq, data_file("orange-haplo.fasta"))


def test_no_seed():
    genotype = microhapulator.sim.sim(["SA004047P"], ["mh07CP-004", "mh14CP-003"])
    assert len(genotype.data["markers"]) == 2
    assert sorted(genotype.data["markers"]) == ["mh07CP-004", "mh14CP-003"]


def test_bad_panel():
    with pytest.raises(ValueError, match=r"invalid panel") as ve:
        microhapulator.sim.sim(["SA004047P"], ["DUUUUDE", "SWEEEET"])


def test_panel_file():
    arglist = ["sim", "Pashtun", "Pashtun", data_file("minipanel.txt")]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.sim.main(args)
