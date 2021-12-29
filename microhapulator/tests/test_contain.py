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
import microhapulator.api as mhapi
from microhapulator.profile import Profile
from microhapulator.tests import data_file
import pytest


@pytest.mark.parametrize(
    "f1,f2,total,contained",
    [
        ("prof/four-brits-sim.json", "prof/one-brit-sim.json", 38, 38),
        ("prof/four-brits-sim.json", "prof/one-american-sim.json", 38, 31),
        ("prof/one-brit-sim.json", "prof/one-american-sim.json", 38, 11),
    ],
)
def test_contain(f1, f2, total, contained):
    profile1 = Profile(data_file(f1))
    profile2 = Profile(data_file(f2))
    c, t = mhapi.contain(profile1, profile2)
    assert t == total
    assert c == contained


def test_contain_cli(capsys):
    arglist = [
        "contain",
        data_file("prof/one-brit-sim.json"),
        data_file("prof/one-italian-sim.json"),
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.contain.main(args)
    terminal = capsys.readouterr()
    assert '"containment": 0.4444' in terminal.out
