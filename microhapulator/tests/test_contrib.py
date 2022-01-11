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
    "pjson,numcontrib",
    [
        ("prof/single-contrib-1.json", 1),
        ("prof/single-contrib-2.json", 1),
        ("prof/single-contrib-3.json", 1),
        ("prof/two-contrib-even.json", 2),
        ("prof/three-contrib-even.json", 3),
        ("prof/three-contrib-log.json", 3),
    ],
)
def test_contrib_json(pjson, numcontrib):
    profile = Profile(fromfile=data_file(pjson))
    n, *data = mhapi.contrib(profile)
    assert n == numcontrib


def test_contrib_main(capsys):
    arglist = ["contrib", data_file("prof/three-contrib-log.json")]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.contrib.main(args)
    terminal = capsys.readouterr()
    assert '"min_num_contrib": 3' in terminal.out
