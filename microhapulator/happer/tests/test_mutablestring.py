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

from microhapulator.happer import MutableString
import pytest


@pytest.mark.parametrize(
    "thestring", [("GATTACA"), ("ATGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGA"), ("MAHERSHALALHASHBAZ")]
)
def test_basic(thestring):
    ms = MutableString(thestring)
    assert str(ms) == thestring
    assert repr(ms) == thestring
    assert ms == thestring
    assert len(ms) == len(thestring)
    assert ms[3] == thestring[3]
    assert ms[2:6] == thestring[2:6]


def test_add_ops():
    ms = MutableString("HOWGOESIT")
    assert ms + "MAN" == "HOWGOESITMAN"
    assert ms + MutableString("MAN") == "HOWGOESITMAN"

    ms += "MANS!"
    assert ms == "HOWGOESITMANS!"
    assert "ITMANS" in ms

    ms[9:9] = "FELLOWHU"
    assert ms == "HOWGOESITFELLOWHUMANS!"


def test_del_ops():
    ms = MutableString("Sesame seed buns")
    del ms[15]
    assert ms == "Sesame seed bun"

    del ms[7:12]
    assert ms == "Sesame bun"
