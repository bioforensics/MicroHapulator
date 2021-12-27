# -------------------------------------------------------------------------------------------------
# Copyright (c) 2021, DHS.
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
from microhapulator.balance import balance
from microhapulator.profile import Profile
from microhapulator.tests import data_file
import pandas
import pytest


def test_balance_basic(capfd):
    profile = Profile(fromfile=data_file("prof/three-contrib-log.json"))
    obs_data = balance(profile)
    exp_data = pandas.read_csv(data_file("three-contrib-log-balance.csv"))
    assert obs_data.equals(exp_data)
    terminal = capfd.readouterr()
    assert "MHDBL000212: ▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 50.00" in terminal.out


def test_balance_cli(tmp_path, capfd):
    outfile = str(tmp_path / "balance.csv")
    arglist = ["balance", "--csv", outfile, data_file("prof/three-contrib-log.json")]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.balance.main(args)
    obs_data = pandas.read_csv(outfile)
    exp_data = pandas.read_csv(data_file("three-contrib-log-balance.csv"))
    assert obs_data.equals(exp_data)
    terminal = capfd.readouterr()
    assert "MHDBL000212: ▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 50.00" in terminal.out


def test_balance_cli_no_discard(capfd):
    arglist = ["balance", "--no-discarded", data_file("prof/three-contrib-log.json")]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.balance.main(args)
    terminal = capfd.readouterr()
    print(terminal.out)
    assert "MHDBL000212: ▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 40.00" in terminal.out
