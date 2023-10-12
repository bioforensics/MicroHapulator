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
import microhapulator.api as mhapi
from microhapulator.profile import Profile
from microhapulator.tests import data_file
import pandas as pd
import pytest
import sys


def test_interlocus_balance_basic(capfd):
    profile = Profile(fromfile=data_file("prof/three-contrib-log.json"))
    chisq, obs_data = mhapi.interlocus_balance(profile)
    exp_data = pd.read_csv(data_file("three-contrib-log-balance.csv"))
    exp_data.to_csv(sys.stdout, index=False)
    obs_data.to_csv(sys.stdout, index=False)
    assert obs_data.equals(exp_data)
    assert chisq == pytest.approx(0.00928395)
    terminal = capfd.readouterr()
    assert "MHDBL000212: ▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 50.00" in terminal.out


def test_locbalance_cli(tmp_path, capfd):
    outfile = str(tmp_path / "balance.csv")
    arglist = ["locbalance", "--csv", outfile, data_file("prof/three-contrib-log.json")]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.locbalance.main(args)
    obs_data = pd.read_csv(outfile)
    exp_data = pd.read_csv(data_file("three-contrib-log-balance.csv"))
    exp_data.to_csv(sys.stdout, index=False)
    obs_data.to_csv(sys.stdout, index=False)
    assert obs_data.equals(exp_data)
    terminal = capfd.readouterr()
    print(terminal.out)
    assert "Extent of imbalance (chi-square statistic): 0.0093" in terminal.out
    assert "MHDBL000212: ▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 50.00" in terminal.out


def test_locbalance_cli_no_discard(capfd):
    arglist = ["locbalance", "--no-discarded", data_file("prof/three-contrib-log.json")]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.locbalance.main(args)
    terminal = capfd.readouterr()
    print(terminal.out)
    assert "Extent of imbalance (chi-square statistic): 0.0221" in terminal.out
    assert "MHDBL000212: ▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 40.00" in terminal.out


def test_heterozygote_balance_basic(tmp_path):
    figfile = tmp_path / "figure.png"
    profile = Profile(fromfile=data_file("prof/single-contrib-2.json"))
    tstat, obs_data = mhapi.heterozygote_balance(profile, tofile=figfile)
    assert tstat == pytest.approx(3.90845)
    exp_data = pd.read_csv(data_file("het-balance.tsv"), sep="\t")
    assert obs_data.equals(exp_data)
    assert figfile.is_file()
