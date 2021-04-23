# -----------------------------------------------------------------------------
# Copyright (c) 2021, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import microhapulator
from microhapulator.balance import balance
from microhapulator.profile import Profile
from microhapulator.tests import data_file
import pandas
import pytest


def test_balance_basic(capfd):
    profile = Profile(fromfile=data_file('three-contrib-log.json'))
    obs_data = balance(profile)
    exp_data = pandas.read_csv(data_file('three-contrib-log-balance.csv'))
    assert obs_data.equals(exp_data)
    terminal = capfd.readouterr()
    assert 'MHDBL000212: ▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 50.00' in terminal.out


def test_balance_cli(tmp_path, capfd):
    outfile = str(tmp_path / 'balance.csv')
    arglist = ['balance', '--csv', outfile, data_file('three-contrib-log.json')]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.balance.main(args)
    obs_data = pandas.read_csv(outfile)
    exp_data = pandas.read_csv(data_file('three-contrib-log-balance.csv'))
    assert obs_data.equals(exp_data)
    terminal = capfd.readouterr()
    assert 'MHDBL000212: ▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 50.00' in terminal.out


def test_balance_cli_no_discard(capfd):
    arglist = ['balance', '--no-discarded', data_file('three-contrib-log.json')]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.balance.main(args)
    terminal = capfd.readouterr()
    print(terminal.out)
    assert 'MHDBL000212: ▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 40.00' in terminal.out
