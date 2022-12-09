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
from microhapulator.tests import data_file
import pandas as pd


def test_repetitive():
    marker_bam = data_file("bam/repetitive.bam")
    fullrefr_bam = data_file("bam/repetitive-fullrefr.bam")
    marker_def = data_file("def/repetitive-offsets.tsv")
    obs_data = mhapi.repetitive_mapping(marker_bam, fullrefr_bam, marker_def)
    exp_data = pd.read_csv(data_file("repetitive-reads.csv"))
    assert obs_data.equals(exp_data)


def test_repetitive_cli(tmp_path):
    outfile = str(tmp_path / "repetitive-reads.csv")
    arglist = [
        "repetitive",
        data_file("bam/repetitive.bam"),
        data_file("bam/repetitive-fullrefr.bam"),
        data_file("def/repetitive-offsets.tsv"),
        "--out",
        outfile,
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.repetitive.main(args)
    obs_data = pd.read_csv(outfile)
    exp_data = pd.read_csv(data_file("repetitive-reads.csv"))
    assert obs_data.equals(exp_data)
