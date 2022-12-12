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


def test_mappingqc(tmp_path):
    marker_counts = data_file("test-mapped-reads.txt")
    full_refr_counts = data_file("test-fullrefr-mapped-reads.txt")
    rep_counts = data_file("repetitive-reads.csv")
    figure = tmp_path / "donut.png"
    obs_data = mhapi.read_mapping_qc(marker_counts, full_refr_counts, rep_counts, figure)
    exp_data = pd.read_csv(data_file("test-mapping-qc.csv"))
    assert obs_data.equals(exp_data)
    with open(marker_counts) as fh:
        total_reads = int(next(fh).split()[2])
    assert total_reads == obs_data.sum(axis="columns")[0]
    assert len(obs_data.columns) == 4
    assert figure.is_file()


def test_no_refr_offsets(tmp_path):
    marker_counts = data_file("test-mapped-reads.txt")
    full_refr_counts = data_file("test-fullrefr-mapped-reads.txt")
    rep_counts = data_file("repetitive-reads-empy.csv")
    figure = tmp_path / "donut-no-rerf-offsets.png"
    obs_data = mhapi.read_mapping_qc(marker_counts, full_refr_counts, rep_counts, figure)
    exp_data = pd.read_csv(data_file("mapping-qc-no-offsets.csv"))
    assert obs_data.equals(exp_data)
    with open(marker_counts) as fh:
        total_reads = int(next(fh).split()[2])
    assert total_reads == obs_data.sum(axis="columns")[0]
    assert len(obs_data.columns) == 3
    assert figure.is_file()


def test_mappingqc_cli(tmp_path):
    arglist = [
        "mappingqc",
        "--marker",
        data_file("test-mapped-reads.txt"),
        "--refr",
        data_file("test-fullrefr-mapped-reads.txt"),
        "--rep",
        data_file("repetitive-reads.csv"),
        "--csv",
        str(tmp_path / "mapping-qc.csv"),
        "--figure",
        str(tmp_path / "donut.png"),
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.mappingqc.main(args)
    obs_data = pd.read_csv(str(tmp_path / "mapping-qc.csv"))
    exp_data = pd.read_csv(data_file("test-mapping-qc.csv"))
    assert obs_data.equals(exp_data)
