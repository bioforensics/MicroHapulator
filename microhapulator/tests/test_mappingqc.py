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
    figure = tmp_path / "donut.png"
    obs_data = mhapi.read_mapping_qc(
        data_file("mapping_workdir/analysis/SRM8398-2/SRM8398-2.bam.stats"),
        data_file(
            "mapping_workdir/analysis/SRM8398-2/fullrefr/SRM8398-2-fullrefr-mapped-reads.txt"
        ),
        data_file("mapping_workdir/analysis/SRM8398-2/SRM8398-2-repetitive-reads.csv"),
        figure,
    )
    exp_data = pd.read_csv(data_file("test-mapping-qc.csv"))
    print(obs_data)
    print(exp_data)
    assert obs_data.equals(exp_data)
    assert figure.is_file()


def test_mappingqc_cli(tmp_path):
    csv = tmp_path / "mapping-qc.csv"
    figure = tmp_path / "donut.png"
    arglist = [
        "mappingqc",
        data_file("mapping_workdir/analysis/SRM8398-2/SRM8398-2.bam.stats"),
        data_file(
            "mapping_workdir/analysis/SRM8398-2/fullrefr/SRM8398-2-fullrefr-mapped-reads.txt"
        ),
        data_file("mapping_workdir/analysis/SRM8398-2/SRM8398-2-repetitive-reads.csv"),
        "--csv",
        str(csv),
        "--figure",
        str(figure),
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.mappingqc.main(args)
    obs_data = pd.read_csv(csv)
    exp_data = pd.read_csv(data_file("test-mapping-qc.csv"))
    assert obs_data.equals(exp_data)
    assert figure.is_file()
