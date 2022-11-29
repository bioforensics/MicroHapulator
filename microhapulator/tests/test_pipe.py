# -------------------------------------------------------------------------------------------------
# Copyright (c) 2022, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from glob import glob
import microhapulator
import microhapulator.api as mhapi
from microhapulator.cli.pipe import validate_sample_input_files
from microhapulator.profile import SimulatedProfile, TypingResult
from microhapulator.tests import data_file
import pandas as pd
import pytest
from shutil import copyfile
from subprocess import run


def test_validate_sample_input_files():
    assert validate_sample_input_files(1, "S1", reads_are_paired=False) is True
    assert validate_sample_input_files(2, "S2", reads_are_paired=True) is True
    with pytest.raises(ValueError, match=r"expected 2 in paired-end mode"):
        validate_sample_input_files(1, "S3", reads_are_paired=True)
    with pytest.raises(ValueError, match=r"expected 1 in single-end mode"):
        validate_sample_input_files(2, "S4", reads_are_paired=False)
    with pytest.raises(FileNotFoundError):
        validate_sample_input_files(0, "S5")
    with pytest.raises(ValueError, match=r"expected 2 in paired-end mode"):
        validate_sample_input_files(4, "S6", reads_are_paired=True)


def test_pipe_missing_files(tmp_path):
    arglist = [
        "pipe",
        data_file("refr/usc10-refr.fna"),
        data_file("def/usc10-offsets.tsv"),
        data_file(""),
        "GBRusc",
        "--workdir",
        str(tmp_path),
        "--threads",
        "1",
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    with pytest.raises(FileNotFoundError, match=r"sample GBRusc: found 0 FASTQ files"):
        microhapulator.cli.pipe.main(args)


def test_pipe_gbr_usc10(tmp_path):
    arglist = [
        "pipe",
        data_file("refr/usc10-refr.fna"),
        data_file("def/usc10-offsets.tsv"),
        data_file(""),
        "gbr-usc",
        "--workdir",
        str(tmp_path),
        "--threads=1",
        "--static=5",
        "--dynamic=0.02",
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.pipe.main(args)
    expected = SimulatedProfile(fromfile=data_file("prof/gbr-usc10-sim.json"))
    observed = TypingResult(fromfile=tmp_path / "analysis" / "gbr-usc" / "gbr-usc-type.json")
    diff = list(mhapi.diff(observed, expected))
    assert len(diff) == 0
    report = tmp_path / "report.html"
    assert report.is_file()
    with open(report, "r") as fh:
        assert "Read lengths, uniform for all samples, are shown below" in fh.read()
    expected = pd.read_csv(data_file("gbr-usc-summary.tsv"), sep="\t")
    observed = pd.read_csv(tmp_path / "analysis" / "summary.tsv", sep="\t")
    assert observed.equals(expected)
    profile = tmp_path / "analysis" / "gbr-usc" / "profiles" / "gbr-usc-quant.csv"
    assert profile.is_file()
    expected = pd.read_csv(data_file("gbr-usc-profile.csv"))
    observed = pd.read_csv(profile)
    assert observed.equals(expected)
    call_pngs = glob(str(tmp_path / "analysis" / "*" / "callplots" / "*.png"))
    assert len(call_pngs) == 10
    donut_pngs = glob(str(tmp_path / "analysis" / "gbr-usc" / "gbr-usc-donut.png"))
    assert len(donut_pngs) == 1


def test_pipe_jpt_usc10_single(tmp_path):
    arglist = [
        "pipe",
        data_file("refr/usc10-refr.fna"),
        data_file("def/usc10-offsets.tsv"),
        data_file(""),
        "jpt-usc10",
        "--workdir",
        str(tmp_path),
        "--threads=1",
        "--static=5",
        "--dynamic=0.02",
        "--single",
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.pipe.main(args)
    expected = SimulatedProfile(fromfile=data_file("prof/jpt-usc10-sim.json"))
    observed = TypingResult(fromfile=tmp_path / "analysis" / "jpt-usc10" / "jpt-usc10-type.json")
    diff = list(mhapi.diff(observed, expected))
    assert len(diff) == 0
    report = tmp_path / "report.html"
    assert report.is_file()
    with open(report, "r") as fh:
        assert "Read Merging" not in fh.read()


@pytest.mark.parametrize(
    "filepath,samplename,plotmarker",
    [
        ("prof/gujarati-ind1-gt-filt.json", "Guj1", True),
        ("prof/gujarati-ind2-gt-filt.json", "Guj2", False),
        ("prof/gujarati-ind3-gt-filt.json", None, True),
    ],
)
def test_plot_haplotype_calls(filepath, samplename, plotmarker, tmp_path):
    result = TypingResult(fromfile=data_file(filepath))
    mhapi.plot_haplotype_calls(result, tmp_path, sample=samplename, plot_marker_name=plotmarker)
    pngs = glob(str(tmp_path / "*.png"))
    assert len(pngs) == 4
