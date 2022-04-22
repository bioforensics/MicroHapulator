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

import microhapulator
import microhapulator.api as mhapi
from microhapulator.profile import SimulatedProfile, TypingResult
from microhapulator.tests import data_file
import pandas as pd
import pytest
from shutil import copyfile
from subprocess import run


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
    with pytest.raises(FileNotFoundError, match=r"sample GBRusc: expected 2 FASTQ files, found 0"):
        microhapulator.cli.pipe.main(args)


def test_pipe_gbr_usc10(tmp_path):
    hg38 = str(tmp_path / "hg38-placeholder.fasta")
    copyfile(data_file("refr/usc10-refr.fna"), hg38)
    run(["bwa", "index", hg38])
    arglist = [
        "pipe",
        data_file("refr/usc10-refr.fna"),
        data_file("def/usc10-offsets.tsv"),
        data_file(""),
        "gbr-usc",
        "--workdir",
        str(tmp_path),
        "--threads",
        "1",
        "--hg38",
        hg38,
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.pipe.main(args)
    expected = SimulatedProfile(fromfile=data_file("prof/gbr-usc10-sim.json"))
    observed = TypingResult(fromfile=tmp_path / "analysis" / "gbr-usc" / "gbr-usc-type.json")
    diff = list(mhapi.diff(observed, expected))
    assert len(diff) == 0
    assert (tmp_path / "report.html").is_file()
    expected = pd.read_csv(data_file("gbr-usc-summary.tsv"), sep="\t")
    observed = pd.read_csv(tmp_path / "analysis" / "summary.tsv", sep="\t")
    assert observed.equals(expected)
