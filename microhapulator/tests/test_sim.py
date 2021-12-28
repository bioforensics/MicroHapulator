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

import filecmp
import microhapulator
from microhapulator.profile import SimulatedProfile
from microhapulator.tests import data_file
import pandas as pd
import pytest


def test_meaning_of_life():
    freqs = pd.read_csv(data_file("freq/ceu50-freq.tsv"), sep="\t")
    observed = microhapulator.sim.sim(freqs, seed=42)
    expected = SimulatedProfile(fromfile=data_file("prof/meaning-of-life.json.gz"))
    assert observed == expected


def test_main(tmp_path):
    outfile = str(tmp_path / "profile.json")
    arglist = [
        "sim",
        "--out",
        outfile,
        "--seed",
        "1985",
        data_file("freq/ceu50-freq.tsv"),
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.sim.main(args)
    observed = SimulatedProfile(fromfile=outfile)
    expected = SimulatedProfile(fromfile=data_file("prof/bitusa-profile.json"))
    assert observed == expected


def test_no_seed():
    freqs = pd.read_csv(data_file("freq/asw2-freq.tsv"), sep="\t")
    genotype = microhapulator.sim.sim(freqs)
    assert len(genotype.data["markers"]) == 2
    assert sorted(genotype.data["markers"]) == ["mh07CP-004", "mh14CP-003"]


def test_main_haplo_seq(tmp_path):
    profile = str(tmp_path / "profile.json")
    hapseq = str(tmp_path / "haplo.fasta")
    arglist = [
        "sim",
        "--seed",
        "293847",
        "--out",
        profile,
        "--haplo-seq",
        hapseq,
        "--sequences",
        data_file("refr/orange-refr.fasta"),
        "--markers",
        data_file("def/orange-offsets.tsv"),
        data_file("freq/asw2-freq.tsv"),
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.sim.main(args)
    observed = SimulatedProfile(fromfile=profile)
    expected = SimulatedProfile(fromfile=data_file("prof/orange-sim-profile.json"))
    assert observed == expected
    assert filecmp.cmp(hapseq, data_file("orange-haplo.fasta"))
