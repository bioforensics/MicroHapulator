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
import microhapulator.api as mhapi
from microhapulator.marker import MicrohapIndex
from microhapulator.profile import Profile
from microhapulator.tests import data_file
import numpy.random
import os
import pytest
from tempfile import NamedTemporaryFile


@pytest.mark.parametrize(
    "n,totalreads,prop,result",
    [
        (3, 100, None, [33, 33, 33]),
        (3, 100, [0.5, 0.4, 0.1], [50, 40, 10]),
        (4, 4000000, [20, 30, 40, 50], [571428, 857142, 1142857, 1428571]),
    ],
)
def test_proportions(n, totalreads, prop, result):
    assert mhapi.calc_n_reads_from_proportions(n, totalreads, prop) == result


def test_proportions_failure_modes():
    message = r"mismatch between contributor number and proportions"
    with pytest.raises(ValueError, match=message):
        mhapi.calc_n_reads_from_proportions(3, 1000, [0.6, 0.4])


def test_even_mixture():
    seed = numpy.random.randint(1, 2**32 - 1)
    print("Seed:", seed)
    numpy.random.seed(seed)
    freqs = microhapulator.load_marker_frequencies(data_file("freq/acb-dozen-freq.tsv"))
    index = MicrohapIndex.from_files(
        data_file("def/acb-dozen-offsets.tsv"), fasta_path=data_file("refr/acb-dozen-refr.fasta")
    )
    profiles = list()
    for _ in range(numpy.random.randint(2, 6)):
        p = mhapi.sim(freqs)
        profiles.append(p)
    sequencer = mhapi.seq(profiles, index, totalreads=1000)
    for n, read1, read2 in sequencer:
        pass
    numfragments = n * 2
    assert numfragments == pytest.approx(1000, abs=50)


def test_complex_genotype(capsys):
    profile = Profile(fromfile=data_file("prof/mixture-genotype.json"))
    index = MicrohapIndex.from_files(
        data_file("def/russ4-offsets.tsv"), fasta_path=data_file("refr/russ4-refr.fasta.gz")
    )
    sequencer = mhapi.seq(list(profile.unmix()), index, totalreads=200)
    for n, read in enumerate(sequencer):
        pass
    terminal = capsys.readouterr()
    assert terminal.err.count("Individual seed=") == 3


def test_uneven_mixture(capsys):
    freqs = microhapulator.load_marker_frequencies(data_file("freq/russ4-freq.tsv"))
    index = MicrohapIndex.from_files(
        data_file("def/russ4-offsets.tsv"), fasta_path=data_file("refr/russ4-refr.fasta.gz")
    )
    profiles = [mhapi.sim(freqs) for _ in range(3)]
    sequencer = mhapi.seq(profiles, index, totalreads=500, proportions=[0.5, 0.3, 0.2])
    for read in sequencer:
        pass
    terminal = capsys.readouterr()
    assert "numreads=250" in terminal.err
    assert "numreads=150" in terminal.err
    assert "numreads=100" in terminal.err


def test_mixture_failure_modes():
    freqs = microhapulator.load_marker_frequencies(data_file("freq/russ4-freq.tsv"))
    index = MicrohapIndex.from_files(
        data_file("def/russ4-offsets.tsv"), fasta_path=data_file("refr/russ4-refr.fasta.gz")
    )
    profiles = [mhapi.sim(freqs) for _ in range(3)]
    with pytest.raises(ValueError, match=r"number of profiles must match number of seeds"):
        for read in mhapi.seq(profiles, index, seeds=[42, 1776]):
            pass
    with pytest.raises(ValueError, match=r"mismatch between contributor number and proportions"):
        for read in mhapi.seq(profiles, index, proportions=[0.5, 0.3, 0.1, 0.1]):
            pass
    message = r"specified proportions result in 0 reads for 1 or more individuals"
    with pytest.raises(ValueError, match=message):
        for read in mhapi.seq(profiles, index, totalreads=500, proportions=[1, 100, 10000]):
            pass


def test_main(tmp_path):
    outfile = str(tmp_path / "reads.fastq")
    arglist = [
        "seq",
        "--out",
        outfile,
        "--seeds",
        "123454321",
        "--num-reads",
        "500",
        "--signature",
        "srd6Sei",
        data_file("def/orange-offsets.tsv"),
        data_file("refr/orange-refr.fasta"),
        data_file("prof/orange-sim-profile.json"),
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.seq.main(args)
    assert filecmp.cmp(outfile, data_file("orange-reads.fastq"))


def test_main_no_seed():
    with NamedTemporaryFile(suffix=".fastq") as outfile:
        arglist = [
            "seq",
            "--out",
            outfile.name,
            "--num-reads",
            "200",
            data_file("def/orange-offsets.tsv"),
            data_file("refr/orange-refr.fasta"),
            data_file("prof/orange-sim-profile.json"),
        ]
        args = microhapulator.cli.get_parser().parse_args(arglist)
        microhapulator.cli.seq.main(args)
        with open(outfile.name, "r") as fh:
            filelines = fh.read().strip().split("\n")
            assert len(filelines) == 800  # 200 reads * 4 lines per read = 800 lines


def test_main_mixture(capsys):
    arglist = [
        "seq",
        "--seeds",
        "42",
        "1776",
        "--proportions",
        "0.8",
        "0.2",
        "--num-reads",
        "500",
        data_file("def/yellow-offsets.tsv"),
        data_file("refr/yellow-refr.fasta.gz"),
        data_file("prof/yellow-mix-gt.json"),
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.seq.main(args)
    terminal = capsys.readouterr()
    outlines = terminal.out.strip().split("\n")
    nrecords = len(outlines) / 4
    assert nrecords == pytest.approx(500, abs=25)
    assert outlines[-3] == (
        "TCAATTCAATTTCTACCCTCAGCATCAAGGCAGGGGTTCATCATAATGGGTATTGGAGGCTCAAAGAAA"
        "ATTTAGGCTCAGCACACACACACACACACACACACACACACAGCGATTTTTAATGCTGGTACAATCACA"
        "GGAGACTGCGACCCAGCCCTCCTCAGCGCCTCGGGTGCTCACGGGCACTCCTGGAGTCTCGGCCACACT"
        "AAGTCCCCCTGGTGGCCACACAGAAGAAGAGGTGGTAAAACTTTCTGGGAGTGAGATCAAAAATTTTAG"
        "GAGTCTAAAAACATACTTTTCTAAG"
    )


def test_main_out_stdout(capsys):
    arglist = [
        "seq",
        "--num-reads",
        "100",
        data_file("def/orange-offsets.tsv"),
        data_file("refr/orange-refr.fasta"),
        data_file("prof/orange-sim-profile.json"),
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.seq.main(args)
    terminal = capsys.readouterr()
    outlines = terminal.out.strip().split("\n")
    nrecords = len(outlines) / 4
    assert nrecords == pytest.approx(100, abs=5)


def test_main_out_one_filename(tmp_path):
    outfile = str(tmp_path / "reads-interleaved.fastq")
    arglist = [
        "seq",
        "--out",
        outfile,
        "--num-reads",
        "100",
        data_file("def/orange-offsets.tsv"),
        data_file("refr/orange-refr.fasta"),
        data_file("prof/orange-sim-profile.json"),
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.seq.main(args)
    assert os.path.isfile(outfile)
    with open(outfile, "r") as fh:
        outlines = fh.read().strip().split("\n")
        nrecords = len(outlines) / 4
        assert nrecords == pytest.approx(100, abs=5)


def test_main_out_two_filenames(tmp_path):
    f1 = str(tmp_path / "reads-R1.fastq")
    f2 = str(tmp_path / "reads-R2.fastq")
    arglist = [
        "seq",
        "--out",
        f1,
        f2,
        "--num-reads",
        "100",
        data_file("def/orange-offsets.tsv"),
        data_file("refr/orange-refr.fasta"),
        data_file("prof/orange-sim-profile.json"),
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.seq.main(args)
    assert os.path.isfile(f1)
    assert os.path.isfile(f2)
    for fn in (f1, f2):
        with open(fn, "r") as fh:
            outlines = fh.read().strip().split("\n")
            nrecords = len(outlines) / 4
            assert nrecords == pytest.approx(50, abs=5)


def test_main_out_three_filenames(capsys):
    arglist = [
        "seq",
        "--out",
        "ONE",
        "TWO",
        "THREE",
        "--",
        data_file("def/orange-offsets.tsv"),
        data_file("refr/orange-refr.fasta"),
        data_file("prof/orange-sim-profile.json"),
    ]
    with pytest.raises(SystemExit):
        _ = microhapulator.cli.get_parser().parse_args(arglist)
    terminal = capsys.readouterr()
    assert "expected 1 or 2 output filenames, got 3" in terminal.err


def test_regression_seq_locus_names_in_refr():
    profile = Profile(fromfile=data_file("prof/mwg4-sim.json"))
    index = MicrohapIndex.from_files(
        data_file("def/mwg4-offsets.tsv"), data_file("refr/mwg4-refr.fasta")
    )
    sequencer = mhapi.seq([profile], index, totalreads=1000)
    for n, read1, read2 in sequencer:
        pass
    numfragments = n * 2
    assert numfragments == pytest.approx(1000, abs=50)


def test_seq_multidef():
    profile = Profile(fromfile=data_file("pashtun-sim/tiny-panel-multidef-bogus-profile.json"))
    index = MicrohapIndex.from_files(
        data_file("pashtun-sim/tiny-panel-multidef.tsv"),
        data_file("pashtun-sim/tiny-panel-multidef.fasta"),
    )
    sequencer = mhapi.seq([profile], index, totalreads=1000)
    msg = "simulating sequencing of profiles with multiple allele definitions at a locus is not permitted"
    with pytest.raises(ValueError, match=msg):
        _ = list(sequencer)
