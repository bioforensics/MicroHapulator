# -------------------------------------------------------------------------------------------------
# Copyright (c) 2018, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from microhapulator import happer
from microhapulator.happer.mutate import PloidyMismatchError
from microhapulator.happer.tests import data_file
import pytest


def test_populate_hap_index():
    with open(data_file("haplo-test1.bed"), "r") as alsfile:
        alleles = happer.allele.parse_alleles(alsfile)
        ploidy, haplotypes = happer.mutate.populate_haplotype_index(alleles)
    print(haplotypes)
    assert ploidy == 2
    assert len(haplotypes) == 1
    assert "chr17" in haplotypes
    assert len(haplotypes["chr17"]) == 2
    assert [len(hap) for hap in haplotypes["chr17"]] == [5, 5]


def test_empty_input():
    with pytest.raises(ValueError) as ve:
        ploidy, hapl = happer.mutate.populate_haplotype_index("")
    assert "no input provided to populate haplotype index" in str(ve)


def test_bad_ploidy():
    with open(data_file("ploidy-mismatch.bed"), "r") as alsfile:
        with pytest.raises(PloidyMismatchError) as pme:
            alleles = happer.allele.parse_alleles(alsfile)
            pldy, hapl = happer.mutate.populate_haplotype_index(alleles)
    assert "ploidy confusion" in str(pme)


def test_mutate_simple():
    with open(data_file("pico-refr.fa"), "r") as seqfile:
        seqdata = [data for data in happer.seqio.parse_fasta(seqfile)]
        sequences = [seq for label, seq in seqdata]
    with open(data_file("pico-hapl-1.bed"), "r") as allelefile:
        mutator = happer.mutate.mutate(seqdata, allelefile)
        haploseqs = [sequence for label, sequence in mutator]
    assert len(haploseqs) == 6
    assert "CAACCTTACGATCTA" in haploseqs[0]
    assert "CAACCTTCCGATCTA" in haploseqs[1]
    assert "GGGATAGACCCGTGG" in haploseqs[0]
    assert "GGGATAGGCCCGTGG" in haploseqs[1]
    assert haploseqs[2] == sequences[1]
    assert haploseqs[3] == sequences[1]
    assert "TGGACCGCATTGCAG" in haploseqs[4]
    assert "TGGACCTCATTGCAG" in haploseqs[5]


@pytest.mark.parametrize(
    "mainmethod",
    [
        happer.mutate.main,
        happer.__main__.main,
    ],
)
def test_mutate_cli(mainmethod, tmp_path):
    outfile = tmp_path / "out.fasta"
    arglist = ["--out", outfile, data_file("pico-refr.fa"), data_file("pico-hapl-1.bed")]
    arglist = map(str, arglist)
    args = happer.get_parser().parse_args(arglist)
    mainmethod(args)
    outseqs = list(happer.seqio.parse_fasta(open(outfile, "r")))
    deflines = [label for label, sequence in outseqs]
    sequences = [sequence for label, sequence in outseqs]
    assert len(deflines) == 6
    assert len(sequences) == 6
    assert deflines == [
        "chr1:hap1",
        "chr1:hap2",
        "chr2:hap1",
        "chr2:hap2",
        "chr3:hap1",
        "chr3:hap2",
    ]
    assert "CAACCTTACGATCTA" in sequences[0]
    assert "TGGACCTCATTGCAG" in sequences[5]
