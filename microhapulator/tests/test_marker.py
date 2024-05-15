# -------------------------------------------------------------------------------------------------
# Copyright (c) 2023, DHS.
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
from microhapulator.marker import MicrohapIndex
from microhapulator.tests import data_file
import pytest


def test_index_basic():
    defn = data_file("def/2x2.tsv")
    refr = data_file("refr/2x2.fasta")
    index = MicrohapIndex.from_files(defn, refr)
    loci = [locus.id for locus in index.loci.values()]
    assert loci == ["mh21FHL-002", "mh21KK-320"]
    marker = index.markers["mh21KK-320.v1"]
    assert marker.offsets_locus == [10, 80, 169, 195]
    assert index.sequence("mh21KK-320").startswith("TAGGATGGGCGACCTTTCCTGTGGGCTAAGGTAGGAAAGCAGAAA")


@pytest.mark.parametrize(
    "tsv,nloci,locusid,offsets",
    [
        ("def/acb-dozen-offsets.tsv", 12, "mh08USC-8qC", [137028438, 137028480, 137028499]),
        ("def/yellow-offsets.tsv", 5, "mh01KK-106", [4167403, 4167500, 4167563]),
    ],
)
def test_index_loading(tsv, nloci, locusid, offsets):
    microhaps = MicrohapIndex.from_files(data_file(tsv))
    assert len(microhaps.loci) == nloci
    mhit = iter(microhaps)
    locus, marker = next(mhit)
    locus, marker = next(mhit)
    assert locus.id == locusid
    print(marker.offsets_chrom[:3])
    assert marker.offsets_chrom[:3] == offsets


def test_load_marker_definitions_missing_column():
    message = r"column\(s\) missing from marker definition file: Marker"
    with pytest.raises(ValueError, match=message):
        MicrohapIndex.from_files(data_file("def/orange-offsets-missing.tsv"))


def test_load_marker_definitions_minimal():
    defn_file = data_file("def/default-panel-offsets-minimal.tsv")
    with pytest.raises(ValueError, match="missing from marker definition file: Chrom, OffsetHg3"):
        MicrohapIndex.from_files(defn_file)


@pytest.mark.parametrize(
    "tsv,fasta,nseq,firstseq",
    [
        ("def/acb-dozen-offsets.tsv", "refr/acb-dozen-refr.fasta", 12, "mh07USC-7qB"),
        ("def/russ4-offsets.tsv", "refr/russ4-refr.fasta.gz", 4, "mh01KK-001"),
    ],
)
def test_load_marker_reference_sequences(tsv, fasta, nseq, firstseq):
    microhaps = MicrohapIndex.from_files(data_file(tsv), data_file(fasta))
    assert len(microhaps.sequences) == nseq
    observed = next(iter(microhaps.sequences.keys()))
    print(firstseq, observed)
    assert firstseq == observed
