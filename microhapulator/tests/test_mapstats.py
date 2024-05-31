# -------------------------------------------------------------------------------------------------
# Copyright (c) 2024, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from microhapulator.mapstats import MappingStats, MappingSummary
from microhapulator.tests import data_file
import pytest


def test_load_mapstats_from_workdir():
    stats = MappingStats.from_workdir("SRM8398-2", workdir=data_file("mapping_workdir"))
    assert stats.total == 123335
    assert stats.mapped == 122058
    assert stats.chisq == pytest.approx(0.1421)
    assert stats.total_reads == "123,335"
    assert stats.mapped_reads == "122,058"
    assert stats.mapping_rate == "99.0%"
    assert stats.chi_square == "0.142"


def test_load_mapsummary_from_workdir():
    samples = ("SRM8398-1", "SRM8398-2", "SRM8398-3")
    summary = MappingSummary.from_workdir(samples, workdir=data_file("mapping_workdir"))
    assert len(summary.markers) == 29
    repetitive = summary.repetitive_reads_by_marker()
    assert repetitive["mh10KK-162.v3"]["SRM8398-1"].repetitive == 108
    assert repetitive["mh10KK-162.v3"]["SRM8398-2"].repetitive == 92
    assert repetitive["mh10KK-162.v3"]["SRM8398-3"].repetitive == 107
    assert repetitive["mh17KK-278.v1"]["SRM8398-1"].repetitive == 0
    assert repetitive["mh17KK-278.v1"]["SRM8398-2"].repetitive == 0
    assert repetitive["mh17KK-278.v1"]["SRM8398-3"].repetitive == 0
