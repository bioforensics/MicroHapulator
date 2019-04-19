#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from io import StringIO
import microhapdb
import microhapulator
import numpy
import pytest


def test_sim_genotype_roundtrip():
    output = StringIO()
    panel = ['MHDBL000011', 'MHDBL000091', 'MHDBL000171', 'MHDBL000159', 'MHDBL000069']
    populations = ['MHDBP000004', 'MHDBP000022']
    genotype = microhapulator.genotype.SimulatedGenotype()
    simulator = microhapulator.panel.sample_panel(populations, panel)
    for haplotype, locusid, allele in simulator:
        genotype.add(haplotype, locusid, allele)
    print(genotype, file=output)
    bedlines = output.getvalue().split('\n')
    testgenotype = microhapulator.genotype.SimulatedGenotype(frombed=bedlines)
    assert testgenotype == genotype
