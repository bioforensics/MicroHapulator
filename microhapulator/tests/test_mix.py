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

import microhapulator
from microhapulator.profile import SimulatedProfile
from microhapulator.tests import data_file
from tempfile import NamedTemporaryFile


def test_mix_main():
    with NamedTemporaryFile(suffix=".json.gz") as outfile:
        arglist = [
            "mix",
            "--out",
            outfile.name,
            data_file("prof/green-sim-gt-1.json.gz"),
            data_file("prof/green-sim-gt-2.json.gz"),
            data_file("prof/green-sim-gt-3.json.gz"),
        ]
        args = microhapulator.cli.get_parser().parse_args(arglist)
        microhapulator.cli.mix.main(args)
        p = SimulatedProfile(fromfile=outfile.name)
        testp = SimulatedProfile(fromfile=data_file("prof/green-sim-gt-combined.json.gz"))
        assert p == testp
