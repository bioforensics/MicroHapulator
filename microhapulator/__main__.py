#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import microhapulator


def main(arglist=None):
    """Entry point for the MicroHapulator CLI.

    Isolated as a method so that the CLI can be called by other Python code
    (e.g. for testing), in which case the arguments are passed to the function.
    If no arguments are passed to the function, parse them from the command
    line.
    """
    args = microhapulator.cli.parse_args(arglist)
    if args.subcmd is None:  # pragma: no cover
        microhapulator.cli.parser().parse_args(['-h'])

    assert args.subcmd in microhapulator.cli.mains
    mainmethod = microhapulator.cli.mains[args.subcmd]
    versionmessage = '[MicroHapulator] running version {}'.format(microhapulator.__version__)
    microhapulator.plog(versionmessage)
    mainmethod(args)
