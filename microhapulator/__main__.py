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


def main(arglist=None):  # pragma: no cover
    """Entry point for the MicroHapulator CLI.

    Isolated as a method so that the CLI can be called by other Python code
    (e.g. for testing), in which case the arguments are passed to the function.
    If no arguments are passed to the function, parse them from the command
    line.
    """
    args = microhapulator.cli.parse_args(arglist)
    if args.subcmd is None:
        microhapulator.cli.get_parser().parse_args(["-h"])

    assert args.subcmd in microhapulator.cli.mains
    mainmethod = microhapulator.cli.mains[args.subcmd]
    versionmessage = "[MicroHapulator] running version {}".format(microhapulator.__version__)
    microhapulator.plog(versionmessage)
    mainmethod(args)
