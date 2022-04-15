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


from argparse import ArgumentParser, RawDescriptionHelpFormatter
from microhapulator import __version__
from . import contain
from . import contrib
from . import convert
from . import diff
from . import dist
from . import filter
from . import hetbalance
from . import locbalance
from . import mix
from . import pipe
from . import prob
from . import seq
from . import sim
from . import type
from . import unite
import sys

mains = {
    "contain": contain.main,
    "contrib": contrib.main,
    "convert": convert.main,
    "diff": diff.main,
    "dist": dist.main,
    "filter": filter.main,
    "hetbalance": hetbalance.main,
    "locbalance": locbalance.main,
    "mix": mix.main,
    "pipe": pipe.main,
    "prob": prob.main,
    "seq": seq.main,
    "sim": sim.main,
    "type": type.main,
    "unite": unite.main,
}

subparser_funcs = {
    "contain": contain.subparser,
    "contrib": contrib.subparser,
    "convert": convert.subparser,
    "diff": diff.subparser,
    "dist": dist.subparser,
    "filter": filter.subparser,
    "hetbalance": hetbalance.subparser,
    "locbalance": locbalance.subparser,
    "mix": mix.subparser,
    "pipe": pipe.subparser,
    "prob": prob.subparser,
    "seq": seq.subparser,
    "sim": sim.subparser,
    "type": type.subparser,
    "unite": unite.subparser,
}


def get_parser():
    # https://patorjk.com/software/taag/, "Small" font
    bubbletext = r"""
≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
  __  __ _            _  _                _      _
 |  \/  (_)__ _ _ ___| || |__ _ _ __ _  _| |__ _| |_ ___ _ _
 | |\/| | / _| '_/ _ \ __ / _` | '_ \ || | / _` |  _/ _ \ '_|
 |_|  |_|_\__|_| \___/_||_\__,_| .__/\_,_|_\__,_|\__\___/_|
                               |_|
≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
Invoke `mhpl8r <subcmd> --help` and replace `<subcmd>` with one of the
subcommands listed below to see instructions for that operation.
"""
    subcommandstr = ", ".join(sorted(list(mains.keys())))
    parser = ArgumentParser(
        description=bubbletext,
        formatter_class=RawDescriptionHelpFormatter,
    )
    parser._positionals.title = "Subcommands"
    parser._optionals.title = "Global arguments"
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="MicroHapulator v{}".format(__version__),
    )
    subparsers = parser.add_subparsers(dest="subcmd", metavar="subcmd", help=subcommandstr)
    for func in subparser_funcs.values():
        func(subparsers)
    return parser


def main(arglist=None):  # pragma: no cover
    """Entry point for the MicroHapulator CLI.

    Isolated as a method so that the CLI can be called by other Python code
    (e.g. for testing), in which case the arguments are passed to the function.
    If no arguments are passed to the function, parse them from the command
    line.
    """
    args = get_parser().parse_args(arglist)
    if args.subcmd is None:
        get_parser().parse_args(["-h"])
    assert args.subcmd in mains
    mainmethod = mains[args.subcmd]
    print("[MicroHapulator] running version", __version__, file=sys.stderr)
    mainmethod(args)
