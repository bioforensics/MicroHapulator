#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from argparse import ArgumentParser, RawDescriptionHelpFormatter
import microhapulator
import sys
from . import contain
from . import contrib
from . import diff
from . import dist
from . import mix
from . import prob
from . import seq
from . import sim
from . import type
from . import unite

mains = {
    'contain': microhapulator.contain.main,
    'contrib': microhapulator.contrib.main,
    'diff': microhapulator.diff.main,
    'dist': microhapulator.dist.main,
    'mix': microhapulator.mix.main,
    'prob': microhapulator.prob.main,
    'seq': microhapulator.seq.main,
    'sim': microhapulator.sim.main,
    'type': microhapulator.type.main,
    'unite': microhapulator.unite.main,
}

subparser_funcs = {
    'contain': contain.subparser,
    'contrib': contrib.subparser,
    'diff': diff.subparser,
    'dist': dist.subparser,
    'mix': mix.subparser,
    'prob': prob.subparser,
    'seq': seq.subparser,
    'sim': sim.subparser,
    'type': type.subparser,
    'unite': unite.subparser,
}


def get_parser():
    # https://patorjk.com/software/taag/, "Small" font
    bubbletext = r'''
≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
  __  __ _            _  _                _      _
 |  \/  (_)__ _ _ ___| || |__ _ _ __ _  _| |__ _| |_ ___ _ _
 | |\/| | / _| '_/ _ \ __ / _` | '_ \ || | / _` |  _/ _ \ '_|
 |_|  |_|_\__|_| \___/_||_\__,_| .__/\_,_|_\__,_|\__\___/_|
                               |_|
≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
Invoke `mhpl8r <subcmd> --help` and replace `<subcmd>` with one of the
subcommands listed below to see instructions for that operation.
'''
    subcommandstr = ', '.join(sorted(list(mains.keys())))
    parser = ArgumentParser(
        description=bubbletext,
        formatter_class=RawDescriptionHelpFormatter,
    )
    parser._positionals.title = 'Subcommands'
    parser._optionals.title = 'Global arguments'
    parser.add_argument('-v', '--version', action='version',
                        version='MicroHapulator v{}'.format(microhapulator.__version__))
    parser.add_argument('-l', '--logfile', metavar='F', help='log file for '
                        'diagnostic messages, warnings, and errors')
    parser.add_argument('--tee', action='store_true', help='write diagnostic '
                        'output to logfile AND terminal (stderr)')
    subparsers = parser.add_subparsers(dest='subcmd', metavar='subcmd',
                                       help=subcommandstr)
    for func in subparser_funcs.values():
        func(subparsers)
    return parser


def parse_args(arglist=None):  # pragma: no cover
    '''Parse arguments from the command line and configure logging.

    Calling this function in a testing context has proven problematic. Please
    use only in a runtime context.
    '''
    args = get_parser().parse_args(arglist)
    microhapulator.logstream = sys.stderr
    if args.logfile and args.logfile != '-':
        microhapulator.logstream = open(args.logfile, 'w')
    microhapulator.teelog = args.tee
    return args
