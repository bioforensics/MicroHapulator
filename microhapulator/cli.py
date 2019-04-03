#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


import argparse
import microhapulator
import microhapdb


def get_parser():
    cli = argparse.ArgumentParser()
    cli.add_argument('--panel', nargs='+', metavar='ID', help='list of '
                     'MicroHapDB locus IDs to simulate; by default, a panel '
                     'of 22 ALFRED microhaplotype loci is used')
    cli.add_argument('refr', help='reference genome file')
    cli.add_argument('popid', nargs='+', help='population ID(s)')
    return cli


def main(args=None):
    """MicroHapulator main method."""
    if args is None:  # pragma: no cover
        args = get_parser().parse_args()

    assert len(args.popid) in (1, 2)

    if args.panel:
        raise NotImplementedError('working on this feature')
    else:
        loci = microhapdb.loci.query('Source == "ALFRED"').\
            sort_values('AvgAe', ascending=False).\
            drop_duplicates('Chrom')

    print(len(args.popid), 'populations and', len(loci), 'microhap loci')
