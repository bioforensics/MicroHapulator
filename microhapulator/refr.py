#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


import microhapulator
from microhapulator.panel import LocusContext, default_panel
import microhapdb
from pyfaidx import Fasta as Fastaidx


def get_seqs(locusids, seqindex, delta=30, minlength=350):
    locusids = microhapdb.standardize_ids(locusids)
    loci = microhapdb.loci[microhapdb.loci.ID.isin(locusids)]
    for i, rowdata in loci.iterrows():
        context = LocusContext(rowdata, mindelta=delta, minlen=minlength)
        yield context.defline(), context.sequence(seqindex)


def main(args=None):
    if args is None:  # pragma: no cover
        args = get_parser().parse_args()
    locusids = args.panel if args.panel else default_panel()
    seqindex = Fastaidx(args.refrfasta)
    with microhapulator.open(args.out, 'w') as fh:
        seqiter = get_seqs(locusids, seqindex, delta=args.delta, minlength=args.min_length)
        for defline, sequence in seqiter:
            print('>', defline, '\n', sequence, sep='', file=fh)
