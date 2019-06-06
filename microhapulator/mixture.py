#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import microhapulator
import numpy
from string import ascii_letters, digits
import sys


def calc_n_reads_from_proportions(n, totalreads, prop):
    if prop is None:
        prop = [1.0 / n for _ in range(n)]
    else:
        if len(prop) != n:
            raise ValueError('mismatch between contributor number and proportions')
    normprop = [x / sum(prop) for x in prop]
    return [int(totalreads * x) for x in normprop]


def mixture(individuals, panel, relaxed=False, hapseeds=None, seqseeds=None,
            seqthreads=1, totalreads=1000000, proportions=None, gtfile=None):
    n = len(individuals)
    if hapseeds is None:
        hapseeds = [numpy.random.randint(1, 2**32 - 1) for _ in range(n)]
    if seqseeds is None:
        seqseeds = [numpy.random.randint(1, 2**32 - 1) for _ in range(n)]
    if n != len(hapseeds) or n != len(seqseeds):
        msg = 'number of individuals must match number of "--hap-seeds" and "--seq-seeds"'
        raise ValueError(msg)
    if gtfile:
        genotypes = list()
        for popids, hapseed in zip(individuals, hapseeds):
            genotype = microhapulator.sim.simulate_genotype(
                popids, panel, hapseed=hapseed, relaxed=relaxed,
            )
            genotypes.append(genotype)
        merged_genotype = microhapulator.genotype.SimulatedGenotype.merge(genotypes)
        merged_genotype.dump(gtfile)
    numreads = calc_n_reads_from_proportions(n, totalreads, proportions)
    if 0 in numreads:
        raise ValueError('specified proportions result in 0 reads for 1 or more individuals')
    readsignature = microhapulator.sim.new_signature()
    reads_sequenced = 0
    for indiv, hapseed, seqseed, nreads in zip(individuals, hapseeds, seqseeds, numreads):
        message = 'Individual population={pop} numreads={n}'.format(pop=','.join(indiv), n=nreads)
        microhapulator.plog('[MicroHapulator::mixture]', message)
        simulator = microhapulator.sim.sim(
            indiv, panel, relaxed=relaxed, hapseed=hapseed, seqseed=seqseed,
            seqthreads=seqthreads, numreads=nreads, readsignature=readsignature,
            readindex=reads_sequenced,
        )
        for data in simulator:
            yield data
        reads_sequenced = data[0]


def main(args=None):
    if args is None:  # pragma: no cover
        args = get_parser().parse_args()
    simulator = mixture(
        args.indiv, args.panel, relaxed=args.relaxed, hapseeds=args.hap_seeds,
        seqseeds=args.seq_seeds, seqthreads=args.seq_threads, totalreads=args.num_reads,
        proportions=args.proportions, gtfile=args.genotype,
    )
    with microhapulator.open(args.out, 'w') as fh:
        for n, defline, sequence, qualities in simulator:
            print(defline, sequence, '+\n', qualities, sep='', end='', file=fh)
