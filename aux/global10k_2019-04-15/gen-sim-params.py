#!/usr/bin/env python3

import argparse
import microhapdb
import numpy
import pandas
import sys

cli = argparse.ArgumentParser()
cli.add_argument('-r', '--num-reads', type=int, default=100000)
cli.add_argument('-s', '--seed', type=int)
cli.add_argument('n', type=int)
args = cli.parse_args()

panel = microhapdb.loci.query('Source == "ALFRED"').sort_values('AvgAe', ascending=False).head(100)
panelstr = ' '.join(panel.ID)

if not args.seed:
    args.seed = numpy.random.randint(2**32-1)
print('Using random seed', args.seed, file=sys.stderr)
numpy.random.seed(args.seed)

populations = microhapdb.populations.query('Source == "ALFRED"')
for popid in populations.ID:
    for _ in range(args.n):
        haploseed = numpy.random.randint(2**32-1)
        seqseed = numpy.random.randint(2**32-1)
        prefix = '{pop:s}-{pop:s}-{hapseed:d}-{seqseed:d}'.format(
            pop=popid, hapseed=haploseed, seqseed=seqseed
        )
        cmd = (
            'mhpl8r sim --panel {panel:s} --hap-seed {hapseed:d} --num-reads {nr:d}'
            ' --seq-seed {seqseed:d} --genotype {pfx:s}.bed --out {pfx:s}.fastq'
            ' hg38.fasta {pop:s} {pop:s}'.format(
                panel=panelstr, hapseed=haploseed, seqseed=seqseed, pfx=prefix,
                pop=popid, nr=args.num_reads,
            )
        )
        print(cmd)
