#!/usr/bin/env python3

import argparse
import numpy
import pandas
import sys

cli = argparse.ArgumentParser()
cli.add_argument('-r', '--num-reads', type=int, default=100000)
cli.add_argument('-s', '--seed', type=int)
cli.add_argument('panel')
cli.add_argument('demo')
cli.add_argument('n', type=int)
args = cli.parse_args()

panel = pandas.read_csv(args.panel, sep='\t')
demo = pandas.read_csv(args.demo, sep='\t')
demo.Weight /= sum(demo.Weight)
if not args.seed:
    args.seed = numpy.random.randint(2**32-1)
print('Using random seed', args.seed, file=sys.stderr)
numpy.random.seed(args.seed)

for _ in range(args.n):
    pop1 = numpy.random.choice(demo.Population, p=demo.Weight)
    if numpy.random.random() < 0.1:
        otherpops = demo[demo.Population != pop1].copy()
        otherpops.Weight /= sum(otherpops.Weight)
        pop2 = numpy.random.choice(otherpops.Population, p=otherpops.Weight)
    else:
        pop2 = pop1
    haploseed = numpy.random.randint(2**32-1)
    seqseed = numpy.random.randint(2**32-1)
    panelstr = ' '.join(panel.ID)
    prefix = '{p1:s}-{p2:s}-{hapseed:d}-{seqseed:d}'.format(
        p1=pop1, p2=pop2, hapseed=haploseed, seqseed=seqseed
    )
    cmd = (
        'mhpl8r sim --panel {panel:s} --hap-seed {hapseed:d} --num-reads {nr:d}'
        ' --seq-seed {seqseed:d} --genotype {pfx:s}.bed --out {pfx:s}.fastq'
        ' hg38.fasta {p1:s}'.format(
            panel=panelstr, hapseed=haploseed, seqseed=seqseed, pfx=prefix,
            p1=pop1, nr=args.num_reads,
        )
    )
    if pop1 != pop2:
        cmd += ' ' + pop2
    print(cmd)
