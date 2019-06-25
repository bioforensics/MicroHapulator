#!/usr/bin/env python

import argparse
from gzip import open as gzopen
import microhapulator
import multiprocessing
import pandas


cli = argparse.ArgumentParser()
cli.add_argument('threads', type=int, default=1)
cli.add_argument('data')
args = cli.parse_args()


def simulate(data):
    index, row = data
    genotypes = [
        microhapulator.genotype.SimulatedGenotype(fromfile=contrib + '-simulated-genotype.json')
        for contrib in row['Contributors'].split(',')
    ]
    seeds = list(map(int, row['SeqSeeds'].split(',')))
    proportions = list(map(float, row['Proportions'].split(',')))
    assert len(genotypes) == len(seeds)
    assert len(genotypes) == len(proportions)
    sequencer = microhapulator.seq.seq(
        genotypes, seeds=seeds, threads=1, totalreads=250000,
        proportions=proportions, sig=row['ID']
    )
    with gzopen(row['ID'] + '-reads.fastq.gz', 'wt') as fh:
        for n, defline, sequence, qualities in sequencer:
            print(defline, sequence, '+\n', qualities, sep='', end='', file=fh)
    print('Simulated', n, 'reads for sample', row['ID'])


column_types = {'Contributors': str, 'SeqSeeds': str, 'Proportions': str}
data = pandas.read_csv(args.data, sep='\t', dtype=column_types)
with multiprocessing.Pool(args.threads) as pool:
    pool.map(simulate, data.iterrows())
