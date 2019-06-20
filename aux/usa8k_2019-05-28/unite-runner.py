#!/usr/bin/env python

import argparse
import microhapulator
import multiprocessing
import numpy
import pandas


cli = argparse.ArgumentParser()
cli.add_argument('threads', type=int, default=1)
cli.add_argument('data')
args = cli.parse_args()


def unite(data):
    index, row = data
    mom = microhapulator.genotype.SimulatedGenotype(fromfile=row['MaternalID'] + '-genotype.json')
    dad = microhapulator.genotype.SimulatedGenotype(fromfile=row['PaternalID'] + '-genotype.json')
    numpy.random.seed(row['Seed'])
    kid = microhapulator.genotype.Genotype.unite(mom, dad)
    filename = row['ID'] + '-genotype.json'
    kid.dump(filename)


data = pandas.read_csv(args.data, sep='\t')
with multiprocessing.Pool(args.threads) as pool:
    pool.map(unite, data.iterrows())
