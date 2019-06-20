#!/usr/bin/env python

import argparse
import microhapulator
import multiprocessing
import pandas


cli = argparse.ArgumentParser()
cli.add_argument('threads', type=int, default=1)
cli.add_argument('data')
args = cli.parse_args()

usapanel = microhapulator.panel.panel_usa()


def simulate(data):
    index, row = data
    populations = [
        row['MaternalHaploPop'],
        row['PaternalHaploPop'],
    ]
    genotype = microhapulator.sim.sim(populations, usapanel, seed=row['HaploSeed'])
    filename = row['ID'] + '-genotype.json'
    genotype.dump(filename)


data = pandas.read_csv(args.data, sep='\t')
with multiprocessing.Pool(args.threads) as pool:
    pool.map(simulate, data.iterrows())
