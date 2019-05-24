#!/usr/bin/env python3

import argparse
import microhapulator
import numpy
import pandas
from string import ascii_letters, digits
import sys


def haplotype_pair(popdata):
    pop1 = numpy.random.choice(popdata.Population, p=popdata.Weight)
    if numpy.random.random() < 0.1:
        otherpops = popdata[popdata.Population != pop1].copy()
        otherpops.Weight /= sum(otherpops.Weight)
        pop2 = numpy.random.choice(otherpops.Population, p=otherpops.Weight)
    else:
        pop2 = pop1
    return pop1, pop2


def simulate(popdata, labelpfx, serial):
    pop1, pop2 = haplotype_pair(popdata)
    haploseed = numpy.random.randint(2**32-1)
    seqseed = numpy.random.randint(2**32-1)
    label = '{pfx:s}{serial:04d}'.format(pfx=labelpfx, serial=serial)
    return label, pop1, pop2, haploseed, seqseed


def proportion(ncontrib, even=False):
    if ncontrib == 1:
        return [1.0]
    if even:
        return [1.0 / ncontrib] * ncontrib
    samp_vals = sorted(numpy.random.random(ncontrib - 1) * 0.9 + 0.05)
    total = 0.0
    props = list()
    for val in samp_vals:
        prop = val - total
        props.append(prop)
        total += prop
    props.append(1.0 - total)
    return props


def get_parser():
    cli = argparse.ArgumentParser()
    cli.add_argument('--ncontrib', type=int, default=1, metavar='NCONTRIB', help='number of contributors in each sample; default is 1')
    cli.add_argument('--even', action='store_true', help='enforce even contribution when there are multiple contributors')
    cli.add_argument('--seed', type=int, help='seed for random number generator')
    cli.add_argument('label', help='label prefix for samples')
    cli.add_argument('n', type=int, help='number of samples to simulate')
    return cli


def gen_params(nsamples, label, ncontrib=1, even=False, seed=None):
    demofile = microhapulator.package_file('usa-demographics.tsv')
    demo = pandas.read_csv(demofile, sep='\t')
    demo.Weight /= sum(demo.Weight)
    if not seed:
        seed = numpy.random.randint(2**32-1)
    print('Using random seed', seed, file=sys.stderr)
    numpy.random.seed(seed)
    for x in range(nsamples):
        sampleid = 'usa8k-sample-' + ''.join([numpy.random.choice(list(ascii_letters + digits)) for _ in range(7)])
        for y, p in zip(range(ncontrib), proportion(ncontrib, even)):
            data = simulate(demo, label, x+1)
            yield [sampleid, *data, '{:.4f}'.format(p)]


def main(args):
    columns = ['ID', 'Label', 'MaternalHaploPop', 'PaternalHaploPop', 'HaploSeed', 'SeqSeed', 'Proportion']
    print(*columns, sep='\t')
    simulator = gen_params(args.n, args.label, ncontrib=args.ncontrib, even=args.even, seed=args.seed)
    for data in simulator:
        print(*data, sep='\t')


if __name__ == '__main__':
    main(get_parser().parse_args())
