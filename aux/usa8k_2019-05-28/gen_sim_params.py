#!/usr/bin/env python3

import argparse
from collections import defaultdict, namedtuple
import microhapulator
import numpy
import pandas
from string import ascii_letters, digits
import sys


def get_id(length=7):
    symbols = list(ascii_letters + digits)
    return ''.join([numpy.random.choice(symbols) for _ in range(length)])


Contributor = namedtuple('Contributor', 'label mompop dadpop haploseed')


class SimulatedSample(object):
    def __init__(self, even=False):
        self.label = get_id()
        self.contributors = list()
        self.seqseeds = list()
        self.evenprops = even
        self._prop = None

    def add_contributor(self, maternalpop, paternalpop, haploseed, seqseed):
        contrib = Contributor(get_id(), maternalpop, paternalpop, haploseed)
        self.contributors.append(contrib)
        self.seqseeds.append(seqseed)

    @property
    def ncontrib(self):
        return len(self.contributors)

    @property
    def popstring(self):
        assert self.ncontrib == 1
        pops = sorted(
            set(
                [
                    self.contributors[0].mompop,
                    self.contributors[0].dadpop,
                ]
            )
        )
        return ','.join(pops)

    @property
    def proportions(self):
        if self._prop is None:
            if self.ncontrib == 1:
                self._prop = [1.0]
            else:
                if self.evenprops:
                    self._prop = [1.0 / self.ncontrib] * self.ncontrib
                else:
                    samp_vals = sorted(
                        numpy.random.random(self.ncontrib - 1) * 0.9 + 0.05
                    )
                    total = 0.0
                    self._prop = list()
                    for val in samp_vals:
                        prop = val - total
                        self._prop.append(prop)
                        total += prop
                    self._prop.append(1.0 - total)
        return self._prop

    @property
    def simparams(self):
        return '\n'.join(['\t'.join(map(str, c)) for c in self.contributors])

    @property
    def seqparams(self):
        contribs = ','.join([c.label for c in self.contributors])
        seqseeds = ','.join(map(str, self.seqseeds))
        props = ','.join(map('{:.4f}'.format, self.proportions))
        return '\t'.join((self.label, contribs, seqseeds, props))


def haplotype_pair(popdata):
    pop1 = numpy.random.choice(popdata.Population, p=popdata.Weight)
    if numpy.random.random() < 0.1:
        otherpops = popdata[popdata.Population != pop1].copy()
        otherpops.Weight /= sum(otherpops.Weight)
        pop2 = numpy.random.choice(otherpops.Population, p=otherpops.Weight)
    else:
        pop2 = pop1
    return pop1, pop2


def simulate(popdata):
    pop1, pop2 = haplotype_pair(popdata)
    haploseed = numpy.random.randint(2**32-1)
    seqseed = numpy.random.randint(2**32-1)
    return pop1, pop2, haploseed, seqseed


def get_parser():
    cli = argparse.ArgumentParser()
    cli.add_argument(
        '--ncontrib', type=int, default=1, metavar='NCONTRIB',
        help='number of contributors in each sample; default is 1'
    )
    cli.add_argument(
        '--even', action='store_true', help='enforce even contribution when '
        'there are multiple contributors'
    )
    cli.add_argument(
        '--seed', type=int, help='seed for random number generator'
    )
    cli.add_argument(
        'n', type=int, help='number of samples to simulate'
    )
    cli.add_argument(
        'contribout', help='output file to store parameters for simulating '
        'sample contributors'
    )
    cli.add_argument(
        'seqout', help='output file to store parameters for simulating sample '
        'sequencing'
    )
    return cli


def gen_params(nsamples, ncontrib=1, even=False, seed=None):
    demofile = microhapulator.package_file('usa-demographics.tsv')
    demo = pandas.read_csv(demofile, sep='\t')
    demo.Weight /= sum(demo.Weight)
    if not seed:
        seed = numpy.random.randint(2**32-1)
    print('Using random seed', seed, file=sys.stderr)
    numpy.random.seed(seed)
    for x in range(nsamples):
        sample = SimulatedSample(even=even)
        for y in range(ncontrib):
            pop1, pop2, haploseed, seqseed = simulate(demo)
            sample.add_contributor(pop1, pop2, haploseed, seqseed)
        yield sample


def main(args):
    samples = defaultdict(list)
    simulator = gen_params(
        args.n, ncontrib=args.ncontrib, even=args.even, seed=args.seed
    )
    with open(args.contribout, 'w') as cout, open(args.seqout, 'w') as sout:
        columns = ['ID', 'Contributors', 'SeqSeeds', 'Proportions']
        print(*columns, sep='\t', file=sout)
        columns = ['ID', 'MaternalHaploPop', 'PaternalHaploPop', 'HaploSeed']
        print(*columns, sep='\t', file=cout)
        for sample in simulator:
            print(sample.seqparams, file=sout)
            print(sample.simparams, file=cout)


if __name__ == '__main__':
    main(get_parser().parse_args())
