#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


import argparse
import numpy
import microhapulator
import microhapdb
import sys


def get_parser():
    cli = argparse.ArgumentParser()
    cli.add_argument('--panel', nargs='+', metavar='ID', help='list of '
                     'MicroHapDB locus IDs to simulate; by default, a panel '
                     'of 22 ALFRED microhaplotype loci is used')
    cli.add_argument('-r', '--relaxed', action='store_true', help='randomly '
                     'draw an allele (from a uniform distribution) for locus '
                     'where no allele frequency data is available for the '
                     'requested population')
    cli.add_argument('-s', '--seed', type=int, default=None, metavar='SEED',
                     help='seed for random number generator')
    cli.add_argument('refr', help='reference genome file')
    cli.add_argument('popid', nargs='+', help='population ID(s)')
    return cli


def validate_populations(popids):
    if len(popids) not in (1, 2):
        message = 'please provide only 1 or 2 population IDs'
        raise ValueError(message)
    haplopops = list()
    invalidids = set()
    for popid in popids:
        try:
            pop = microhapdb.id_xref(popid)
            haplopops.append(list(pop.ID)[0])
        except StopIteration:
            invalidids.add(popid)
    if len(invalidids) > 0:
        message = 'invalid population ID(s) "{}"'.format(','.join(invalidids))
        raise ValueError(message)
    if len(haplopops) == 1:
        haplopops = haplopops * 2
    return haplopops


def validate_loci(popids, panel=None, relaxed=False):
    if panel is None:
        loci = microhapdb.loci.query('Source == "ALFRED"').\
            sort_values('AvgAe', ascending=False).\
            drop_duplicates('Chrom')
    else:
        loci = microhapdb.loci[microhapdb.loci.ID.isin(panel)]
    if relaxed:
        return list(loci.ID)

    loci_to_keep = list()
    for locusid in list(loci.ID):
        keep = True
        for popid in popids:
            f = microhapdb.frequencies
            allelefreqs = f[(f.Population == popid) & (f.Locus == locusid)]
            if len(allelefreqs) == 0:
                keep = False
                message = 'no allele frequencies available'
                message += ' for population "{pop}"'.format(pop=popid)
                message += ' at locus "{loc}"'.format(loc=locusid)
                message += '; excluding from simulation'
                print('WARNING:', message, file=sys.stderr)
        if keep:
            loci_to_keep.append(locusid)
    return loci_to_keep


def sample_panel(popids, loci):
    for haplotype, popid in enumerate(popids):
        for locusid in loci:
            f = microhapdb.frequencies
            allelefreqs = f[(f.Population == popid) & (f.Locus == locusid)]
            if len(allelefreqs) == 0:
                message = 'no allele frequencies available'
                message += ' for population "{pop}"'.format(pop=popid)
                message += ' at locus "{loc}"'.format(loc=locusid)
                message += '; in "relaxed" mode, drawing an allele uniformly'
                print('WARNING:', message, file=sys.stderr)
                alleles = list(f[f.Locus == locusid].Allele.unique())
                sampled_allele = numpy.random.choice(alleles)
            else:
                alleles = list(allelefreqs.Allele)
                freqs = list(allelefreqs.Frequency)
                freqs = [x / sum(freqs) for x in freqs]
                sampled_allele = numpy.random.choice(alleles, p=freqs)
            yield haplotype, locusid, sampled_allele


def main(args=None):
    if args is None:  # pragma: no cover
        args = get_parser().parse_args()

    haplopops = validate_populations(args.popid)
    loci = validate_loci(haplopops, panel=args.panel, relaxed=args.relaxed)
    genotype = microhapulator.Genotype()
    if args.seed:
        numpy.random.seed(args.seed)
    for haplotype, locus, allele in sample_panel(haplopops, loci):
        genotype.add(haplotype, locus, allele)
    print(genotype)
    print(len(args.popid), 'populations and', len(loci), 'microhap loci')
