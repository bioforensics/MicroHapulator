#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


import argparse
import happer
from happer.mutate import mutate
import numpy
import microhapulator
import microhapdb
import os
import pyfaidx
import shutil
import subprocess
import sys
import tempfile


def get_parser():
    cli = argparse.ArgumentParser()
    cli._positionals.title = 'Input configuration'
    cli._optionals.title = 'Miscellaneous'

    hapargs = cli.add_argument_group('Haplotype simulation')
    hapargs.add_argument(
        '--panel', nargs='+', metavar='ID', help='list of MicroHapDB locus '
        'IDs for which to simulate data; by default, a panel of 22 ALFRED '
        'microhaplotype loci is used'
    )
    hapargs.add_argument(
        '-r', '--relaxed', action='store_true', help='if a locus in the panel '
        'has no frequency data for a requested population, randomly draw an '
        'allele (from a uniform distribution) from all possible alleles; by '
        'default, these loci are exluded from simulation'
    )
    hapargs.add_argument(
        '--hap-seed', type=int, default=None, metavar='INT', help='random '
        'seed for simulating haplotypes'
    )

    seqargs = cli.add_argument_group('Targeted sequencing')
    seqargs.add_argument(
        '-n', '--num-reads', type=int, default=500000, metavar='N',
        help='number of reads to simulate; default is 500000'
    )
    seqargs.add_argument(
        '--seq-seed', type=int, default=None, metavar='INT', help='random '
        'seed for simulated sequencing'
    )
    seqargs.add_argument(
        '--seq-threads', type=int, default=None, metavar='INT', help='number '
        'of threads to use when simulating targeted amplicon sequencing'
    )
    outargs = cli.add_argument_group('Output configuration')
    outargs.add_argument(
        '-o', '--out', metavar='FILE', default='-', required=True,
        help='write simulated MiSeq reads in FASTQ format to FILE; use '
        '`/dev/stdout` to write reads to standard output'
    )
    outargs.add_argument(
        '--genotype', metavar='FILE', help='write simulated genotype data in '
        'BED format to FILE'
    )
    outargs.add_argument(
        '--haploseq', metavar='FILE', help='write simulated haplotype '
        'sequences in FASTA format to FILE'
    )

    cli.add_argument('refr', help='reference genome file')
    cli.add_argument('popid', nargs='+', help='population ID(s)')
    cli._action_groups[1], cli._action_groups[-1] = cli._action_groups[-1], cli._action_groups[1]
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


def optional_outfile(outfile):
    if outfile:
        return open(outfile, 'w')
    else:
        return tempfile.NamedTemporaryFile(suffix='.fasta')


def main(args=None):
    if args is None:  # pragma: no cover
        args = get_parser().parse_args()

    haplopops = validate_populations(args.popid)
    loci = validate_loci(haplopops, panel=args.panel, relaxed=args.relaxed)
    genotype = microhapulator.Genotype()
    if args.hap_seed:
        numpy.random.seed(args.hap_seed)
    for haplotype, locus, allele in sample_panel(haplopops, loci):
        genotype.add(haplotype, locus, allele)
    if args.genotype:
        with open(args.genotype, 'w') as fh:
            print(genotype, file=fh)

    message = 'simulated microhaplotype variation at {loc:d} loci'.format(loc=len(loci))
    print('[MicroHapulator]', message, file=sys.stderr)

    seqindex = pyfaidx.Fasta(args.refr)
    mutator = mutate(genotype.seqstream(seqindex), genotype.bedstream)
    with optional_outfile(args.haploseq) as fh:
        for defline, sequence in mutator:
            print('>', defline, '\n', sequence, sep='', file=fh)
        fh.flush()
        os.fsync(fh.fileno())
        fqdir = tempfile.mkdtemp()
        isscmd = [
            'iss', 'generate', '--n_reads', str(args.num_reads * 2), '--draft', fh.name,
            '--model', 'MiSeq', '--output', fqdir + '/seq'
        ]
        if args.seq_seed:
            isscmd.extend(['--seed', str(args.seq_seed)])
        if args.seq_threads:
            isscmd.extend(['--cpus', str(args.seq_threads)])
        subprocess.check_call(isscmd)
        with open(fqdir + '/seq_R1.fastq', 'r') as infh, open(args.out, 'w') as outfh:
            shutil.copyfileobj(infh, outfh)
        shutil.rmtree(fqdir)
