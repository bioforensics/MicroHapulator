#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


from argparse import ArgumentParser
from happer.mutate import mutate
from numpy.random import seed
import microhapulator
from microhapulator.genotype import Genotype
from microhapulator.locus import default_panel, validate_loci, sample_panel
from microhapulator.population import validate_populations, check_loci_for_population
from microhapulator.population import exclude_loci_missing_data
import microhapdb
from os import fsync
from pyfaidx import Fasta as Fastaidx
from shutil import rmtree
from subprocess import check_call
from sys import stderr
from tempfile import NamedTemporaryFile, mkdtemp


def get_parser():
    cli = ArgumentParser()
    cli._positionals.title = 'Input configuration'
    cli._optionals.title = 'Miscellaneous'
    cli.add_argument(
        '-v', '--version', action='version',
        version='MicroHapulator version ' + microhapulator.__version__
    )

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


def optional_outfile(outfile):
    if outfile:
        return open(outfile, 'w')
    else:
        return NamedTemporaryFile(mode='wt', suffix='.fasta')


def main(args=None):
    if args is None:  # pragma: no cover
        args = get_parser().parse_args()

    haplopops = validate_populations(args.popid)
    loci = args.panel if args.panel else default_panel()
    loci = validate_loci(loci)
    if not args.relaxed:
        loci = exclude_loci_missing_data(loci, haplopops)
    genotype = Genotype()
    if args.hap_seed:
        seed(args.hap_seed)
    for haplotype, locus, allele in sample_panel(haplopops, loci):
        genotype.add(haplotype, locus, allele)
    if args.genotype:
        with open(args.genotype, 'w') as fh:
            print(genotype, file=fh)

    message = 'simulated microhaplotype variation at {loc:d} loci'.format(loc=len(loci))
    print('[MicroHapulator]', message, file=stderr)

    seqindex = Fastaidx(args.refr)
    mutator = mutate(genotype.seqstream(seqindex), genotype.bedstream)
    with optional_outfile(args.haploseq) as fh:
        for defline, sequence in mutator:
            print('>', defline, '\n', sequence, sep='', file=fh)
        fh.flush()
        fsync(fh.fileno())
        fqdir = mkdtemp()
        try:
            isscmd = [
                'iss', 'generate', '--n_reads', str(args.num_reads * 2), '--draft', fh.name,
                '--model', 'MiSeq', '--output', fqdir + '/seq'
            ]
            if args.seq_seed:
                isscmd.extend(['--seed', str(args.seq_seed)])
            if args.seq_threads:
                isscmd.extend(['--cpus', str(args.seq_threads)])
            check_call(isscmd)
            with open(fqdir + '/seq_R1.fastq', 'r') as infh, open(args.out, 'w') as outfh:
                nreads = 0
                for line in infh:
                    if line.startswith('@MHDBL'):
                        nreads += 1
                        prefix = '@read{:d} MHDBL'.format(nreads)
                        line = line.replace('@MHDBL', prefix, 1)
                    print(line, end='', file=outfh)
        finally:
            rmtree(fqdir)
