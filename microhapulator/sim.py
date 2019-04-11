#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


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
