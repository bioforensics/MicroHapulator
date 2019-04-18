#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

# Core library imports
from collections import defaultdict
from io import StringIO
from os import fsync
from shutil import rmtree
from string import ascii_letters, digits
from subprocess import check_call
from tempfile import NamedTemporaryFile, mkdtemp

# Third-party library imports
from happer.mutate import mutate
import microhapdb
from pyfaidx import Fasta as Fastaidx
from numpy.random import seed, choice

# Internal imports
import microhapulator
from microhapulator.genotype import SimulatedGenotype
from microhapulator.panel import LocusContext, panel_loci, sample_panel
from microhapulator.panel import validate_populations, exclude_loci_missing_data


def optional_outfile(outfile):
    if outfile:
        return open(outfile, 'w')
    else:
        return NamedTemporaryFile(mode='wt', suffix='.fasta')


def main(args=None):
    if args is None:  # pragma: no cover
        args = get_parser().parse_args()

    haplopops = validate_populations(args.popid)
    loci = panel_loci(args.panel)
    if not args.relaxed:
        loci = exclude_loci_missing_data(loci, haplopops)
    if loci in (None, list()):
        raise ValueError('invalid panel: {}'.format(args.panel))
    genotype = SimulatedGenotype()
    if args.hap_seed:
        seed(args.hap_seed)
    for haplotype, locus, allele in sample_panel(haplopops, loci):
        genotype.add(haplotype, locus, allele)
    if args.genotype:
        with open(args.genotype, 'w') as fh:
            print(genotype, file=fh)

    message = 'simulated microhaplotype variation at {loc:d} loci'.format(loc=len(loci))
    microhapulator.plog('[MicroHapulator::sim]', message)

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
            microhapulator.logstream.flush()
            try:
                fsync(microhapulator.logstream.fileno())
            except OSError:  # pragma: no cover
                pass
            check_call(isscmd, stderr=microhapulator.logstream)
            with open(fqdir + '/seq_R1.fastq', 'r') as infh, open(args.out, 'w') as outfh:
                signature = ''.join([choice(list(ascii_letters + digits)) for _ in range(7)])
                nreads = 0
                for line in infh:
                    if line.startswith('@MHDBL'):
                        nreads += 1
                        prefix = '@{sig:s}_read{n:d} MHDBL'.format(sig=signature, n=nreads)
                        line = line.replace('@MHDBL', prefix, 1)
                    print(line, end='', file=outfh)
        finally:
            rmtree(fqdir)
