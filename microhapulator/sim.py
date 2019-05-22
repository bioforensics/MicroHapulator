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


def new_signature():
    return ''.join([choice(list(ascii_letters + digits)) for _ in range(7)])


def simulate_genotype(popids, panel, hapseed=None, relaxed=False, outfile=None):
    haplopops = validate_populations(popids)
    loci = panel_loci(panel)
    if not relaxed:
        loci = exclude_loci_missing_data(loci, haplopops)
    if loci in (None, list()):
        raise ValueError('invalid panel: {}'.format(panel))
    genotype = SimulatedGenotype()
    if hapseed:
        seed(hapseed)
    for haplotype, locus, allele in sample_panel(haplopops, loci):
        genotype.add(haplotype, locus, allele)
    if outfile:
        with microhapulator.open(outfile, 'w') as fh:
            print(genotype, file=fh)
    message = 'simulated microhaplotype variation at {loc:d} loci'.format(loc=len(loci))
    microhapulator.plog('[MicroHapulator::sim]', message)
    return genotype


def sim(popids, panel, relaxed=False, hapseed=None, gtfile=None, hapfile=None,
        seqseed=None, seqthreads=2, numreads=500000, readsignature=None, readindex=0, debug=False):
    genotype = simulate_genotype(
        popids, panel, hapseed=hapseed, relaxed=relaxed, outfile=gtfile
    )
    seqindex = Fastaidx(microhapulator.package_file('hg38.fasta'))
    mutator = mutate(genotype.seqstream(seqindex), genotype.bedstream)
    with optional_outfile(hapfile) as fh:
        for defline, sequence in mutator:
            print('>', defline, '\n', sequence, sep='', file=fh)
        fh.flush()
        fsync(fh.fileno())
        fqdir = mkdtemp()
        try:
            isscmd = [
                'iss', 'generate', '--n_reads', str(numreads * 2), '--draft', fh.name,
                '--model', 'MiSeq', '--output', fqdir + '/seq'
            ]
            if seqseed:
                isscmd.extend(['--seed', str(seqseed)])
            if seqthreads:
                isscmd.extend(['--cpus', str(seqthreads)])
            try:
                microhapulator.logstream.flush()
                fsync(microhapulator.logstream.fileno())
            except (AttributeError, OSError):  # pragma: no cover
                pass
            if debug:
                check_call(isscmd, stderr=microhapulator.logstream)
            else:
                check_call(isscmd)
            with open(fqdir + '/seq_R1.fastq', 'r') as infh:
                if readsignature is None:
                    readsignature = new_signature()
                linebuffer = list()
                for line in infh:
                    if line.startswith('@MHDBL'):
                        readindex += 1
                        prefix = '@{sig:s}_read{n:d} MHDBL'.format(sig=readsignature, n=readindex)
                        line = line.replace('@MHDBL', prefix, 1)
                    linebuffer.append(line)
                    if len(linebuffer) == 4:
                        yield readindex, linebuffer[0], linebuffer[1], linebuffer[3]
                        linebuffer = list()
        finally:
            rmtree(fqdir)


def main(args=None):
    if args is None:  # pragma: no cover
        args = get_parser().parse_args()
    simulator = sim(
        args.popid, args.panel, relaxed=args.relaxed, hapseed=args.hap_seed,
        gtfile=args.genotype, hapfile=args.haploseq, seqseed=args.seq_seed,
        seqthreads=args.seq_threads, numreads=args.num_reads,
    )
    with microhapulator.open(args.out, 'w') as fh:
        for n, defline, sequence, qualities in simulator:
            print(defline, sequence, '+\n', qualities, sep='', end='', file=fh)
