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
from microhapulator.locus import LocusContext
from microhapulator.locus import default_panel, panel_alpha, validate_loci, sample_panel
from microhapulator.population import validate_populations, check_loci_for_population
from microhapulator.population import exclude_loci_missing_data


class Genotype(object):
    """Genotype represented by phased alleles at a number of specified loci.

    Alleles are stored in a dictionary, with microhap locus ID/name as the key
    and a list as the value. Each list contains 2 items, the haplotype/phase 0
    allele and the hap/phase 1 allele.

    >>> gt = Genotype()
    >>> gt.add(0, 'mh21KK-315', 'G,C,T')
    >>> gt.add(1, 'mh21KK-315', 'A,T,C')
    >>> gt.add(0, 'mh21KK-316', 'A,C,G,T')
    >>> gt.add(1, 'mh21KK-316', 'A,T,G,C')
    >>> print(gt)
    mh21KK-315 102     103     G|A
    mh21KK-315 207     208     C|T
    mh21KK-315 247     248     T|C
    mh21KK-316 108     109     A|A
    mh21KK-316 132     133     C|T
    mh21KK-316 179     180     G|G
    mh21KK-316 242     243     T|C
    """
    def __init__(self):
        self._data = defaultdict(lambda: [None] * 2)
        self._contexts = dict()

    def add(self, hapid, locusid, allele):
        assert hapid in (0, 1)
        self._data[locusid][hapid] = allele
        if locusid not in self._contexts:
            locus = microhapdb.id_xref(locusid).iloc[0]
            context = LocusContext(locus)
            self._contexts[locusid] = context

    def seqstream(self, seqindex, prechr=False):
        for locusid, context in sorted(self._contexts.items()):
            yield context.defline(), context.sequence(seqindex, prechr=prechr)

    @property
    def bedstream(self):
        for locusid in sorted(self._data):
            context = self._contexts[locusid]
            alleles_0 = self._data[locusid][0].split(',')
            alleles_1 = self._data[locusid][1].split(',')
            coords = microhapdb.allele_positions(locusid)
            for a0, a1, coord in zip(alleles_0, alleles_1, coords):
                localcoord = context.global_to_local(coord)
                allelestr = a0 + '|' + a1
                yield '\t'.join(
                    (locusid, str(localcoord), str(localcoord + 1), allelestr)
                )

    def __str__(self):
        out = StringIO()
        for line in self.bedstream:
            print(line, file=out)
        return out.getvalue()


def optional_outfile(outfile):
    if outfile:
        return open(outfile, 'w')
    else:
        return NamedTemporaryFile(mode='wt', suffix='.fasta')


def main(args=None):
    if args is None:  # pragma: no cover
        args = get_parser().parse_args()

    haplopops = validate_populations(args.popid)
    if args.panel in (None, ['default']):
        loci = default_panel()
    elif args.panel == ['alpha']:
        loci = panel_alpha()
    else:
        loci = args.panel
    loci = validate_loci(loci)
    if not args.relaxed:
        loci = exclude_loci_missing_data(loci, haplopops)
    if loci in (None, list()):
        raise ValueError('invalid panel: {}'.format(args.panel))
    genotype = Genotype()
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
            except OSError:
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
