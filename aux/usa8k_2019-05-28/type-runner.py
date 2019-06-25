#!/usr/bin/env python

import argparse
import microhapulator
import multiprocessing
import pandas
import subprocess

cli = argparse.ArgumentParser()
cli.add_argument('threads', type=int, default=1)
cli.add_argument('data')
args = cli.parse_args()


def getpanel():
    microhapulator.refr.main(
        microhapulator.cli.get_parser().parse_args(
            ['refr', '--out', 'panel-usa.fasta', 'usa']
        )
    )
    subprocess.check_call(['bwa', 'index', 'panel-usa.fasta'])


def filenames(prefix):
    return (
        '{}-inferred-genotype.json'.format(prefix),
        '{}-reads.fastq.gz'.format(prefix),
        '{}-reads.bam'.format(prefix),
    )


def infergenotype(data):
    index, row = data
    gtfile, readsfile, bamfile = filenames(row['ID'])
    bwacmd = 'bwa mem panel-usa.fasta {fq:s} | samtools view -bS | samtools sort -o {bam:s} - && samtools index {bam:s}'.format(fq=readsfile, bam=bamfile)
    subprocess.check_call(bwacmd, shell=True)
    genotype = microhapulator.type.type(bamfile, 'panel-usa.fasta')
    genotype.dump(gtfile)


getpanel()
data = pandas.read_csv(args.data, sep='\t')
with multiprocessing.Pool(args.threads) as pool:
    pool.map(infergenotype, data.iterrows())
