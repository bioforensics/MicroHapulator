#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from hashlib import sha1
import microhapulator
import pyfaidx
from shutil import copyfileobj
from subprocess import check_call
from urllib.request import urlopen


def download_refr():
    url = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
    outfilegz = microhapulator.package_file('hg38.fasta.gz')
    outfile = microhapulator.package_file('hg38.fasta')
    microhapulator.plog('[MicroHapulator::getrefr] Downloading GRCh38 reference')
    with urlopen(url) as instream:
        with open(outfilegz, 'wb') as outstream:
            copyfileobj(instream, outstream)
    microhapulator.plog('[MicroHapulator::getrefr] Decompressing reference')
    check_call(['gunzip', '-f', outfilegz])


def install_refr(infile):
    outfile = microhapulator.package_file('hg38.fasta')
    microhapulator.plog('[MicroHapulator::getrefr] Copying GRCh38 reference')
    with microhapulator.open(infile, 'r') as instream:
        with microhapulator.open(outfile, 'w') as outstream:
            copyfileobj(instream, outstream)


def compute_refr_hash():
    refrfile = microhapulator.package_file('hg38.fasta')
    sha = sha1()
    with microhapulator.open(refrfile, 'r') as fh:
        while True:
            block = fh.read(2**10)
            if not block:
                break
            sha.update(block.encode('utf-8'))
    return sha.hexdigest()


def getrefr(infile=None, debug=False):
    if infile is None:
        download_refr()
    else:
        install_refr(infile)
    hashes = [
        '09b3f6ab110124cb230d028c3689f803f2a95fba',
        'c44a85ab746fe98ce1945022c643a4c289d2f3ce',
    ]
    if compute_refr_hash() not in hashes:
        raise ValueError('checksum failure')
    microhapulator.plog('[MicroHapulator::getrefr] Indexing reference')
    refr = pyfaidx.Fasta(microhapulator.package_file('hg38.fasta'))
    if debug:
        microhapulator.plog(
            '[MicroHapulator::getrefr] Reference genome now in place at',
            microhapulator.package_file('hg38.fasta')
        )


def main(args):
    if args.path:
        print(microhapulator.package_file('hg38.fasta'))
        return
    getrefr(args.filepath, debug=args.debug)
