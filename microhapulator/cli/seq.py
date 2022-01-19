# -------------------------------------------------------------------------------------------------
# Copyright (c) 2019, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------


from argparse import Action, SUPPRESS
import microhapulator.api as mhapi
from microhapulator.parsers import load_marker_definitions, load_marker_reference_sequences
from microhapulator.profile import Profile
import sys


class OutfilesAction(Action):
    def __call__(self, parser, namespace, values, option_string=None):
        numargs = len(values)
        if numargs < 1 or numargs > 2:
            message = f'expected 1 or 2 output filenames, got {numargs}: {",".join(values)}'
            parser.error(message)
        filehandles = [open(filename, "w") for filename in values]
        setattr(namespace, self.dest, filehandles)


def subparser(subparsers):
    desc = "Simulate paired-end Illumina MiSeq sequencing of the given profile(s)"
    cli = subparsers.add_parser("seq", description=desc)
    cli.add_argument(
        "-o",
        "--out",
        nargs="+",
        default=[sys.stdout],
        action=OutfilesAction,
        help="write simulated paired-end MiSeq reads in FASTQ format to the specified file(s); if "
        "one filename is provided, reads are interleaved and written to the file; if two "
        "filenames are provided, reads are written to paired files; by default, reads are "
        "interleaved and written to the terminal (standard output)",
    )
    cli.add_argument(
        "-n",
        "--num-reads",
        type=int,
        default=500000,
        metavar="N",
        help="number of reads to simulate; default is 500000",
    )
    cli.add_argument(
        "-p",
        "--proportions",
        type=float,
        nargs="+",
        metavar="P",
        help="simulated mixture samples with multiple contributors at the "
        "specified proportions; by default even proportions are used",
    )
    cli.add_argument(
        "-s",
        "--seeds",
        nargs="+",
        type=int,
        default=None,
        metavar="INT",
        help="seeds for random number generator, 1 per profile",
    )
    cli.add_argument("--signature", default=None, help=SUPPRESS)
    cli.add_argument("--threads", type=int, default=1, metavar="INT", help=SUPPRESS)
    cli.add_argument("tsv", help="microhaplotype marker definitions in tabular (TSV) format")
    cli.add_argument("refrseqs", help="microhaplotype reference sequences in FASTA format")
    cli.add_argument(
        "profiles", nargs="+", help="one or more simple or complex profiles " "(JSON files)"
    )


def resolve_profiles(gtfiles):
    profiles = list()
    for gtfile in gtfiles:
        profile = Profile(fromfile=gtfile)
        for p in profile.unmix():
            profiles.append(p)
    return profiles


def main(args):
    if len(args.out) not in (1, 2):
        raise ValueError(f"expected 1 or 2 output files, found {len(args.out)}")
    fh1 = fh2 = args.out[0]
    if len(args.out) == 2:
        fh2 = args.out[1]
    profiles = resolve_profiles(args.profiles)
    markers = load_marker_definitions(args.tsv)
    refrseqs = load_marker_reference_sequences(args.refrseqs)
    sequencer = mhapi.seq(
        profiles,
        markers,
        refrseqs,
        seeds=args.seeds,
        threads=args.threads,
        totalreads=args.num_reads,
        proportions=args.proportions,
        sig=args.signature,
    )
    for n, read1, read2 in sequencer:
        print(read1.identifier, read1.sequence, "+\n", read1.quality, sep="", end="", file=fh1)
        print(read2.identifier, read2.sequence, "+\n", read2.quality, sep="", end="", file=fh2)
    for fh in args.out:
        if fh != sys.stdout:
            fh.close()
