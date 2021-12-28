# -------------------------------------------------------------------------------------------------
# Copyright (c) 2018, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

# Core libraries
import builtins
from contextlib import contextmanager
from gzip import open as gzopen
import os
from pkg_resources import resource_filename
import sys

# Third-party libraries
from Bio import SeqIO
import pandas as pd

# Internal modules
from microhapulator import profile

# Internal perations and command-line interface
from microhapulator.op import *
from microhapulator import __main__
from microhapulator import cli

# Version handling
from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions


logstream = None
teelog = False


def package_file(path):
    return resource_filename("microhapulator", path)


@contextmanager
def open(filename, mode):
    if mode not in ("r", "w"):
        raise ValueError('invalid mode "{}"'.format(mode))
    if filename in ["-", None]:
        filehandle = sys.stdin if mode == "r" else sys.stdout
        yield filehandle
    else:
        openfunc = builtins.open
        if filename.endswith(".gz"):
            openfunc = gzopen
            mode += "t"
        with openfunc(filename, mode) as filehandle:
            yield filehandle


def plog(*args, **kwargs):
    """Print logging output."""
    global logstream
    global teelog
    if logstream is not None:
        print(*args, **kwargs, file=logstream)
    if logstream is None or teelog:
        print(*args, **kwargs, file=sys.stderr)


def load_marker_frequencies(tsvfile):
    frequencies = pd.read_csv(tsvfile, sep="\t")
    missing = set(["Marker", "Haplotype", "Frequency"]) - set(frequencies.columns)
    if len(missing) > 0:
        message = "column(s) missing from marker frequency file: " + ", ".join(sorted(missing))
        raise ValueError(message)
    return frequencies


def load_marker_definitions(tsvfile):
    markers = pd.read_csv(tsvfile, sep="\t")
    missing = set(["Marker", "Offset"]) - set(markers.columns)
    if len(missing) > 0:
        message = "column(s) missing from marker definition file: " + ", ".join(sorted(missing))
        raise ValueError(message)
    return markers


def load_marker_reference_sequences(fastafile):
    with open(fastafile, "r") as fh:
        sequences = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
        sequences = {seqid: record.seq for seqid, record in sequences.items()}
    return sequences


def cross_check_marker_ids(set1, set2, label1, label2):
    set1 = set(set1)
    set2 = set(set2)
    if set1 == set2:
        return
    uniq1 = set1 - set2
    uniq2 = set2 - set1
    message = f"discrepancy between {label1} (set1) and {label2} (set2):"
    if uniq1:
        ustr1 = ", ".join(sorted(uniq1))
        message += f" marker IDs unique to set1={{{ustr1}}};"
    if uniq2:
        ustr2 = ", ".join(sorted(uniq2))
        message += f" marker IDs unique to set2={{{ustr2}}};"
    raise ValueError(message)
