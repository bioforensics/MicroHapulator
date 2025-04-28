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


class InvalidGenomicCoordinateError(ValueError):
    pass


class Allele(object):
    """Class for handling alleles"""

    def __init__(self, seqid, start, end, alleleseq, refrseq=None):
        """Constructor for allele objects.

        Genomic intervals should use half-open 0-indexed interval notation:
        that is, the first nucleotide of a sequence has an index of 0, the
        first 10 nucleotides are designated by the interval [0, 10), and the
        next 10 nucleotides are designated by the interval [10, 20).

        For alleles associated with indels or structural variants, the length
        of the alleleseq can be different from the start and end coordinates of
        the allele, which always correspond to the reference sequence.
        """
        self.seqid = seqid
        self.start = start
        self.end = end
        self.seq = alleleseq
        self.refr = refrseq
        if self.end < self.start:
            message = "allele end cannot be before allele start; "
            message += "{:s}:{:d}-{:d}".format(seqid, start, end)
            raise InvalidGenomicCoordinateError(message)

    def transform(self, offset):
        """Apply an offset to this allele's coordinates.

        This can be helpful when converting between global (chromosome)
        coordinates and local (locus-relative) coordinates.
        """
        newstart = self.start + offset
        newend = self.end + offset
        if self.start + offset < 0:
            message = "invalid allele transformation; "
            message += "{:s}:{:d}-{:d}".format(self.seqid, newstart, newend)
            raise InvalidGenomicCoordinateError(message)
        self.start = newstart
        self.end = newend

    @property
    def slug(self):
        return "{:s}:{:d}-{:d}".format(self.seqid, self.start, self.end)

    @property
    def refrlength(self):
        return self.end - self.start

    def __len__(self):
        return len(self.seq)

    def __eq__(self, other):
        return (
            self.seqid == other.seqid
            and self.start == other.start
            and self.end == other.end
            and self.seq == other.seq
        )

    def __lt__(self, other):
        if self.seqid != other.seqid:
            return self.seqid < other.seqid
        if self.start != other.start:
            return self.start < other.start
        if self.end != other.end:
            return self.end < other.end
        if len(self) != len(other):
            return len(self) < len(other)
        return self.seq < other.seq


def parse_alleles(instream):
    """Read allele info from a BED-formatted input stream.

    The first 4 columns should be tab-delimited, and formatted as follows.

    chr1    7324977  7324978   A|T
    chr1    7325021  7325022   C|T
    chr1    7325106  7325107   G|C

    Any other columns, if present, are ignored.
    """
    for line in instream:
        if line.startswith("#") or line.strip() == "":
            continue
        seqid, start, end, alleles, *remainder = line.strip().split()
        start, end = int(start), int(end)
        yield tuple([Allele(seqid, start, end, a) for a in alleles.split("|")])
