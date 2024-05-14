# -------------------------------------------------------------------------------------------------
# Copyright (c) 2023, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from Bio import SeqIO
from collections import defaultdict
from microhapdb.nomenclature import Identifier
from microhapulator import open as mhopen
import pandas as pd
from warnings import warn


class MicrohapIndex:
    def __init__(self, loci=None, seqmap=None):
        self.sequences = dict()
        if seqmap is not None:
            for seqid, sequence in seqmap.items():
                self.sequences[seqid] = sequence
        self.loci = dict()
        self.markers = dict()
        if loci is not None:
            for locus in loci:
                if locus.id in self.loci:
                    raise ValueError(f"duplicate locus ID {locus.id}")
                self.loci[locus.id] = locus
                if locus.id in self.sequences:
                    locus.sequence = self.sequences[locus.id]
                for marker in locus:
                    if marker.id in self.markers:
                        raise ValueError(f"duplicate marker ID {marker.id}")
                    self.markers[marker.id] = marker

    @classmethod
    def from_files(cls, csv_path, fasta_path=None):
        definitions = cls.parse_definitions_from_csv(csv_path)
        sequences = cls.parse_sequences_from_fasta(fasta_path) if fasta_path else None
        return cls(loci=definitions, seqmap=sequences)

    @staticmethod
    def parse_definitions_from_csv(path):
        table = pd.read_csv(path, sep=None, engine="python")
        MicrohapIndex.check_definition_columns(table)
        loci = defaultdict(Locus)
        for markerid, rowgroup in table.groupby("Marker"):
            offsets = rowgroup.Offset
            chrom = MicrohapIndex.marker_chromosome(rowgroup)
            chrom_offsets = MicrohapIndex.chromosome_offsets(rowgroup)
            marker = MarkerDefinition(markerid, offsets, chrom_offsets, chrom)
            loci[marker.locus].add(marker)
        return loci.values()

    @staticmethod
    def check_definition_columns(table):
        missing = set(["Marker", "Offset", "Chrom", "OffsetHg38"]) - set(table.columns)
        if len(missing) > 0:
            missing = ", ".join(sorted(missing))
            message = f"column(s) missing from marker definition file: {missing}"
            raise ValueError(message)

    @staticmethod
    def marker_chromosome(rowgroup):
        if len(rowgroup.Chrom.unique()) > 1:
            chroms = ", ".join(sorted(rowgroup.Chrom.unique()))
            markerid = rowgroup.Marker.iloc[0]
            message = f"ambiguous chromosome definition for {markerid}: {chroms}"
            raise ValueError(message)
        return rowgroup.Chrom.iloc[0]

    @staticmethod
    def chromosome_offsets(rowgroup):
        if rowgroup.OffsetHg38.isnull.any():
            markerid = rowgroup.Marker.iloc[0]
            message = f"incomplete marker definition for {markerid}, includes null values"
            raise ValueError(message)
        return rowgroup.OffsetHg38

    @staticmethod
    def parse_sequences_from_fasta(path):
        with mhopen(path, "r") as fh:
            sequences = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
            sequences = {seqid: record.seq for seqid, record in sequences.items()}
        return sequences

    def __getitem__(self, key):
        if key in self.loci:
            return self.loci[key]
        elif key in self.markers:
            return self.markers[key]
        else:
            raise KeyError(key)

    def sequence(self, key):
        for locus in self.loci.values():
            if locus.id == key:
                return self.sequences[locus.id]
            for marker in locus:
                if marker.id == key:
                    return self.sequences[locus.id]
        raise KeyError(key)

    def __iter__(self):
        for locusid, locus in sorted(self.loci.items()):
            for marker in locus:
                yield locus, marker

    def validate(self, refrids=None, symmetric=False):
        seqnames = set(refrids) if refrids is not None else set(self.sequences.keys())
        locusnames = set(self.loci.keys())
        if len(locusnames) == 0:
            warn("no marker definitions in microhap index")
            return
        missing = locusnames - seqnames
        if len(missing) > 0:
            missing = ", ".join(sorted(missing))
            message = f"reference sequences missing for the following loci: {missing}"
            raise ValueError(message)
        if symmetric:
            extra = seqnames - locusnames
            if len(extra) > 0:
                extra = ", ".join(sorted(extra))
                message = f"reference sequences with no marker definition: {extra}"
                raise ValueError(message)


class Locus:
    def __init__(self, ident=None, markers=None, sequence=None):
        self.identifier = None
        if ident:
            self.identifier = Identifier(ident)
        self.marker_definitions = list()
        if markers is not None:
            for marker in markers:
                self.marker_definitions.add(marker)
        self.sequence = sequence

    def __iter__(self):
        yield from self.marker_definitions

    def add(self, marker):
        if self.identifier is None:
            self.identifier = Identifier(marker.locus)
        else:
            locusid = str(self.identifier)
            if locusid != marker.locus:
                warn(f"locus ID mismatch: {locusid} vs {marker.locus}")
        self.marker_definitions.append(marker)

    @property
    def id(self):
        if self.identifier.valid:
            if self.identifier.suffix is None:
                return str(self.identifier)
            else:
                raise ValueError(f"locus ID cannot include a suffix: {self.identifier._raw}")
        else:
            return self.identifier._raw

    @property
    def chrom(self):
        chroms = [m.chrom for m in self]
        if len(set(chroms)) > 1:
            raise ValueError(f"ambiguous chromosome annotation for locus {self.id}")
        return chroms[0]

    def span(self, chrom=False):
        start, end = float("Inf"), -1
        for marker in self.marker_definitions:
            mstart, mend = marker.span(chrom=chrom)
            start = min(mstart, start)
            end = max(mend, end)
        return start, end

    def offsets(self, chrom=False):
        combined = set()
        for marker in self.marker_definitions:
            if chrom:
                off = marker.offsets_chrom
            else:
                off = marker.offsets_locus
            combined |= set(off)
        return combined


class MarkerDefinition:
    def __init__(self, ident, offsets_locus, offsets_chrom, chrom):
        self.identifier = Identifier(ident)
        self.offsets_locus = sorted(map(int, offsets_locus))
        self.offsets_chrom = sorted(map(int, offsets_chrom))
        self._chrom = chrom

    @property
    def id(self):
        if self.identifier.valid:
            return str(self.identifier)
        else:
            return self.identifier._raw

    @property
    def locus(self):
        if self.identifier.valid:
            return self.identifier.locus
        else:
            return self.identifier._raw

    @property
    def chrom(self):
        if self.identifier.valid:
            return self.identifier.chrom
        else:
            return self._chrom

    def span(self, chrom=False):
        if chrom:
            return min(self.offsets_chrom), max(self.offsets_chrom) + 1
        else:
            return min(self.offsets_locus), max(self.offsets_locus) + 1
