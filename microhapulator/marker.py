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

from collections import defaultdict
from microhapdb.nomenclature import Identifier
import pandas as pd


class DefinitionIndex:
    def __init__(self):
        self.offsets = defaultdict(lambda: defaultdict(list))
        self.chrom_offsets = defaultdict(lambda: defaultdict(list))
        self.marker_chroms = dict()

    @classmethod
    def from_csv(cls, path, strict=False):
        table = pd.read_csv(path, sep=None, engine="python")
        missing = set(["Marker", "Offset"]) - set(table.columns)
        if len(missing) > 0:
            missing = ", ".join(sorted(missing))
            message = f"column(s) missing from marker definition file: {missing}"
            raise ValueError(message)
        index = cls()
        if "OffsetHg38" not in table.columns:
            index.chrom_offsets = None
        for i, entry in table.iterrows():
            chrom = None
            if not strict and hasattr(entry, "Chrom") and not pd.isna(entry.Chrom):
                chrom = entry.Chrom
            if index.chrom_offsets is not None and pd.isna(entry.OffsetHg38):
                raise ValueError("empty OffsetHg38 value in marker definition file")
            o38 = entry.OffsetHg38 if hasattr(entry, "OffsetHg38") else None
            index.add(entry.Marker, entry.Offset, offset38=o38, chrom=chrom, strict=strict)
        return index

    def add(self, markerid, offset, offset38=None, chrom=None, strict=False):
        if strict:
            ident = Identifier(markerid)
            if not ident.valid:
                msg = f"invalid identifier {entry.Marker}: {ident.errors}; {ident.warnings}"
                raise ValueError(msg)
            self.offsets[ident.locus][str(ident)].append(offset)
            if offset38:
                self.chrom_offsets[ident.locus][str(ident)].append(offset38)
            self.marker_chroms[ident.locus] = ident.chrom
            self.marker_chroms[markerid] = ident.chrom
        else:
            self.offsets[markerid][markerid].append(offset)
            if offset38:
                self.chrom_offsets[markerid][markerid].append(offset38)
            if chrom:
                if markerid in self.marker_chroms and self.marker_chroms[markerid] != chrom:
                    raise ValueError(f"chrom mismatch: {chrom} vs {self.marker_chroms[markerid]}")
                self.marker_chroms[markerid] = chrom

    def __getitem__(self, key):
        return self.offsets[key]

    def __iter__(self):
        for locusid, markers in self.offsets.items():
            for markerid, offsets in markers.items():
                yield markerid, locusid, set(offsets)

    @property
    def loci(self):
        return sorted(self.offsets.keys())

    @property
    def has_chrom_offsets(self):
        return self.chrom_offsets is not None and len(self.chrom_offsets) > 0

    def offset_span(self, locus, chrom=False):
        positions = self.chrom_offsets if chrom else self.offsets
        start, end = float("Inf"), -1
        for marker, offsets in positions[locus].items():
            start = min(min(offsets), start)
            end = max(max(offsets) + 1, end)
        return start, end

    def all_offsets(self, locus, chrom=False):
        positions = self.chrom_offsets if chrom else self.offsets
        combined = set()
        for marker, offsets in positions[locus].items():
            combined |= set(offsets)
        return combined
