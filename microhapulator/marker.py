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


class DefinitionIndex(defaultdict):
    @classmethod
    def from_csv(cls, path, strict=False):
        table = pd.read_csv(path, sep=None, engine="python")
        missing = set(["Marker", "Offset"]) - set(table.columns)
        if len(missing) > 0:
            missing = ", ".join(sorted(missing))
            message = f"column(s) missing from marker definition file: {missing}"
            raise ValueError(message)
        index = cls(lambda: defaultdict(list))
        for i, entry in table.iterrows():
            if strict:
                ident = Identifier(entry.Marker)
                if not ident.valid:
                    msg = f"invalid identifier {entry.Marker}: {ident.errors}; {ident.warnings}"
                    raise ValueError(msg)
                index[ident.locus][str(ident)].append(entry.Offset)
            else:
                index[entry.Marker][entry.Marker].append(entry.Offset)
        return index

    def __iter__(self):
        for locusid, markers in self.items():
            for markerid, offsets in markers.items():
                yield markerid, locusid, set(offsets)

    @property
    def loci(self):
        return sorted(self.keys())
