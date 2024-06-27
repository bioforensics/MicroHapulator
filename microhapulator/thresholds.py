# -------------------------------------------------------------------------------------------------
# Copyright (c) 2024, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

import pandas as pd


class ThresholdIndex:
    def __init__(self, static=5, dynamic=0.035, ambiguous=0.2, min_read_length=50, workdir="."):
        self.global_static = static
        self.global_dynamic = dynamic
        self.marker_static = {}
        self.marker_dynamic = {}
        self.ambiguous = ambiguous
        self.min_read_length = min_read_length

    def set(self, marker, static=None, dynamic=None):
        if static:
            self.marker_static[marker] = static
        if dynamic:
            self.marker_dynamic[marker] = dynamic

    def get(self, marker):
        static = self.global_static
        dynamic = self.global_dynamic
        if marker in self.marker_static:
            static = self.marker_static[marker]
        if marker in self.marker_dynamic:
            dynamic = self.marker_dynamic[marker]
        return static, dynamic

    @classmethod
    def load(
        cls,
        configfile=None,
        global_static=5,
        global_dynamic=0.02,
        ambiguous=0.2,
        min_read_length=50,
    ):
        index = cls(
            static=global_static,
            dynamic=global_dynamic,
            ambiguous=ambiguous,
            min_read_length=min_read_length,
        )
        if configfile:
            config = pd.read_csv(configfile, sep=None, engine="python")
            config = config.fillna(0).astype({"Static": "Int64"})
            missing = set(["Marker", "Static", "Dynamic"]) - set(config.columns)
            if len(missing) > 0:
                missingstr = ",".join(sorted(missing))
                raise ValueError(f"filter config file missing column(s): {missingstr}")
            if len(config.Marker) != len(config.Marker.unique()):
                raise ValueError("filter config file contains duplicate entries for some markers")
            config = config.set_index("Marker", drop=True)
            for marker, row in config.iterrows():
                index.set(marker, static=row.Static, dynamic=row.Dynamic)
        return index
