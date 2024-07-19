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


class MarkerDetails:
    def __init__(self, locus, marker):
        self.locus = locus
        self.marker = marker

    @property
    def identifier(self):
        return self.marker.id

    @property
    def seq_length(self):
        return len(self.sequence)

    @property
    def sequence(self):
        return self.locus.sequence.strip().upper()

    @property
    def gc_content(self):
        seq = self.sequence
        return (seq.count("G") + seq.count("C")) / len(seq) * 100

    @property
    def marker_offsets(self):
        return ", ".join([str(o) for o in sorted(self.marker.offsets_locus)])

    @property
    def chrom_offsets(self):
        return ", ".join([str(o) for o in sorted(self.marker.offsets_chrom)])

    @property
    def chromosome(self):
        return self.marker.chrom

    @staticmethod
    def from_index(index):
        for locus, marker in index:
            yield MarkerDetails(locus, marker)
