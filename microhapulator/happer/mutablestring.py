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


class MutableString(object):
    """Mutable string class

    A string class that supports editing of the string contents without
    creating a new copy of the string with each edit. The string is stored
    internally as a list of chararacters. Despite the overhead in converting
    strings to lists and back to strings when edits are done, the overall
    approach provides substantial performance improvements when, for example,
    applying mutations to a long DNA sequence.

    Borrows heavily from https://stackoverflow.com/a/10572792/459780.
    """

    def __init__(self, data):
        self.data = list(data)

    def __str__(self):
        return "".join(self.data)

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return str(self) == str(other)

    def __add__(self, chars):
        """Addition operator

        For consistency with other Python objects, the addition operator
        creates a new object rather than appending in place. However, this
        makes a copy of the data which goes against the performance this data
        structure optimizes for.
        """
        newdata = list(self.data) + list(str(chars))
        return MutableString("".join(newdata))

    def __iadd__(self, chars):
        self.data.extend(list(str(chars)))
        return self

    def __contains__(self, teststr):
        return teststr in str(self)

    def __setitem__(self, index, value):
        self.data[index] = value

    def __getitem__(self, index):
        if type(index) == slice:
            return "".join(self.data[index])
        return self.data[index]

    def __delitem__(self, index):
        del self.data[index]

    def __len__(self):
        return len(self.data)
