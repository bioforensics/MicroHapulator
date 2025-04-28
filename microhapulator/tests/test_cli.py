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
# -----------------------------------------------------------------------------

from microhapulator import open as mhopen
from microhapulator.cli import main
from microhapulator.tests import data_file
import pytest


def test_microhapulator_open():
    thefile = data_file("refr/three-loci-refr.fasta")
    with mhopen(thefile, "r") as filehandle:
        filecontents = filehandle.read()
        assert len(filecontents.strip().split("\n")) == 6
    with pytest.raises(ValueError, match=r'invalid mode "p"'):
        with mhopen(thefile, "p") as filehandle:
            pass


# These parameters and the following two test functions are mostly just smoke tests to make sure
# the code doesn't go up in flames when data is written to stdout.
# @pytest.fixture(params=[
#     ("diff", [data_file("prof/euramer-sim-gt.json"), data_file("prof/euramer-inf-gt.json")]),
#     ("dist", [data_file("prof/euramer-sim-gt.json"), data_file("prof/euramer-inf-gt.json")]),
#     ("sim", [data_file("pashtun-sim/tiny-panel-multidef-freq.tsv")]),
#     ("unite", [data_file("prof/swedish-mom.json"), data_file("prof/swedish-dad.json")]),
# ])
# def stdout_regress_param(request):
#     return request.param
#
#
# def test_regression_stdout(stdout_regress_param, capsys, tmp_path):
#     command, infiles = stdout_regress_param
#     main([command, *infiles, "-o", outfile])
#     terminal = capsys.readouterr()
#     assert terminal.out == ""
#     main([command, *infiles])
#     assert terminal.out != ""


@pytest.mark.parametrize(
    "command, infiles",
    [
        ("diff", [data_file("prof/euramer-sim-gt.json"), data_file("prof/euramer-inf-gt.json")]),
        ("dist", [data_file("prof/euramer-sim-gt.json"), data_file("prof/euramer-inf-gt.json")]),
        ("sim", [data_file("freq/korea-5loc-freq.tsv")]),
        ("unite", [data_file("prof/swedish-mom.json"), data_file("prof/swedish-dad.json")]),
    ],
)
def test_regression_output_file(command, infiles, capsys, tmp_path):
    # This is mostly just a smoke test to make sure the code doesn't go up in flames when data is
    # written to stdout.
    outfile = str(tmp_path / "output.data")
    main([command, *infiles, "-o", outfile])
    terminal = capsys.readouterr()
    assert terminal.out == ""
    main([command, *infiles])
    terminal = capsys.readouterr()
    assert terminal.out != ""
