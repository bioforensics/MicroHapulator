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

from setuptools import setup
import versioneer

with open("README.md", "r") as infile:
    longdesc = infile.read()

setup(
    name="microhapulator",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Software package for simulating and analyzing microhaplotype sequence data",
    long_description=longdesc,
    long_description_content_type="text/markdown",
    url="https://github.com/bioforensics/microhapulator",
    author="Daniel Standage",
    author_email="daniel.standage@nbacc.dhs.gov",
    packages=[
        "microhapulator",
        "microhapulator.cli",
        "microhapulator.happer",
        "microhapulator.pipe",
        "microhapulator.tests",
    ],
    package_data={
        "microhapulator": [
            "microhapulator/Snakefile",
            "microhapulator/data/*",
            "microhapulator/tests/data/*",
            "microhapulator/tests/data/*/*",
        ]
    },
    include_package_data=True,
    install_requires=[
        "biopython",
        "insilicoseq>=1.5.4,<2.0",
        "jsonschema>=4.0",
        "matplotlib>=3.0",
        "microhapdb>=0.12",
        "multiqc>=1.14",
        "nbformat>=5.0,<5.6",
        "numpy>=1.19",
        "pandas>1.0",
        "pulp==2.3.1",
        "scipy>=1.7",
        "seaborn>=0.13.2",
        "snakemake>=8",
        "termgraph>=0.5",
        "tqdm>=4.0",
    ],
    entry_points={
        "console_scripts": [
            "mhpl8r = microhapulator.cli:main",
            "happer = microhapulator.happer.__main__:main",
        ]
    },
    classifiers=[
        "Environment :: Console",
        "Framework :: IPython",
        "Framework :: Jupyter",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Legal Industry",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    zip_safe=True,
)
