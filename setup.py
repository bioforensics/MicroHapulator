#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from setuptools import setup
import versioneer

with open('README.md', 'r') as infile:
    longdesc = infile.read()

setup(
    name='microhapulator',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Simulator for microhaplotype sequences',
    long_description=longdesc,
    long_description_content_type='text/markdown',
    url='https://github.com/bioforensics/microhapdb',
    author='Daniel Standage',
    author_email='daniel.standage@nbacc.dhs.gov',
    packages=['microhapulator'],
    # package_data={
    #     'microhapulator': ['microhapulator/data/*']
    # },
    include_package_data=True,
    install_requires=['pyfaidx', 'microhapdb', 'happer'],
    entry_points={
        'console_scripts': ['mhpl8r = microhapulator.cli:main']
    },
    classifiers=[
        'Environment :: Console',
        'Framework :: IPython',
        'Framework :: Jupyter',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 3',
        'Topic : Scientific/Engineering :: Bio-Informatics',
    ],
    zip_safe=True,
)
