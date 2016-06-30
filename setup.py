#!/usr/bin/env python

from setuptools import setup
from os import path

with open(path.join(path.abspath(path.dirname(__file__)), 'README.md')) as f:
    README = f.read()

setup(
    name='biseqt',
    version='0.0.1',
    description='Biseqt is a biological sequence analysis tool',
    long_description=README,
    url='https://github.com/amirkdv/biseqt',
    author='Amir Kadivar',
    author_email='amir@amirkdv.ca',
    license='BSD',
    packages=['biseqt'],
    # NOTE numpy has some weirdness with setuptools (cf. open issue
    # https://github.com/numpy/numpy/issues/2434). We need it because it's a
    # dependency of scipy but has to be installed either manually before our
    # setup.py or as part of setup_requires (cf.
    # https://github.com/numpy/numpy/issues/2434#issuecomment-65252402)
    # However, it is still much faster to install numpy separately via pip
    # because the following way wants to compile all the fortran and C/C++ code.
    setup_requires=[ # these are the packages that must be installed or
        'numpy',
    ],
    install_requires=[
        'scipy',
        'matplotlib',
        'biopython',    # for sequence IO
        'pysqlite',     # for kmer handling with sqlite3
        'termcolor',    # for colored text output
        'cffi',         # for the C component
        'python-igraph',# for bindings to igraph
        'cairocffi',    # for igraph plots
    ],
    extras_require={
        'docs': [
            'sphinx',
            'pycparser', # used to generate C component docs index
            'sphinx_rtd_theme',
            'sphinxcontrib-wiki',
            'breathe',
            'mock', # only used on rtfd.org to mock packages with binary dependencies
        ],
        'tests': [
            'flake8',
            'pytest',
            'tox'
        ]
    }
)
