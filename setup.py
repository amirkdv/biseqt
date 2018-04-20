#!/usr/bin/env python

from setuptools import setup
from os import path, environ

try:
    import numpy
except ImportError:
    raise ImportError('numpy must be already installed')

with open(path.join(path.abspath(path.dirname(__file__)), 'README.md')) as f:
    README = f.read()

install_requires = [
    'scipy',
    'matplotlib',
    'apsw',         # for speaking to SQLite databases
    'termcolor',    # for colored text output
    'cffi',         # for the C component
]

# so we don't install anything on readthedocs.io
if environ.get('READTHEDOCS', None) == 'True':
    install_requires = [
        'alabaster', # our theme
        'sphinxcontrib-wiki', # for breaking down doc pages to sections
        'breathe',   # for doxygen integration, doxygen is executed in conf.py
    ]

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
    # NOTE numpy has some weirdness with setuptools; cf. open issue:
    #   https://github.com/numpy/numpy/issues/2434
    #
    # We need it because it's a dependency of scipy but has to be installed
    # either manually before our setup.py or declared as setup_requires, cf.
    #   https://github.com/numpy/numpy/issues/2434#issuecomment-65252402.
    #
    # However, it is still much faster to install numpy separately via pip
    # because installing by setup_requires wants to compile all the fortran and
    # C/C++ code.
    setup_requires=[ # for setup.py to work at all
        'numpy',
    ],
    install_requires=install_requires,
    extras_require={
        'docs': [
            'sphinx',
            'pycparser', # for generating C component docs index, cf. cdocs.py
            'alabaster', # our theme
            'sphinxcontrib-wiki>=0.5.0', # for breaking down doc pages to sections
            'breathe',   # for doxygen integration
            'mock',      # for rtfd.org to skip installing system dependencies
        ],
        'tests': [
            'mock',      # for mocking external functions
            'flake8',    # for style enforcement
            'pytest',    # to run and manage tests
            'pytest-cov',# to get coverage reports from tests
            'coverage-badge', # to generate the coverage badge, cf. circle.yml
        ]
    }
)
