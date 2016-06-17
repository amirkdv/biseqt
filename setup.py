#!/usr/bin/env python

from setuptools import setup
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md')) as f:
    long_description = f.read()

# FIXME manual instructions:
# * cffi:
#   * backend: apt-get install python-cffi libffi-dev
#   * frontend: pip install cffi==A.B.C where A.B.C comes from: apt-cache policy python-cffi
# * igraph:
#   * apt-get install python-igraph (via pip requires c dependencies)
# * matplotlib:
#   * apt-get install python-matplotlib (via pip requires many dependencies)
#   * apt-get install dvipng (for latex in images)

setup(
    name='biseqt',
    version='0.0.1',
    description='Oval is an overlap alignment search tool',
    long_description=long_description,
    url='https://github.com/amirkdv/biseqt',
    author='Amir Kadivar',
    author_email='amir@amirkdv.ca',
    #license='MIT',
    #keywords='sequence DNA overlap alignment fragment assembly',
    packages=['biseqt'],
    setup_requires=[
        'cffi',
        'biopython',
    ],
    install_requires=[
        'biopython',
        'termcolor',
        'cffi',
        'python-igraph',
        # 'pycairo', can't be done yet; see https://bugs.freedesktop.org/show_bug.cgi?id=58772
    ],
    extras_require={
        'DOCS': [
            'pycparser', # used to generate C component docs index
            'sphinx',
            'sphinx.ext.napoleon',
            'sphinx.rdt.theme',
            'sphinx.ext.mathjax',
            'breathe',
            'mock' # only used on rtfd.org to mock packages with binary dependencies
        ]
    }
)
