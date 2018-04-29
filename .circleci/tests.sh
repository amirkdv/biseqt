#!/bin/bash
set -e # fail on any failing line
set -u # fail on any undefined variable

pip install numpy
pip install pysam
pip install -e .
pip install -e .[docs]
pip install -e .[tests]

make docs
make tests
