#!/bin/bash
set -e

pip install numpy
pip install pysam
pip install -e .
pip install -e .[docs]
pip install -e .[tests]
