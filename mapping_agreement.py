#!/usr/bin/env python

from biseqt.mapping import Mapping
from matplotlib import pyplot as plt
from itertools import combinations

def load_mappings(path):
    with open(path) as f:
        return eval(f.read())

bwa = load_mappings('leishmania/bwa.mappings.txt')
blasr = load_mappings('leishmania/blasr.mappings.txt')
mapped = set(bwa.keys()).intersection(set(blasr.keys()))
overlaps = {}
agree = 0
for i,j in combinations(mapped, 2):
    bwa_overlap = bwa[i].strand == bwa[j].strand and not (bwa[i].ref_from > bwa[j].ref_to or bwa[i].ref_to > bwa[j].ref_from)
    blasr_overlap = blasr[i].strand == blasr[j].strand and not (blasr[i].ref_from > blasr[j].ref_to or blasr[i].ref_to > blasr[j].ref_from)
    agree += 1 if bwa_overlap == blasr_overlap else 0

agree = (200.0*agree/(len(mapped)*(len(mapped)-1)))
print '%%%.2f' % agree, 'of pairs agree among ', len(mapped)
