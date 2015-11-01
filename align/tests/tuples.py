#!/usr/bin/env python
import sys
import os

from .. import align, tuples, seq, assembly

wordlen = 10
A = seq.Alphabet('ACGT')
DB = tuples.TuplesDB(db='genome.db', wordlen=wordlen, alphabet=A)

params = {
    'm':        3, # match score
    'tr':      -2, # transition score (A<->G and C<->T)
    'tv':      -2, # transversion score
    'go_score':      -3, # gap open score
    'ge_score':      -2, # gap extend score
    'band':    -1, # band width if positive
    # FIXME max_decr and decr_def should somehow be related statistically to
    # guarantee something.
    'max_decr': 3, # max number of decreasing windows to drop a seed
    'decr_def': -10, # what consistutes a "decrease" in windows of seed extension
}
with open('data/dna.mtrtv.matrix') as f:
    subst_scores = eval(f.read().strip(), params)

C = align.AlignParams(alphabet=A, subst_scores=subst_scores,
    go_score=params['go_score'], ge_score=params['ge_score'],
    max_diversion=params['band'])

print 'creating overlap graph by tuple extension'
G = assembly.overlap_graph_by_tuple_extension(
    DB, C, max_decr=params['max_decr'], decr_def=params['decr_def']
)
assembly.save_overlap_graph(G, 'overlap_tuple.svg')

print 'creating overlap graph by brute force overlap alignment'
G = assembly.overlap_graph_by_alignment(DB, C, min_score=60)
assembly.save_overlap_graph(G, 'overlap_align.svg')
