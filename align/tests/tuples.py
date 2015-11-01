#!/usr/bin/env python
import sys
import os

from .. import align, tuples, seq
from Bio import SeqIO

wordlen = 10
limit = -1

DB = 'genome.db'

A = seq.Alphabet('ACGT')
B = tuples.TuplesDB(db=DB, wordlen=wordlen, alphabet=A)

params = {
    'm':        3, # match score
    'tr':      -2, # transition score (A<->G and C<->T)
    'tv':      -2, # transversion score
    'go_score':      -3, # gap open score
    'ge_score':      -2, # gap extend score
    'band':    -1, # band width if positive
}
with open('data/dna.mtrtv.matrix') as f:
    subst_scores = eval(f.read().strip(), params)

C = align.AlignParams(alphabet=A, subst_scores=subst_scores,
    go_score=params['go_score'], ge_score=params['ge_score'],
    max_diversion=params['band'])

import networkx as nx
import matplotlib.pyplot as plt
G = nx.Graph()

seqids = B.seqids()
#seqids = seqids[:2]
for idx in range(len(seqids)):
    for t_idx in range(idx+1, len(seqids)):
        S = B.loadseq(seqids[idx])
        T = B.loadseq(seqids[t_idx])
        F = tuples.OverlapFinder(S, T, C)
        exacts = B.exactly_matching_segments(seqids[idx], seqids[t_idx])
        if exacts:
            segments = F.extend(exacts)
            if segments:
                print 'found %d overlap candidates for seqs (%d, %d)\n      %s' % (len(segments), idx, t_idx, segments[0])
                G.add_edge(idx, t_idx, weight=segments[0][1])

pos = nx.circular_layout(G)
plt.figure(figsize=(50,50))
nx.draw_networkx_nodes(G, pos, node_size=2000, node_color='k')
nx.draw_networkx_labels(G, pos, font_size=30, font_color='w')
nx.draw_networkx_edges(G, pos, width=5)
nx.draw_networkx_edge_labels(G, pos, font_size=26, edge_labels={(f,t):w['weight'] for f,t,w in G.edges(data=True)})
plt.savefig('overlap.svg')
