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
    # FIXME max_decr and decr_def should somehow be related statistically to
    # guarantee something.
    'max_decr': 3, # max number of decreasing windows to drop a seed
    'decr_def': -5, # what consistutes a "decrease" in windows of seed extension
}
with open('data/dna.mtrtv.matrix') as f:
    subst_scores = eval(f.read().strip(), params)

C = align.AlignParams(alphabet=A, subst_scores=subst_scores,
    go_score=params['go_score'], ge_score=params['ge_score'],
    max_diversion=params['band'])

import networkx as nx
import matplotlib.pyplot as plt
G = nx.DiGraph()

score_to_report = 80

seqids = B.seqids()
#seqids = seqids[:2]
for idx in range(len(seqids)):
    for t_idx in range(idx+1, len(seqids)):
        S = B.loadseq(seqids[idx])
        T = B.loadseq(seqids[t_idx])
        P = align.AlignProblem(S=S, T=T, params=C, align_type=align.ALIGN_OVERLAP)
        score = P.solve()
        if score > 60:
            sys.stdout.write('* S=%s, T=%s: ' % (str(idx).rjust(2), str(t_idx).rjust(2)))
            transcript = P.traceback()
            if transcript.idx_T == 0:
                sys.stdout.write(' S --(%.2f)--> T ' % score)
            if transcript.idx_S == 0:
                sys.stdout.write(' T --(%.2f)--> S ' % score)
            sys.stdout.write('\n')
        F = tuples.OverlapFinder(S, T, C)
        exacts = B.exactly_matching_segments(seqids[idx], seqids[t_idx])
        if exacts:
            segments = F.extend(exacts, max_decr=params['max_decr'], decr_def=params['decr_def'])
            if segments:
                sys.stdout.write('  found potential overlap:')
                if segments[0][0].idx_T == 0:
                    G.add_edge(idx, t_idx, weight=segments[0][1])
                    sys.stdout.write(' S --(%0.2f)--> T' % segments[0][1])
                if segments[0][0].idx_S == 0:
                    G.add_edge(t_idx, idx, weight=segments[0][1])
                    sys.stdout.write(' T --(%0.2f)--> S' % segments[0][1])
                sys.stdout.write('\n')
                print '  ' + str(segments[0])
                #print 'len(S) = %d, len(T) = %d, alignment length: %d' % (len(S), len(T), len(transcript.opseq)-1)
                #transcript.pretty_print(S, T, sys.stdout, margin=0)

pos = nx.circular_layout(G)
plt.figure(figsize=(50,50))
nx.draw_networkx_nodes(G, pos, node_size=2000, node_color='k')
nx.draw_networkx_labels(G, pos, font_size=30, font_color='w')
nx.draw_networkx_edges(G, pos, width=5)
nx.draw_networkx_edge_labels(G, pos, font_size=26, edge_labels={(f,t):w['weight'] for f,t,w in G.edges(data=True)})
plt.savefig('overlap.svg')
