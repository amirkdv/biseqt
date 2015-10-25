#!/usr/bin/env python
import sys
import os

from .. import align, tuples, utils, seq

#DATA_DIR = os.path.join(os.path.dirname(__file__), '../../../pacbio.c-elegans')
#SRC = os.path.join(DATA_DIR, 'm140928_184123_42139_c100719602550000001823155305141590_s1_p0.1.subreads.fasta')
#QUERY = os.path.join(DATA_DIR, 'celegans.query')
#DB = 'celegans.db'

wordlen = 5
limit = -1

SRC = './test.fa'
QUERY = './query.fa'
DB = 'test.db'

A = seq.Alphabet('ACGT')
B = tuples.TuplesDB(db=DB, wordlen=wordlen, alphabet=A)
B.initdb()
B.populate(SRC, lim=limit)
B.index()

params = {
    'm':        3, # match score
    'tr':      -2, # transition score (A<->G and C<->T)
    'tv':      -2, # transversion score
    'go':      -3, # gap open score
    'ge':      -2, # gap extend score
    'band':    -1, # band width if positive
    'show_dp':  0, # whether to print the DP table
    'type': align.ALIGN_GLOBAL, # type of alignments
}
with open('data/dna.mtrtv.matrix') as f:
    subst_scores = eval(f.read().strip(), params)

window = 10
num_seeds = 10 # number of seeds to pursue
C = align.AlignParams(alphabet=A,subst_scores=subst_scores,
    gap_open_score=params['go'], gap_extend_score=params['ge'],
    max_diversion=params['band'])

with open(QUERY) as f:
    query = tuples.Query(f.read().strip(), tuplesdb=B, align_params=C)

hits = query.hitsummary()
cands = hits.keys()
cands.sort(key=lambda k: hits[k], reverse=True)

# purse num_seeds for the best candidate
            # seeds.sort(key=lambda k: k.len, reverse=True)
for seed in query.seeds(seqid=cands[0]):
    if not num_seeds:
        break
    num_seeds -= 1
    transcript = query.expand_seed(seed)
    print transcript
