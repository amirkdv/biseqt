#!/usr/bin/env python
import sys
from .. import align, utils, seq

params = {
    'go_prob': 0.1, # gap open score
    'ge_prob': 0.3, # gap extend score
    'band':    -1, # band width if positive
    'show_dp':  0, # whether to print the DP table
    'type': align.ALIGN_GLOBAL, # type of alignments
    'subst_probs': [[0.7 if k==i else 0.1 for k in range(4)] for i in range(4)]
}
A = seq.Alphabet('ACGT')
subst_scores = align.AlignParams.subst_scores_from_probs(params['subst_probs'], A)
go_score, ge_score = align.AlignParams.gap_scores_from_probs(params['go_prob'], params['ge_prob'])

for i in params['subst_probs']:
    print i
for i in subst_scores:
    print i

print '---'
print params['go_prob'], params['ge_prob']
print go_score, ge_score

#S = align.Sequence('ACCCGT', A)
#T = align.Sequence('CC', A)
S = A.randseq(100)
T, m_transcript = S.mutate(go_prob=params['go_prob'], ge_prob=params['ge_prob'],
    subst_probs=params['subst_probs'])
C = align.AlignParams(subst_scores=subst_scores, alphabet=A,
    go_score=go_score, ge_score=ge_score, max_diversion=params['band'])
P = align.AlignProblem(S=S, T=T, params=C,
    align_type=params['type'])
transcript = P.solve(print_dp_table=params['show_dp'])

print '\n' + transcript + '\n'
if transcript[:3] != 'Err':
    utils.print_alignment(S, T, transcript, sys.stdout, margin=10)

print 'mutation opseq score: %.2f' % P.score(m_transcript)
