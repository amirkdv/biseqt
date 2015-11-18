#!/usr/bin/env python
import sys
from math import erf, sqrt, log

from .. import pw, seq
from ..homopolymeric import HpCondenser

params = {
    'gap_prob': 0.05, # gap open score
    'hp_gap_prob': 0.1, # homopolymeric gap probability (linear model)
    'hp_maxlen': 5, # maxlen of the HpCondenser
    'band':    -1, # band width if positive
    'show_condensed_probs': 10,
    'show_dp':  0, # whether to print the DP table
    'type': pw.GLOBAL, # type of alignments
    'subst_probs': [[0.91 if k==i else 0.03 for k in range(4)] for i in range(4)],
}
A = seq.Alphabet('ACGT')

subst_scores = pw.AlignParams.subst_scores_from_probs(A, **params)
go_score, ge_score = pw.AlignParams.gap_scores_from_probs(params['gap_prob'], params['gap_prob'])
# hp_go_score, hp_ge_score = pw.AlignParams.gap_scores_from_probs(params['hp_go_prob'], params['hp_ge_prob'])

Tr = HpCondenser(A, maxlen=params['hp_maxlen'])

subst_probs_d = Tr.condense_subst_probs(**params)
A_d = Tr.dst_alphabet
if params['show_condensed_probs']:
    for idx, row in enumerate(subst_probs_d):
        print A_d.letters[idx], round(sum(row),3), [round(f, 5) for f in row]

S = A.randseq(300, hp_prob=0.1)
T, _ = S.mutate(go_prob=params['gap_prob'], ge_prob=params['gap_prob'], subst_probs=params['subst_probs'])
C = pw.AlignParams(alphabet=A, subst_scores=subst_scores, go_score=go_score, ge_score=ge_score)
with pw.AlignProblem(S=S, T=T, params=C, align_type=params['type']) as P:
    print 'Alignment in original alphabet:'
    P.solve(print_dp_table=params['show_dp'])
    transcript = P.traceback()
    print transcript
    transcript.pretty_print(S, T, sys.stdout)
print

S_d, T_d = Tr.condense_sequence(S), Tr.condense_sequence(T)
subst_scores_d = pw.AlignParams.subst_scores_from_probs(A_d, subst_probs=subst_probs_d, **{k:params[k] for k in params if k != 'subst_probs'})
C_d = pw.AlignParams(alphabet=A_d, subst_scores=subst_scores_d, go_score=go_score, ge_score=ge_score)
with pw.AlignProblem(S=S_d, T=T_d, params=C_d, align_type=params['type']) as P:
    print 'Alignment in condensed alphabet:'
    P.solve(print_dp_table=params['show_dp'])
    transcript_d = P.traceback()
    print transcript_d
    transcript_d.pretty_print(S_d, T_d, sys.stdout)

    if transcript_d:
        expanded_transcript_d = Tr.expand_transcript(S, T, transcript_d)
        print expanded_transcript_d
        print 'Expanded condensed alignment: '
        expanded_transcript_d.pretty_print(S, T, sys.stdout)
        print 'Score of condensed alignment in original alphabet:', C.score(S,T,expanded_transcript_d.opseq)
