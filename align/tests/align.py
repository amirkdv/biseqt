#!/usr/bin/env python
import sys
from .. import align, seq

params = {
    'length': 500,
    'go_prob': 0.1, # gap open score
    'ge_prob': 0.3, # gap extend score
    'band':    -1, # band width if positive
    'show_dp':  0, # whether to print the DP table
    'type': align.ALIGN_GLOBAL, # type of alignments
    'subst_probs': [[0.91 if k==i else 0.03 for k in range(4)] for i in range(4)]
}
A = seq.Alphabet('ACGT')
subst_scores = align.AlignParams.subst_scores_from_probs(params['subst_probs'], A)
go_score, ge_score = align.AlignParams.gap_scores_from_probs(params['go_prob'], params['ge_prob'])

print 'Substitution probabilities:'
for i in params['subst_probs']:
    print i
print 'Substitution scores:'
for i in subst_scores:
    print [round(f,2) for f in i]

print 'Pr(go) = %.2f, Pr(ge) = %.2f +----> Score(go)=%.2f, Score(ge)=%.2f' % \
    (params['go_prob'], params['ge_prob'], go_score, ge_score)


#S = seq.Sequence('AGTA', A)
#T = seq.Sequence('GTCGAGT', A)

S = A.randseq(params['length'])
T, m_opseq = S.mutate(go_prob=params['go_prob'], ge_prob=params['ge_prob'],
    subst_probs=params['subst_probs'])
C = align.AlignParams(subst_scores=subst_scores, alphabet=A,
    go_score=go_score, ge_score=ge_score, max_diversion=params['band'])
P = align.AlignProblem(S=S, T=T, params=C,
    align_type=params['type'])
score = P.solve(print_dp_table=params['show_dp'])
transcript = P.traceback()
print 'optimal alignment has score %.2f' % score

print '\n--> optimal alignment:\n%s\n' % str(transcript)
if transcript:
    transcript.pretty_print(S, T, sys.stdout, margin=10)

m_transcript = align.Transcript('(0,0),%.2f:%s' % (P.score(m_opseq), m_opseq))
print '\n--> mutation transcript:\n%s\n' % str(m_transcript)
m_transcript.pretty_print(S, T, sys.stdout, margin=10)
