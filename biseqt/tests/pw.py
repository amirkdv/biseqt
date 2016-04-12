#!/usr/bin/env python
from .. import pw, seq

params = {
    'go_prob': 0.05, # gap open score
    'ge_prob': 0.1, # gap extend score
    'subst_probs': [[0.91 if k==i else 0.03 for k in range(4)] for i in range(4)]
}
import os
A = seq.Alphabet('ACGT')
subst_scores = pw.AlignScores.subst_scores_from_probs(A, gap_prob=params['go_prob'], **params)
#subst_scores = [[+1 if k == i else float(os.environ['SCORE']) for k in range(4)] for i in range(4)]
go_score, ge_score = pw.AlignScores.gap_scores_from_probs(params['go_prob'], params['ge_prob'])

print 'Substitution probabilities:'
for i in params['subst_probs']:
    print i
print 'Substitution scores:'
for i in subst_scores:
    print [round(f,2) for f in i]

print 'Pr(go) = %.2f, Pr(ge) = %.2f +----> Score(go)=%.2f, Score(ge)=%.2f' % \
    (params['go_prob'], params['ge_prob'], go_score, ge_score)


#S = seq.Sequence('GGCCGAACCCAGTCGTACGTCTCCTTGTGAAGTATAAGCTGCATATAGATCATTGATAAAGATTTAGGTAGCTGACAAGCCCCGGAGGCATGGCTTGCAT', A)
#T = seq.Sequence('CAGAGACAGCGCCAGAGAGACCAGCCAGAA', A)
S = A.randseq(100)
T, m_opseq = S.mutate(**params)
C = pw.AlignScores(subst_scores=subst_scores, alphabet=A,
    go_score=go_score, ge_score=ge_score)
F = pw.AlignFrame(S, T)
with pw.AlignTable(F, C, alnmode=pw.STD_MODE, alntype=pw.GLOBAL) as P:
    score = P.solve()
    transcript = P.traceback()

print '\n--> optimal global alignment:\n%s\n' % str(transcript)
if transcript:
    transcript.pretty_print(S, T)

S = A.randseq(8000)
T, m_opseq = seq.Sequence(S[4000:] + A.randstr(1000), A).mutate(**params)
F = pw.AlignFrame(S, T)
with pw.AlignTable(F, C, alnmode=pw.BANDED_MODE, alntype=pw.B_OVERLAP, dmin=3800, dmax=4200) as P:
    score = P.solve()
    transcript = P.traceback()

print '\n--> optimal overlap alignment (banded):\n%s\n' % str(transcript)
if transcript:
    transcript.pretty_print(S, T)
