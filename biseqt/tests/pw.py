#!/usr/bin/env python
from .. import pw, seq

params = {
    'go_prob': 0.1, # gap open score
    'ge_prob': 0.25, # gap extend score
    'subst_probs': [[0.94 if k==i else 0.02 for k in range(4)] for i in range(4)]
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
S = A.randseq(1000)
T, m_opseq = S.mutate(**params)
C = pw.AlignScores(subst_scores=subst_scores, alphabet=A,
    go_score=go_score, ge_score=ge_score)
F = pw.AlignFrame(S, T)
with pw.AlignTable(F, C, alnmode=pw.STD_MODE, alntype=pw.GLOBAL) as P:
    score = P.solve()
    transcript = P.traceback()
    #P.rasterplot('debug.png', transcript)
    #P.rasterplot('debug_true.png', pw.Transcript(S_idx=0, T_idx=0, opseq=m_opseq, score=0))
    #raise

print '\n--> optimal global alignment:\n%s\n' % str(transcript)
if transcript:
    transcript.pretty_print(S, T)

S = A.randseq(8000)
T, _ = seq.Sequence(S[4000:] + A.randstr(1000), A).mutate(**params)
#T = A.randseq(3000)
F = pw.AlignFrame(S, T)
with pw.AlignTable(F, C, alnmode=pw.BANDED_MODE, alntype=pw.B_OVERLAP, dmin=2800, dmax=5200) as P:
#with pw.AlignTable(F, C, alnmode=pw.STD_MODE, alntype=pw.OVERLAP) as P:
    score = P.solve()
    transcript = P.traceback()
    #P.rasterplot('debug.png', transcript)
    #raise

print '\n--> optimal overlap alignment (banded):\n%s\n' % str(transcript)
if transcript:
    transcript.pretty_print(S, T)
