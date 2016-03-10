#!/usr/bin/env python
from .. import pw, seq

params = {
    'length': 500,
    'go_prob': 0.1, # gap open score
    'ge_prob': 0.3, # gap extend score
    'show_dp':  0, # whether to print the DP table
    'alntype': pw.GLOBAL, # type of alignments
    'subst_probs': [[0.91 if k==i else 0.03 for k in range(4)] for i in range(4)]
}
A = seq.Alphabet('ACGT')
subst_scores = pw.AlignScores.subst_scores_from_probs(A, gap_prob=params['go_prob'], **params)
go_score, ge_score = pw.AlignScores.gap_scores_from_probs(params['go_prob'], params['ge_prob'])

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
T, m_opseq = S.mutate(**params)
C = pw.AlignScores(subst_scores=subst_scores, alphabet=A,
    go_score=go_score, ge_score=ge_score)
F = pw.AlignFrame(S, T)
with pw.AlignTable(F, C, alntype=params['alntype']) as P:
    score = P.solve(print_dp_table=params['show_dp'])
    transcript = P.traceback()

print '\n--> optimal alignment:\n%s\n' % str(transcript)
if transcript:
    transcript.pretty_print(S, T)

m_transcript = pw.Transcript(S_idx=0, T_idx=0, score=P.score(m_opseq), opseq=m_opseq)
print '\n--> mutation transcript:\n%s\n' % str(m_transcript)
m_transcript.pretty_print(S, T)
