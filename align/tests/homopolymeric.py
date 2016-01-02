#!/usr/bin/env python
from math import erf, sqrt, log

from .. import pw, seq
from ..homopolymeric import HpCondenser


params = {
    'go_prob': 0.1,        # gap open probability
    'ge_prob': 0.15,         # gap extend probability
    'subst_probs': [[0.97 if k==i else 0.01 for k in range(4)] for i in range(4)],
    'type': pw.OVERLAP, # type of alignments
    # -------- Hp configuration ----------
    'hp_maxlen': 5,
    'hp_gap_score': -0.5,
    # ------ Not used, broken model ------------
    'hp_gap_prob': 0.4, # homopolymeric gap probability (linear model)
    # -----
    'band':    -1, # band width if positive
    'show_dp':  0, # whether to print the DP table
    'show_condensed_probs': 0,
}
A = seq.Alphabet('ACGT')

subst_scores = pw.AlignParams.subst_scores_from_probs(A, **params)
print subst_scores
go_score, ge_score = pw.AlignParams.gap_scores_from_probs(params['go_prob'], params['ge_prob'])
print 'go_score: ', round(go_score,2), 'ge_score: ', round(ge_score,2), 'hp ge_score: ', params['hp_gap_score']
# hp_go_score, hp_ge_score = pw.AlignParams.gap_scores_from_probs(params['hp_go_prob'], params['hp_ge_prob'])

Tr = HpCondenser(A, maxlen=params['hp_maxlen'])

# subst_probs_d = Tr.condense_subst_probs(**params)
# A_d = Tr.dst_alphabet
# subst_scores_d = pw.AlignParams.subst_scores_from_probs(A_d, subst_probs=subst_probs_d, **{k:params[k] for k in params if k != 'subst_probs'})
# C_d = pw.AlignParams(alphabet=A_d, subst_scores=subst_scores_d, go_score=go_score, ge_score=ge_score)

# ================================
from Bio import SeqIO
lookup_S = 'R8139da7f'
lookup_T = 'R12c42542'
# + R12c42542 217456-221530 #97                        +--[2706.00]-->                             R8139da7f 217322-231365 #83
for rec in SeqIO.parse('leishmania/reads.annotated.fa', 'fasta'):
    if lookup_S in rec.name:
        S = str(rec.seq).upper()
        assert(set(S) == set('ACGT'))
        break
    continue
    if lookup_T in rec.name:
        T = str(rec.seq).upper()
        assert(set(T) == set('ACGT'))

#print len(S)
#raise
S = seq.Sequence(S[:6000], A)
#T = seq.Sequence(T[1368:], A)
T = seq.Sequence(str(SeqIO.parse('leishmania/genome.fa', 'fasta').next().seq).upper()[216391:216391+6000].replace('N', 'A'), A)
# S = seq.Sequence(S[S_idx+offset:S_idx+offset+L].replace('N', ''), A)
# T = seq.Sequence(T[offset:offset+L].replace('N', 'A'), A)
#====================================

# S = A.randseq(300, hp_prob=0.1)
# T, _ = S.mutate(go_prob=params['go_prob'], ge_prob=params['ge_prob'], subst_probs=params['subst_probs'])
C = pw.AlignParams(alphabet=A, subst_scores=subst_scores, go_score=go_score, ge_score=ge_score)
for row in C.subst_scores:
    print ', '.join(('%+.2f' % f).rjust(6) for f in row)
with pw.AlignProblem(S=S, T=T, params=C, align_type=params['type']) as P:
    print 'Alignment in original alphabet:'
    P.solve(print_dp_table=params['show_dp'])
    transcript = P.traceback()
    print transcript
    transcript.pretty_print(S, T)
print

S_d, T_d = Tr.condense_sequence(S), Tr.condense_sequence(T)
C_d = Tr.condense_align_params(C, hp_gap_score=params['hp_gap_score'])
for row in C_d.subst_scores:
    print ', '.join(('%+.2f' % f).rjust(6) for f in row)
with pw.AlignProblem(S=S_d, T=T_d, params=C_d, align_type=params['type']) as P:
    print 'Alignment in condensed alphabet:'
    P.solve(print_dp_table=params['show_dp'])
    transcript_d = P.traceback()
    print transcript_d
    transcript_d.pretty_print(S_d, T_d)
    # transcript_d.pretty_print(S_d, T_d)

    print

    if transcript_d:
        expanded_transcript_d = Tr.expand_transcript(S, T, transcript_d)
        expanded_transcript_d.score = C.score(S, T, expanded_transcript_d.opseq)
        print 'Alignment in condensed alphabet expanded to original alphabet: '
        print expanded_transcript_d
        expanded_transcript_d.pretty_print(S, T)
