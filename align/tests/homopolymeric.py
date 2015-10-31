#!/usr/bin/env python
import sys
from math import erf, sqrt, log

from .. import align, seq
from ..homopolymeric import HpCondensor

maxlen = 5

Tr = HpCondensor(maxlen=5)

params = {
    'm':        3, # match score
    'tr':      -2, # transition score (A<->G and C<->T)
    'tv':      -2, # transversion score
    'go_score':      -3, # gap open score
    'ge_score':      -2, # gap extend score
    'hp_go_score':    0, # homopolymeric gap open score
    'hp_ge_score': -0.5, # homopolymeric gap extend score
    'band':    -1, # band width if positive
    'show_dp':  0, # whether to print the DP table
    'type': align.ALIGN_GLOBAL, # type of alignments
}

with open('data/dna.mtrtv.matrix') as f:
    subst_scores = eval(f.read().strip(), params)

A = seq.Alphabet('ACGT')
scores = Tr.translate_subst_scores(subst_scores, go_score=params['go_score'],
    ge_score=params['ge_score'], hp_go_score=params['hp_go_score'],
    hp_ge_score=params['hp_ge_score'])
alphabet = [x for _,_,x in Tr.translate_alphabet('ACGT')]
A_d = seq.Alphabet(alphabet)

s = 'AACCCCGGGGGGGGGGGGT'
t = 'AATGGGGGGGGGTTT'
s_d = Tr.condense(s)
t_d = Tr.condense(t)
S = seq.Sequence(s, A)
T = seq.Sequence(t, A)
print S
print T
opseq_d = 'BMSSS'
opseq = Tr.expand_opseq(s, t, opseq_d)
print opseq_d
print opseq
seq.Transcript('(0,0),0:%s' % opseq).pretty_print(S, T, sys.stdout)
S_d = seq.Sequence(s_d, A_d)
T_d = seq.Sequence(t_d, A_d)
C = align.AlignParams(alphabet=A_d,subst_scores=scores,
    go_score=params['go_score'], ge_score=params['ge_score'],
    max_diversion=params['band'])
P = align.AlignProblem(S=S_d, T=T_d, params=C,
    align_type=params['type'])
P.solve(print_dp_table=params['show_dp'])
transcript_d = P.traceback()
print "opseq score: %.2f" % P.score(opseq_d)

if transcript_d:
    print '\n' + str(transcript_d) + '\n'
    opseq = Tr.expand_opseq(S, T, transcript_d.opseq)
    transcript = seq.Transcript('(%d,%d),0:%s' % ( \
        transcript_d.S_idx * A_d.letter_length,
        transcript_d.T_idx * A_d.letter_length,
        opseq_d)
    )
    transcript.pretty_print(S, T, sys.stdout, margin=10)
