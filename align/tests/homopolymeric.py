#!/usr/bin/env python
import sys
from math import erf, sqrt, log

from .. import align, seq
from ..homopolymeric import HpCondensor

maxlen = 5

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
C = align.AlignParams(alphabet=A,subst_scores=subst_scores,
    go_score=params['go_score'], ge_score=params['ge_score'],
    max_diversion=params['band'])
Tr = HpCondensor(A, maxlen=5)

S = seq.Sequence('AACCCCCCCCGGGT', A)
T = seq.Sequence('AATCCGGGTTT', A)
S_d = Tr.condense_sequence(S)
T_d = Tr.condense_sequence(T)
print S, S_d
print T, T_d
transcript_d = align.Transcript(raw_transcript='(0,0),0:MISMS')
transcript = Tr.expand_transcript(S, T, transcript_d)
print 'Condensed transcript: ' , transcript_d
transcript_d.pretty_print(S_d, T_d, sys.stdout)
print 'Expanded  transcript: ', transcript
transcript.pretty_print(S, T, sys.stdout)

C_d = Tr.condense_align_params(C, hp_go_score=params['hp_go_score'],
    hp_ge_score=params['hp_ge_score'])
with align.AlignProblem(S=S_d, T=T_d, params=C_d, align_type=params['type']) as P:
    P.solve(print_dp_table=params['show_dp'])
    correct_transcript_d = P.traceback()

if transcript_d:
    print 'Condensed alignment transcript: ', transcript_d
    transcript = Tr.expand_transcript(S, T, transcript_d)
    print 'Expanded correct alignment: '
    transcript.pretty_print(S, T, sys.stdout, margin=10)
