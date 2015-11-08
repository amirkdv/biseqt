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

S = seq.Sequence('AACCCCGGGGGGGGGGGGT', A)
T = seq.Sequence('AATGGGGGGGGGTTT', A)
print S
print T
S_d = Tr.condense(S)
T_d = Tr.condense(T)
opseq_d = 'MSSS'
opseq = Tr.expand_opseq(S, T, opseq_d)
print opseq_d
print opseq
align.Transcript(raw_transcript='(0,0),0:%s' % opseq).pretty_print(S, T, sys.stdout)

C_d = Tr.translate_align_params(C, hp_go_score=params['hp_go_score'],
    hp_ge_score=params['hp_ge_score'])
with align.AlignProblem(S=S_d, T=T_d, params=C_d, align_type=params['type']) as P:
    P.solve(print_dp_table=params['show_dp'])
    transcript_d = P.traceback()
    print "opseq score: %.2f" % P.score(opseq_d)

if transcript_d:
    print '\n' + str(transcript_d) + '\n'
    opseq = Tr.expand_opseq(S, T, transcript_d.opseq)
    transcript = align.Transcript(
        idx_S=transcript_d.idx_S * Tr.dst_alphabet.letter_length,
        idx_T=transcript_d.idx_T * Tr.dst_alphabet.letter_length, score=0, opseq=opseq_d)
    transcript.pretty_print(S, T, sys.stdout, margin=10)
