#!/usr/bin/env python
import sys
from .. import align
from .. import utils
from .. import randseq

params = {
    'm':        3, # match score
    'tr':      -2, # transition score (A<->G and C<->T)
    'tv':      -2, # transversion score
    'go':      -3, # gap open score
    'ge':      -2, # gap extend score
    'band':    -1, # band width if positive
    'show_dp':  0, # whether to print the DP table
    'type': align.ALIGN_END_ANCHORED, # type of alignments
}
mutation_rates = {i:{k:0.7 if k==i else 0.1 for k in 'ACGT'} for i in 'ACGT'}

with open('data/dna.mtrtv.matrix') as f:
    subst_scores = eval(f.read().strip(), params)

#S = align.Sequence('')
#T = align.Sequence('')
S = align.Sequence(randseq.gen(1000))
T, m_transcript = randseq.mutate(str(S), gap_open=0.1, gap_continue=0.3, rates=mutation_rates)
T = align.Sequence(T)

C = align.AlignParams(alphabet='ACGT',subst_scores=subst_scores,
    gap_open_score=params['go'], gap_extend_score=params['ge'],
    max_diversion=params['band'])
P = align.AlignProblem(S=S, T=T, params=C,
    align_type=params['type'])
transcript = P.solve(print_dp_table=params['show_dp'])

print '\n' + transcript + '\n'
if transcript[:3] != 'Err':
    utils.print_alignment(S, T, transcript, sys.stdout, margin=10)
