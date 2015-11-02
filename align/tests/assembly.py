#!/usr/bin/env python
import sys
import os

from .. import align, tuples, seq, assembly

A = seq.Alphabet('ACGT')

params = {
    'wordlen': 3,           # tuple word lengths
    'genome_length': 50,    # length of randomly generated genome
    'coverage': 3,          # coverage of random sequencing reads
    'read_len_mean': 20,    # average length of sequencing read
    'read_len_var': 1,      # variance of sequencing read length
    'min_align_score': 5,   # minimum overlap alignment score to constitue an edge
    'go_prob': 0.0001,      # gap open score
    'ge_prob': 0.0003,      # gap extend score
    'subst_probs': [[0.999 if k==i else 0.0001 for k in range(4)] for i in range(4)],
    # FIXME max_decr and decr_def should somehow be related statistically to
    # guarantee something.
    'drop_threshold': -10   # maximum total score drop allowed before dropping a seed
}
subst_scores = align.AlignParams.subst_scores_from_probs(params['subst_probs'], A)
go_score, ge_score = align.AlignParams.gap_scores_from_probs(params['go_prob'], params['ge_prob'])
C = align.AlignParams(
    alphabet=A, subst_scores=subst_scores,
    go_score=go_score, ge_score=ge_score
)

def show_params():
    print 'Substitution probabilities:'
    for i in params['subst_probs']:
        print i
    print 'Substitution scores:'
    for i in subst_scores:
        print [round(f,2) for f in i]
    print 'Pr(go) = %.2f, Pr(ge) = %.2f +----> Score(go)=%.2f, Score(ge)=%.2f' % \
        (params['go_prob'], params['ge_prob'], go_score, ge_score)

def create_example(db):
    seq.make_sequencing_fixture('genome.fa', 'reads.fa',
        genome_length=params['genome_length'],
        coverage=params['coverage'],
        len_mean=params['read_len_mean'],
        len_var=params['read_len_var'],
        subst_probs=params['subst_probs'],
        ge_prob=params['ge_prob'],
        go_prob=params['go_prob']
    )
    B = tuples.TuplesDB(db, wordlen=params['wordlen'], alphabet=A)
    B.initdb()
    B.populate('reads.fa');
    B.index()

def overlap_by_tuple_extension(db, path):
    B = tuples.TuplesDB(db, wordlen=params['wordlen'], alphabet=A)
    G = assembly.overlap_graph_by_tuple_extension(B, C, params['drop_threshold'])
    assembly.save_overlap_graph(G, path)

def overlap_by_alignment(db, path):
    B = tuples.TuplesDB(db, wordlen=params['wordlen'], alphabet=A)
    G = assembly.overlap_graph_by_alignment(B, C, min_score=params['min_align_score'])
    assembly.save_overlap_graph(G, path)

def overlap_by_known_order(db, path):
    B = tuples.TuplesDB(db, wordlen=params['wordlen'], alphabet=A)
    G = assembly.overlap_graph_by_known_order(B)
    assembly.save_overlap_graph(G, path)
