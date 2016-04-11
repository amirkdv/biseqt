#!/usr/bin/env python
import sys
import os
import igraph
from .. import pw, words, seq, overlap, overlap, mapping, ProgressIndicator

params = {
    'show_params': False,   # print a summary of parameters
    'rw_collect': False,      # whether to dump score random walks to 'scores.txt'
    'wordlen': int(os.environ['WORDLEN']), # 10,          # tuple word length for seed extension
    'go_prob': 0.1,        # gap open probability
    'ge_prob': 0.15,         # gap extend probability
    'subst_probs': [[0.97 if k==i else 0.01 for k in range(4)] for i in range(4)],
    # ------------ Assembly ---------------
    'min_margin': 500,       # minimum margin required for the direction to be reliable.
    'window': 50,            # rolling window length for tuple extension.
    'max_new_mins': 5,       # how many consecutive drops are allowed.
    'min_overlap_score': 500,# minimum score required for an overlap to be reported.
    'shift_rolling_sum_width': 100,
    'lower_log_pvalue_cutoff': -200, # FIXME
    'upper_log_pvalue_cutoff': 0, # FIXME
    # ------------- Index ----------------
    'min_seeds_for_homology': 1, # minimum number of seeds for two reads to be considered.
    'min_word_log_pvalue': 0, # minimum p-value allowed for a word to be considered a seed
    # ---------------- simulations ------------------
    'genome_length': 70000, # length of randomly generated genome
    'coverage': 10,          # coverage of random sequencing reads
    'read_len_mean': 2000,  # average length of sequencing read
    'read_len_sd':  200,    # variance of sequencing read length
}

A = seq.Alphabet('ACGT')

go_score, ge_score = pw.AlignScores.gap_scores_from_probs(params['go_prob'], params['ge_prob'])
subst_scores = pw.AlignScores.subst_scores_from_probs(A, **params)
C = pw.AlignScores(
    alphabet=A,
    subst_scores=subst_scores,
    go_score=go_score,
    ge_score=ge_score
)
od_params = overlap.OverlapDiscoveryParams(
    seed_ext_params=overlap.SeedExtensionParams(
        window=params['window'],
        min_score=params['min_overlap_score'],
        max_new_mins=params['max_new_mins'],
        scores=C
    ),
    shift_rolling_sum_width=params['shift_rolling_sum_width'],
)
db_id = lambda G, vid: int(G.vs[vid]['name'].split('#')[1])

def show_params():
    print 'Substitution probs:'
    for i in params['subst_probs']:
        print [round(f,2) for f in i]
    print '\nSubstitution scores:'
    for i in subst_scores:
        print [round(f,2) for f in i]
    print '\nGap probs\n  go = %.2f, ge = %.2f\n\nGap scores\n  go=%.2f, ge=%.2f' % \
        (params['go_prob'], params['ge_prob'], go_score, ge_score)

    print '\nmax_new_mins = %d, window = %d' % \
        (params['max_new_mins'], params['window'])

def create_example(db, reads='reads.fa'):
    seq.make_sequencing_fixture('genome.fa', reads,
        genome_length=params['genome_length'],
        coverage=params['coverage'],
        len_mean=params['read_len_mean'],
        len_sd=params['read_len_sd'],
        subst_probs=params['subst_probs'],
        ge_prob=params['ge_prob'],
        go_prob=params['go_prob'],
    )

def true_overlaps(true_path):
    G = igraph.read(true_path)
    return [set([G.vs[u]['name'], G.vs[v]['name']]) for u, v in G.get_edgelist()]

def create_denovo_db(db, reads):
    B = seq.SeqDB(db, alphabet=A)
    B.initdb()
    B.populate(reads, seq_type=seq.READ)

    Idx = words.Index(seqdb=B, **params)
    Idx.initdb()
    Idx.index()

def build_denovo_overlap_graph(db, path, true_path):
    B = seq.SeqDB(db, alphabet=A)
    Idx = words.Index(seqdb=B, **params)
    if params['show_params']:
        show_params()
    G = overlap.discovery.overlap_graph(Idx, od_params, min_margin=params['min_margin'], rw_collect=params['rw_collect'])
    G.save(path)

def plot_word_pvalues(db, path):
    B = seq.SeqDB(db, alphabet=A)
    Idx = words.Index(seqdb=B, **params)
    words.plot_word_pvalues(Idx, path)

def plot_shift_pvalues(db, path, true_path):
    B = seq.SeqDB(db, alphabet=A)
    Idx = words.Index(seqdb=B, **params)
    G = igraph.read(true_path)
    overlap.plot_shift_signifiance_discrimination(
        path,
        Idx,
        true_overlaps(true_path),
        num_bins=500
    )

def plot_num_seeds(db, path, true_path):
    B = seq.SeqDB(db, alphabet=A)
    Idx = words.Index(seqdb=B, **params)
    overlap.plot_num_seeds_discrimination(path, Idx, true_overlaps(true_path))

def plot_seeds(db, path, true_path, mappings):
    B = seq.SeqDB(db, alphabet=A)
    Idx = words.Index(seqdb=B, **params)
    G = igraph.read(true_path)
    true_overlaps = [set(G.vs[u]['name'], G.vs[v]['name']) for u, v in G.get_edgelist() if G]
    with open(mappings) as f:
        mappings = eval(f.read())
    overlap.plot_all_seeds(
        Idx,
        params['shift_rolling_sum_width'],
        basedir=path,
        true_overlaps=true_overlaps,
        mappings=mappings
    )

def plot_rw(db, path, true_path):
    B = seq.SeqDB(db, alphabet=A)
    overlap.plot_seed_extension_rws(
        path, B.seqinfo(), max_rws=225, draw_type='-+',
        logfile='scores.txt', true_overlaps=true_overlaps(true_path)
    )
