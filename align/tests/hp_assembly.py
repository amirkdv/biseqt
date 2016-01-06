#!/usr/bin/env python
import sys
import os
import igraph

from .. import pw, tuples, seq, overlap, homopolymeric, overlap

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
    'lower_log_pvalue_cutoff': -200, # -10000,
    'upper_log_pvalue_cutoff': 0, # -800,
    # ------------- HP / Index ----------------
    'min_seeds_for_homology': 1, # minimum number of seeds for two reads to be considered.
    'min_word_log_pvalue': -5, # minimum p-value allowed for a word to be considered a seed
    'hp_gap_score': -0.5,   # HpCondenser Hp gap score
    'hp_maxlen': 5,         # HpCondenser maxlen for seed extension
    # ---------------- simulations ------------------
    'genome_length': 70000, # length of randomly generated genome
    'coverage': 10,          # coverage of random sequencing reads
    'read_len_mean': 2000,  # average length of sequencing read
    'read_len_sd':  200,    # variance of sequencing read length
    'hp_prob': 0.15,        # Parameter for geometric distributions of hp stretches
}

A = seq.Alphabet('ACGT')
Tr = homopolymeric.HpCondenser(A, maxlen=params['hp_maxlen'])

go_score, ge_score = pw.AlignParams.gap_scores_from_probs(params['go_prob'], params['ge_prob'])
subst_scores = pw.AlignParams.subst_scores_from_probs(A, **params)
C = pw.AlignParams(
    alphabet=A,
    subst_scores=subst_scores,
    go_score=go_score,
    ge_score=ge_score
)
C_d = Tr.condense_align_params(C, hp_gap_score=params['hp_gap_score'])
subst_scores_d = C_d.subst_scores
db_id = lambda G, vid: int(G.vs[vid]['name'].split('#')[1])

def show_params():
    if not params['show_params']:
        return
    print 'Substitution probabilities:'
    print 'Substitution scores:'
    for i in subst_scores_d:
        print [round(f,2) for f in i]
    print 'Pr(go) = %.2f, Pr(ge) = %.2f +----> Score(go)=%.2f, Score(ge)=%.2f' % \
        (params['go_prob'], params['ge_prob'], go_score, ge_score)

    print 'max_new_mins = %d, window = %d' % \
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
        #hp_prob=params['hp_prob']
    )

def true_overlaps(true_path):
    G = igraph.read(true_path)
    return [set([db_id(G, u), db_id(G, v)]) for u, v in G.get_edgelist()]

def create_db(db, reads):
    B = tuples.TuplesDB(db, alphabet=A)
    Idx = tuples.Index(tuplesdb=B, **params)
    B.initdb()
    B.populate(reads);
    Idx.initdb()
    Idx.index()

def plot_word_pvalues(db, path):
    B = tuples.TuplesDB(db, alphabet=A)
    Idx = tuples.Index(tuplesdb=B, **params)
    Idx.plot_word_pvalues(path)

def plot_shift_pvalues(db, path, true_path):
    B = tuples.TuplesDB(db, alphabet=A)
    Idx = tuples.Index(tuplesdb=B, **params)
    G = igraph.read(true_path)
    overlap.plot_shift_signifiance_discrimination(
        path,
        Idx,
        params['shift_rolling_sum_width'],
        true_overlaps(true_path),
        num_bins=500
    )

def plot_num_seeds(db, path, true_path):
    B = tuples.TuplesDB(db, alphabet=A)
    Idx = tuples.Index(tuplesdb=B, **params)
    overlap.plot_num_seeds_discrimination(path, Idx, true_overlaps(true_path))

def plot_seeds(db, path, true_path):
    B = tuples.TuplesDB(db, alphabet=A)
    Idx = tuples.Index(tuplesdb=B, **params)
    G = igraph.read(true_path)
    true_overlaps = [set([db_id(G, u), db_id(G, v)]) for u, v in G.get_edgelist()]
    overlap.plot_all_seeds(
        Idx,
        params['shift_rolling_sum_width'],
        basedir=path,
        true_overlaps=true_overlaps
    )

def plot_rw(db, path, true_path):
    B = tuples.TuplesDB(db, alphabet=A)
    overlap.plot_seed_extension_rws(
        path, B.seqinfo(), max_rws=225, draw_type='-+',
        logfile='scores.txt', true_overlaps=true_overlaps(true_path)
    )

def build_overlap_graph(db, path, true_path):
    B = tuples.TuplesDB(db, alphabet=A)
    Idx = tuples.Index(tuplesdb=B, **params)
    show_params()
    # FIXME:
    G = overlap.build_overlap_graph(index=Idx, true_overlaps=true_overlaps(true_path),
        align_params=C_d, hp_condenser=Tr, **params)
    # G = overlap.build_overlap_graph(Idx, true_overlaps=true_overlaps(true_path)
    #     align_params=C, **params)
    G.save(path)

def overlap_graph_by_known_order(db):
    """Builds the *correct* weighted, directed graph by using hints left in
    reads databse by ``seq.make_sequencing_fixture()``.

    Args:
        tuplesdb (tuples.TuplesDB): The tuples database.

    Returns:
        networkx.DiGraph
    """
    B = tuples.TuplesDB(db, alphabet=A)
    seqinfo = B.seqinfo()
    seqids = seqinfo.keys()
    vs = set()
    es, ws = [], []
    for sid_idx in range(len(seqids)):
        for tid_idx in range(sid_idx + 1, len(seqids)):
            S_id, T_id = seqids[sid_idx], seqids[tid_idx]
            S_start, T_start = seqinfo[S_id]['start'], seqinfo[T_id]['start']
            S_end = S_start + seqinfo[S_id]['length']
            T_end = T_start + seqinfo[T_id]['length']
            S_name = '%s #%d' % (seqinfo[S_id]['name'], S_id)
            T_name = '%s #%d' % (seqinfo[T_id]['name'], T_id)

            vs = vs.union(set([S_name, T_name]))
            overlap_len = min(S_end, T_end) - max(S_start, T_start)
            # FIXME use overlap.overlap_direction
            if overlap_len > 0:
                if S_start < T_start:
                    es += [(S_name, T_name)]
                    ws += [overlap_len]
                elif S_start > T_start:
                    es += [(T_name, S_name)]
                    ws += [overlap_len]
                # if start is equal, edge goes from shorter read to longer read
                elif S_end < T_end:
                    es += [(S_name, T_name)]
                    ws += [overlap_len]
                elif S_end > T_end:
                    es += [(T_name, S_name)]
                    ws += [overlap_len]
                # if start and end is equal, reads are identical, ignore.

    G = overlap.OverlapGraph()
    G.iG.add_vertices(list(vs))
    es = [(G.iG.vs.find(name=u), G.iG.vs.find(name=v)) for u,v in es]
    G.iG.add_edges(es)
    G.iG.es['weight'] = ws
    return G
