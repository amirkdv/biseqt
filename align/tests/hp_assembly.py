#!/usr/bin/env python
import sys
import os
import igraph
from .. import pw, words, seq, overlap, homopolymeric, overlap, mapping, ProgressIndicator
from ..mapping import Mapping

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
od_params = overlap.OverlapDiscoveryParams(
    seed_ext_params=overlap.SeedExtensionParams(
        window=params['window'],
        min_overlap_score=params['min_overlap_score'],
        max_new_mins=params['max_new_mins'],
        align_params=C
    ),
    shift_rolling_sum_width=params['shift_rolling_sum_width'],
    hp_condenser=Tr,
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
        #hp_prob=params['hp_prob']
    )

def map_reads(db, reference, reads):
    B = words.SeqDB(db, alphabet=A)
    M = assembly.BwaMappingAssembler(B)
    M.initialize(reference, reads)
    M.save(M.mappings(num_reads=os.environ.get('NUM_READS', -1)))

def true_overlaps(true_path):
    G = igraph.read(true_path)
    return [set([db_id(G, u), db_id(G, v)]) for u, v in G.get_edgelist()]

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
        params['shift_rolling_sum_width'],
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
    true_overlaps = [set([db_id(G, u), db_id(G, v)]) for u, v in G.get_edgelist()]
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

def build_true_overlap_graph(db, mappings):
    B = seq.SeqDB(db, alphabet=A)
    seqinfo = B.seqinfo()
    indicator = ProgressIndicator(
        'Building true overlap graph from mappings at %s' % mappings,
        len(seqinfo)*(len(seqinfo) - 1)/2.0,
    )
    indicator.start()
    seqids = seqinfo.keys()
    with open(mappings) as f:
        mappings = eval(f.read())
    vs = set()
    es, ws = [], []
    for sid_idx in range(len(seqids)):
        for tid_idx in range(sid_idx + 1, len(seqids)):
            indicator.progress()
            S_id, T_id = seqids[sid_idx], seqids[tid_idx]
            S_name, T_name = seqinfo[S_id]['name'], seqinfo[T_id]['name']
            S_start, S_end = mappings[S_name].ref_from, mappings[S_name].ref_to
            T_start, T_end = mappings[T_name].ref_from, mappings[T_name].ref_to

            if mappings[S_name].strand == '-':
                S_start, S_end = S_end, S_start
            if mappings[T_name].strand == '-':
                T_start, T_end = T_end, T_start

            S_name = '%s #%d' % (S_name, S_id)
            T_name = '%s #%d' % (T_name, T_id)

            vs = vs.union(set([S_name, T_name]))
            overlap_len = min(S_end, T_end) - max(S_start, T_start)
            if overlap_len > 0:
                if S_start < T_start or (S_start == T_start and S_end < T_end):
                    es += [(S_name, T_name)]
                else:
                    es += [(T_name, S_name)]
                ws += [overlap_len]

    indicator.finish()
    G = overlap.OverlapGraph()
    G.iG.add_vertices(list(vs))
    es = [(G.iG.vs.find(name=u), G.iG.vs.find(name=v)) for u,v in es]
    G.iG.add_edges(es)
    G.iG.es['weight'] = ws
    return G
