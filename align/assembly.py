from . import tuples
from . import align
import re
import networkx as nx
import matplotlib.pyplot as plt

def overlap_graph_by_alignment(tuplesdb, align_params, min_score=80):
    G = nx.DiGraph()
    seqids = tuplesdb.seqids()
    for idx_of_S in range(len(seqids)):
        for idx_of_T in range(idx_of_S + 1, len(seqids)):
            S = tuplesdb.loadseq(seqids[idx_of_S])
            T = tuplesdb.loadseq(seqids[idx_of_T])
            P = align.AlignProblem(
                S=S, T=T,
                params=align_params,
                align_type=align.ALIGN_OVERLAP
            )
            score = P.solve()
            if score >= min_score:
                transcript = P.traceback()
                if transcript.idx_T == 0:
                    G.add_edge(seqids[idx_of_S], seqids[idx_of_T], score=score)
                if transcript.idx_S == 0:
                    G.add_edge(seqids[idx_of_T], seqids[idx_of_S], score=score)
    return G

def overlap_graph_by_tuple_extension(tuplesdb, align_params, max_decr=3, decr_def=0):
    G = nx.DiGraph()
    seqids = tuplesdb.seqids()
    for idx_of_S in range(len(seqids)):
        for idx_of_T in range(idx_of_S + 1, len(seqids)):
            S = tuplesdb.loadseq(seqids[idx_of_S])
            T = tuplesdb.loadseq(seqids[idx_of_T])
            F = tuples.OverlapFinder(S, T, align_params)
            exacts = tuplesdb.exactly_matching_segments(
                seqids[idx_of_S], seqids[idx_of_T]
            )
            if not exacts:
                continue
            segments = F.extend(exacts, max_decr, decr_def)
            if not segments:
                continue
            overlap, score = segments[0]
            if overlap.idx_T == 0:
                G.add_edge(seqids[idx_of_S], seqids[idx_of_T], score=score)
            if overlap.idx_S == 0:
                G.add_edge(seqids[idx_of_T], seqids[idx_of_S], score=score)
    return G

def overlap_graph_by_known_order(tuplesdb):
    G = nx.DiGraph()
    seqs = tuplesdb.seqids(info=True)
    for sid in seqs:
        for tid in seqs:
            if sid == tid:
                continue
            overlap = seqs[sid]['start'] + seqs[sid]['length'] - seqs[tid]['start']
            if seqs[tid]['start'] >= seqs[sid]['start'] and overlap > 0:
                G.add_edge(sid, tid)

    return G


def save_overlap_graph(G, path, figsize=(50,50)):
    pos = nx.circular_layout(G)
    plt.figure(figsize=figsize)
    nx.draw_networkx_nodes(G, pos, node_size=2000, node_color='w')
    nx.draw_networkx_labels(G, pos, font_size=30)
    nx.draw_networkx_edges(G, pos, width=5)
    edge_data = G.edges(data=True)
    if edge_data and 'score' in edge_data[0]:
        nx.draw_networkx_edge_labels(G, pos, font_size=26, edge_labels={(f,t):w['score'] for f,t,w in G.edges(data=True)})
    plt.xticks([])
    plt.yticks([])
    plt.savefig(path)
