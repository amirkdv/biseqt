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
            G.add_node(seqids[idx_of_S])
            G.add_node(seqids[idx_of_T])
            S = tuplesdb.loadseq(seqids[idx_of_S])
            T = tuplesdb.loadseq(seqids[idx_of_T])
            with align.AlignProblem(S=S, T=T, params=align_params,
                align_type=align.ALIGN_OVERLAP) as P:
                score = P.solve()
                if score >= min_score:
                    transcript = P.traceback()
                    if transcript.idx_T == 0:
                        G.add_edge(seqids[idx_of_S], seqids[idx_of_T], score=score)
                    if transcript.idx_S == 0:
                        G.add_edge(seqids[idx_of_T], seqids[idx_of_S], score=score)
    return G

def overlap_graph_by_tuple_extension(tuplesdb, align_params, drop_threshold):
    G = nx.DiGraph()
    seqids = tuplesdb.seqids()
    for idx_of_S in range(len(seqids)):
        for idx_of_T in range(idx_of_S + 1, len(seqids)):
            print idx_of_S+1, idx_of_T+1
            G.add_node(seqids[idx_of_S])
            G.add_node(seqids[idx_of_T])
            S = tuplesdb.loadseq(seqids[idx_of_S])
            T = tuplesdb.loadseq(seqids[idx_of_T])
            F = tuples.OverlapFinder(S, T, align_params, tuplesdb=tuplesdb)
            exacts = F.exactly_matching_segments(
                seqids[idx_of_S], seqids[idx_of_T]
            )
            if not exacts:
                continue
            segments = F.extend(exacts, drop_threshold)
            if not segments:
                continue
            overlap = segments[0]
            if overlap.tx.idx_T == 0:
                G.add_edge(seqids[idx_of_S], seqids[idx_of_T], score=overlap.tx.score)
            if overlap.tx.idx_S == 0:
                G.add_edge(seqids[idx_of_T], seqids[idx_of_S], score=overlap.tx.score)
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
                G.add_edge(sid, tid, score=overlap)

    return G


def save_overlap_graph(G, path, figsize=None):
    pos = nx.circular_layout(G)
    if figsize is None:
        n = G.number_of_nodes()
        figsize = (n*5,n*5)
    plt.figure(figsize=figsize)
    nx.draw_networkx_nodes(G, pos, node_size=10000, node_color='w')
    nx.draw_networkx_labels(G, pos, font_size=100)
    nx.draw_networkx_edges(G, pos, width=5)
    edge_data = G.edges(data=True)
    if edge_data and 'score' in edge_data[0][2]:
        nx.draw_networkx_edge_labels(G, pos, font_size=26, edge_labels={(f,t):'%.2f' % a['score'] for f,t,a in edge_data})
    plt.xticks([])
    plt.yticks([])
    plt.savefig(path, bbox_inches='tight')
