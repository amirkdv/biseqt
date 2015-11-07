import re
import sys
import networkx as nx
import matplotlib.pyplot as plt
from termcolor import colored

from . import tuples
from . import align

def overlap_graph_by_alignment(tuplesdb, align_params=None, min_score=80):
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

def overlap_graph_by_tuple_extension(tuplesdb, align_params=None, window=20,
    drop_threshold=None):
    G = nx.DiGraph()
    seqinfo = tuplesdb.seqinfo()
    seqids = seqinfo.keys()
    sys.stdout.write('finding adjacent reads for sequence: ')
    for sid_idx in range(len(seqids)):
        sys.stdout.write('%d ' % seqids[sid_idx])
        sys.stdout.flush()
        for tid_idx in range(sid_idx + 1, len(seqids)):
            S_id, T_id = seqids[sid_idx], seqids[tid_idx]
            S_info, T_info = seqinfo[S_id], seqinfo[T_id]
            S_name = '%s_P%d' % (S_info['name'], S_info['start'])
            T_name = '%s_P%d' % (T_info['name'], T_info['start'])
            G.add_node(S_id, name=S_name)
            G.add_node(T_id, name=T_name)

            # do they overlap?
            S, T = tuplesdb.loadseq(S_id), tuplesdb.loadseq(T_id)
            F = tuples.OverlapFinder(S, T, align_params, tuplesdb=tuplesdb)
            exacts = F.exactly_matching_segments(S_id, T_id)
            if not exacts:
                continue
            segments = F.extend(exacts, window=window, drop_threshold=drop_threshold)
            if not segments:
                continue
            overlap = segments[0]
            S_len, T_len = F._S_len(overlap.tx.opseq), F._T_len(overlap.tx.opseq)
            if overlap.tx.idx_S == 0 and overlap.tx.idx_T == 0:
                if S_len < T_len:
                    G.add_edge(S_id, T_id, weight=overlap.tx.score)
                elif S_len > T_len:
                    G.add_edge(T_id, S_id, weight=overlap.tx.score)
            elif overlap.tx.idx_T == 0:
                G.add_edge(S_id, T_id, weight=overlap.tx.score)
            elif overlap.tx.idx_S == 0:
                G.add_edge(T_id, S_id, weight=overlap.tx.score)

    sys.stdout.write('\n')
    return G

def overlap_graph_by_known_order(tuplesdb):
    G = nx.DiGraph()
    seqinfo = tuplesdb.seqinfo()
    seqids = seqinfo.keys()
    for sid_idx in range(len(seqids)):
        for tid_idx in range(sid_idx + 1, len(seqids)):
            S_id, T_id = seqids[sid_idx], seqids[tid_idx]
            S_info, T_info = seqinfo[S_id], seqinfo[T_id]
            S_name = '%s_P%d' % (S_info['name'], S_info['start'])
            T_name = '%s_P%d' % (T_info['name'], T_info['start'])
            G.add_node(S_id, name=S_name)
            G.add_node(T_id, name=T_name)
            intersect_min = max(S_info['start'], T_info['start'])
            intersect_max = min(S_info['start'] + S_info['length'], T_info['start'] + T_info['length'])
            if intersect_min < intersect_max:
                overlap = intersect_max - intersect_min
                if S_info['start'] < T_info['start']:
                    G.add_edge(S_id, T_id, weight=overlap)
                elif S_info['start'] > T_info['start']:
                    G.add_edge(T_id, S_id, weight=overlap)
                # if start is equal, edge goes from shorter read to longer read
                elif S_info['start'] + S_info['length'] < T_info['start'] + T_info['length']:
                    G.add_edge(S_id, T_id, weight=overlap)
                elif S_info['start'] + S_info['length'] > T_info['start'] + T_info['length']:
                    G.add_edge(T_id, S_id, weight=overlap)

    return G

def save_graph(G, fname):
    nx.write_gml(G, fname)

# given path will be highlighted
def draw_graph(G, fname, figsize=None, path=[]):
    #pos = nx.circular_layout(G)
    #pos = nx.spring_layout(G)
    pos = nx.fruchterman_reingold_layout(G, k=10)
    if figsize is None:
        n = G.number_of_nodes()
        figsize = (n*2,n*2)
    plt.figure(figsize=figsize)
    # Vertices and their labels
    node_color = ['gray' if u in path else 'white' for u in G.nodes()]
    nx.draw_networkx_nodes(G, pos, node_size=8000, node_color=node_color)
    nx.draw_networkx_labels(G, pos, nx.get_node_attributes(G, 'name'), font_size=14)

    edge_data = G.edges(data=True)
    edge_in_path = lambda u,v: u in path and v in path and path[path.index(u) +1] == v
    edge_color = ['green' if edge_in_path(u,v) else 'black'  for u,v,_ in edge_data]
    edge_width = [4 if edge_in_path(u,v) else 0.1 for u,v,_ in edge_data]
    nx.draw_networkx_edges(G, pos, edge_color=edge_color, width=edge_width)
    if edge_data and 'weight' in edge_data[0][2]:
        nx.draw_networkx_edge_labels(G, pos, font_size=11,
            edge_labels={(f,t):'%.2f' % a['weight'] for f,t,a in edge_data})
    plt.xticks([])
    plt.yticks([])
    plt.savefig(fname, bbox_inches='tight')

def compare_graphs(G1, G2, f):
    E1, E2 = set(G1.edges()), set(G2.edges())
    diff = [('-', edge) for edge in E1 - E2] + [('+', edge) for edge in E2 - E1]
    N1 = nx.get_node_attributes(G1, 'name')
    N2 = nx.get_node_attributes(G2, 'name')
    for edge in sorted(diff, cmp=lambda x, y: cmp(x[1], y[1])):
        if edge[0] == '-':
            src, dst = N1[edge[1][0]], N1[edge[1][1]]
            line = '- [%s]--(%.2f)-->[%s]\n' % (src, G1.get_edge_data(*edge[1])['weight'], dst)
            f.write(colored(line, color='red'))
        else:
            src, dst = N2[edge[1][0]], N2[edge[1][1]]
            line = '+ [%s]--(%.2f)-->[%s]\n' % (src, G2.get_edge_data(*edge[1])['weight'], dst)
            f.write(colored(line, color='green'))

def layout(G):
    path = nx.algorithms.dag.dag_longest_path(G)
    V = dict(G.nodes(data=True))
    return [V[nid]['name'] for nid in path]

def draw_layout(G, fname):
    draw_graph(G, fname, path=nx.algorithms.dag.dag_longest_path(G))
