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
    drop_threshold=0, max_succ_drops=3):
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
            S_min_idx, S_max_idx = S_info['start'], S_info['start'] + S_info['length']
            T_min_idx, T_max_idx = T_info['start'], T_info['start'] + T_info['length']
            S_name = '%s %d-%d #%d' % (S_info['name'], S_min_idx, S_max_idx, S_id)
            T_name = '%s %d-%d #%d' % (T_info['name'], T_min_idx, T_max_idx, T_id)
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
            #print set(['S->T' if x.tx.idx_S < x.tx.idx_T else 'T->S' for x in segments])
            #overlap.tx.pretty_print(tuplesdb.loadseq(S_id), tuplesdb.loadseq(T_id), sys.stdout)
            S_len, T_len = F._S_len(overlap.tx.opseq), F._T_len(overlap.tx.opseq)
            if abs(overlap.tx.idx_S - overlap.tx.idx_T) < window or \
                abs(overlap.tx.idx_S + S_len - (overlap.tx.idx_T + T_len)) < window:
                continue
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

    # break cycles
    V, E = _dict_VE_from_graph(G)
    fn_cycle_es = lambda x: set([(x[i-1], x[i]) for i in range(1, len(x))] + [(x[-1], x[0])])
    if not nx.algorithms.dag.is_directed_acyclic_graph(G):
        sys.stdout.write('Graph is not acyclic, breaking cycles: ')
        cycles = nx.algorithms.cycles.simple_cycles(G)
        candidates = fn_cycle_es(cycles.next())
        fn_weakest = lambda cands: sorted(candidates, key=lambda x: E[x]['weight'])[0]
        for c in cycles:
            new_cands = candidates.intersection(fn_cycle_es(c))
            if not new_cands:
                e = fn_weakest(candidates)
                sys.stdout.write('%d --x--> %d ' % e)
                G.remove_edge(*e)
            candidates = new_cands
        if candidates:
            e = fn_weakest(candidates)
            sys.stdout.write('%d --x--> %d ' % e)
            G.remove_edge(*e)
        sys.stdout.write('\n')
    assert(nx.algorithms.dag.is_directed_acyclic_graph(G))
    return G

def overlap_graph_by_known_order(tuplesdb):
    G = nx.DiGraph()
    seqinfo = tuplesdb.seqinfo()
    seqids = seqinfo.keys()
    for sid_idx in range(len(seqids)):
        for tid_idx in range(sid_idx + 1, len(seqids)):
            S_id, T_id = seqids[sid_idx], seqids[tid_idx]
            S_info, T_info = seqinfo[S_id], seqinfo[T_id]
            S_min_idx, S_max_idx = S_info['start'], S_info['start'] + S_info['length']
            T_min_idx, T_max_idx = T_info['start'], T_info['start'] + T_info['length']
            S_name = '%s %d-%d #%d' % (S_info['name'], S_min_idx, S_max_idx, S_id)
            T_name = '%s %d-%d #%d' % (T_info['name'], T_min_idx, T_max_idx, T_id)
            G.add_node(S_id, name=S_name)
            G.add_node(T_id, name=T_name)
            overlap = min(S_max_idx, T_max_idx) - max(S_min_idx, T_min_idx)
            if overlap > 0:
                if S_min_idx < T_min_idx:
                    G.add_edge(S_id, T_id, weight=overlap)
                elif S_min_idx > T_min_idx:
                    G.add_edge(T_id, S_id, weight=overlap)
                # if start is equal, edge goes from shorter read to longer read
                elif S_max_idx < T_max_idx:
                    G.add_edge(S_id, T_id, weight=overlap)
                elif S_max_idx > T_max_idx:
                    G.add_edge(T_id, S_id, weight=overlap)
                # if start and end is equal, reads are identical, ignore.

    return G

def save_graph(G, fname):
    nx.write_gml(G, fname)

def _dict_VE_from_graph(G):
    V = dict(G.nodes(data=True))
    E = dict([((u,v), attrs) for u, v, attrs in G.edges(data=True)])
    return V, E

def draw_digraph(G, fname, figsize=None, longest_path=False, pos=None,
    edge_colors=None, edge_width=None):
    V, E = _dict_VE_from_graph(G)
    if longest_path:
        if nx.algorithms.is_directed_acyclic_graph(G):
            path = nx.algorithms.dag.dag_longest_path(G)
            edge_highlight = 'green'
        else:
            sys.stdout.write('Graph is not acyclic, longest cycle is highlighted instead of the longest path\n')
            cycles = nx.algorithms.cycles.simple_cycles(G)
            path = sorted(cycles, key=lambda x: len(x), reverse=True)[0]
            path += [path[0]]
            edge_highlight = 'red'
    else:
        path = []

    if pos is None:
        pos = nx.circular_layout(G)
    if figsize is None:
        n = G.number_of_nodes()
        figsize = (n*2,n*2)
    plt.figure(figsize=figsize)

    # Vertices and their labels
    node_color = ['gray' if u in path else 'white' for u in V]
    nx.draw_networkx_nodes(G, pos, node_size=8000, node_color=node_color)
    node_labels = nx.get_node_attributes(G, 'name')
    node_labels = {k: node_labels[k].replace(' ', '\n') for k in node_labels}
    nx.draw_networkx_labels(G, pos, node_labels, font_size=14)

    # Edges and their labels:
    edge_data = G.edges(data=True)
    mod = len(path) - 1
    in_path = lambda u,v: u in path and v in path and path.index(v) == (path.index(u) + 1) % mod
    if edge_width is None:
        edge_width = [2 if in_path(u,v) else 0.7 for u,v,_ in edge_data]
    if edge_colors is None:
        edge_colors = [edge_highlight if in_path(u,v) else 'black'  for u,v,_ in edge_data]
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_width)
    if edge_data and 'weight' in edge_data[0][2]:
        nx.draw_networkx_edge_labels(G, pos, font_size=11,
            edge_labels={(f,t):'%.2f' % a['weight'] for f,t,a in edge_data})

    plt.xticks([])
    plt.yticks([])
    plt.savefig(fname, bbox_inches='tight')

def layout_graph(G):
    assert(nx.algorithms.is_directed_acyclic_graph(G))
    V, E = _dict_VE_from_graph(G)
    path = nx.algorithms.dag.dag_longest_path(G)
    L = nx.DiGraph()
    for nid in V:
        L.add_node(nid, **V[nid])
    for nid_idx in range(1, len(path)):
        attrs = E[(path[nid_idx-1], path[nid_idx])]
        L.add_edge(path[nid_idx-1], path[nid_idx], attrs)

    return L

def diff_graph(G1, G2, fname, figsize=None):
    G = nx.DiGraph()
    V1, E1 = _dict_VE_from_graph(G1)
    V2, E2 = _dict_VE_from_graph(G2)
    for edge in set(G1.edges()).union(set(G2.edges())):
        G.add_node(edge[0], **V1[edge[0]])
        G.add_node(edge[1], **V1[edge[1]])
        kw = E1[edge] if edge in E1 else E2[edge]
        G.add_edge(edge[0], edge[1], **kw)

    V, E = _dict_VE_from_graph(G)
    sE1, sE2 = set(G1.edges()), set(G2.edges())
    both, missing, added = sE1.intersection(sE2), sE1 - sE2, sE2 - sE1
    edge_colors = []
    for edge in G.edges():
        if edge in both:
            edge_colors += ['black']
        elif edge in missing:
            edge_colors += ['red']
        elif edge in added:
            edge_colors += ['green']

    pos = nx.fruchterman_reingold_layout(G, k=2.5)
    draw_digraph(G, fname, pos=pos, edge_colors=edge_colors, edge_width=2)

def compare_graphs(G1, G2, f):
    E1, E2 = set(G1.edges()), set(G2.edges())
    missing, added = E1 - E2, E2 - E1
    f.write('G1 (%d edges) --> G2 (%d edges): %%%.2f lost, %%%.2f added\n' %
        (len(E1), len(E2), len(missing)*100.0/len(E1),
         len(added)*100.0/len(E1)))
    diff = [('-', edge) for edge in missing] + [('+', edge) for edge in added]
    N1 = nx.get_node_attributes(G1, 'name')
    N2 = nx.get_node_attributes(G2, 'name')
    for edge in sorted(diff, cmp=lambda x, y: cmp(x[1], y[1])):
        if edge[0] == '-':
            src, dst = N1[edge[1][0]], N1[edge[1][1]]
            line = '- [%s]--(%.2f)-->[%s]\n' % (src,
             G1.get_edge_data(*edge[1])['weight'], dst)
            f.write(colored(line, color='red'))
        else:
            src, dst = N2[edge[1][0]], N2[edge[1][1]]
            line = '+ [%s]--(%.2f)-->[%s]\n' % (src,
             G2.get_edge_data(*edge[1])['weight'], dst)
            f.write(colored(line, color='green'))
