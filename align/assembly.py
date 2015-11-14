import re
import sys
import networkx as nx
import matplotlib.pyplot as plt
from termcolor import colored

from . import tuples
from . import pw

class OverlapGraph(object):
    """Wraps a networkx directed graph object with additional methods to build
    and process an overlap graph.

    Attributes:
        nxG (networkx.DiGraph): The graph object.

    Args:
        G (Optional[networkx.DiGraph]): The graph object to initialize with, no
        processing is done and if the object is ``None`` a new directed graph is
        instantiated.
    """
    def __init__(self, G=None):
        self.nxG = G if G else nx.DiGraph()
        assert(isinstance(self.nxG, nx.DiGraph))

    def V(self):
        """Returns a ``dict`` of vertex metadata keyed by vertex id."""
        return dict(self.nxG.nodes(data=True))

    def E(self):
        """Returns a ``dict`` of edge metadata keyed by the tuple of endpoint
        vertex ids."""
        return dict([((u,v), attrs) for u, v, attrs in self.nxG.edges(data=True)])

    def layout(self):
        """Finds the heaviest path of the directed graph and creates a new
        :class:`OverlapGraph` containing only this layout path.

        Returns:
            networkx.DiGraph: A linear subgraph (the heaviest path).

        Riases:
            AssertionError: if :attr:`nxG` is not acyclic.
        """
        assert(nx.algorithms.is_directed_acyclic_graph(self.nxG))
        V, E = self.V(), self.E()
        path = nx.algorithms.dag.dag_longest_path(self.nxG)
        L = OverlapGraph()
        for nid in V:
            L.nxG.add_node(nid, **V[nid])
        for nid_idx in range(1, len(path)):
            attrs = E[(path[nid_idx-1], path[nid_idx])]
            L.nxG.add_edge(path[nid_idx-1], path[nid_idx], attrs)

        return L

    def draw(self, fname, figsize=None, longest_path=False, pos=None,
        edge_colors=None, edge_width=None):
        """Draws the graph to a PDF file. If ``longest_path`` is truthy,
        the longest path in the graphs is found and and highlighted. If,
        however, the graph contains cycles the shortes cycle is highlighted in
        red."""
        V, E = self.V(), self.E()
        if longest_path:
            if nx.algorithms.is_directed_acyclic_graph(self.nxG):
                path = nx.algorithms.dag.dag_longest_path(self.nxG)
                edge_highlight = 'green'
            else:
                sys.stdout.write('Graph is not acyclic, shortest cycle is highlighted instead of the longest path\n')
                cycles = nx.algorithms.cycles.simple_cycles(self.nxG)
                path = sorted(cycles, key=lambda x: len(x))[0]
                path += [path[0]]
                edge_highlight = 'red'
        else:
            path = []

        if pos is None:
            pos = nx.circular_layout(self.nxG)
        if figsize is None:
            n = self.nxG.number_of_nodes()
            figsize = (n*2,n*2)
        plt.figure(figsize=figsize)

        # Vertices and their labels
        node_color = ['#b5ffb5' if u in path else '#ff9a9a' if path else 'w' for u in V]
        nx.draw_networkx_nodes(self.nxG, pos, node_size=8000, node_color=node_color)
        node_labels = nx.get_node_attributes(self.nxG, 'name')
        node_labels = {k: node_labels[k].replace(' ', '\n') for k in node_labels}
        nx.draw_networkx_labels(self.nxG, pos, node_labels, font_size=14)

        # Edges and their labels:
        edge_data = self.nxG.edges(data=True)
        mod = len(path) - 1
        in_path = lambda u,v: u in path and v in path and path.index(v) % mod == (path.index(u) + 1) % mod
        if edge_width is None:
            edge_width = [2 if in_path(u,v) else 0.7 for u,v,_ in edge_data]
        if edge_colors is None:
            edge_colors = [edge_highlight if in_path(u,v) else 'black'  for u,v,_ in edge_data]
        nx.draw_networkx_edges(self.nxG, pos, edge_color=edge_colors, width=edge_width)
        if edge_data and 'weight' in edge_data[0][2]:
            nx.draw_networkx_edge_labels(self.nxG, pos, font_size=11,
                edge_labels={(f,t):'%.2f' % a['weight'] for f,t,a in edge_data})

        plt.xticks([])
        plt.yticks([])
        plt.savefig(fname, bbox_inches='tight')

    def diff_text(self, OG, f):
        """Prints a diff-style comparison of our :attr:`nxG` against another
        given :class:`OverlapGraph` and writes the output to the given file
        handle. Missing edges are printed in red with a leading '-' and added
        edges are printed in green with a leading '+'.

        Args:
            OG (OverlapGraph): The "to" directed graph ("from" is us).
            f (file): File handle to write output to.
        """
        E1, E2 = set(self.nxG.edges()), set(OG.nxG.edges())
        missing, added = E1 - E2, E2 - E1
        f.write('G1 (%d edges) --> G2 (%d edges): %%%.2f lost, %%%.2f added\n' %
            (len(E1), len(E2), len(missing)*100.0/len(E1),
             len(added)*100.0/len(E1)))
        diff = [('-', edge) for edge in missing] + [('+', edge) for edge in added]
        N1 = nx.get_node_attributes(self.nxG, 'name')
        N2 = nx.get_node_attributes(OG.nxG, 'name')
        for edge in sorted(diff, cmp=lambda x, y: cmp(x[1], y[1])):
            if edge[0] == '-':
                src, dst = N1[edge[1][0]], N1[edge[1][1]]
                line = '- [%s]--(%.2f)-->[%s]\n' % (src,
                 self.nxG.get_edge_data(*edge[1])['weight'], dst)
                f.write(colored(line, color='red'))
            else:
                src, dst = N2[edge[1][0]], N2[edge[1][1]]
                line = '+ [%s]--(%.2f)-->[%s]\n' % \
                    (src, OG.nxG.get_edge_data(*edge[1])['weight'], dst)
                f.write(colored(line, color='green'))

    def diff_draw(self, OG, fname, figsize=None):
        """Draws the difference between our :attr:`nxG` against another
        given :class:`OverlapGraph`. Shared edges are in black, missing edges
        (from ours to ``OG``) are in red and added edges are in green.

        Args:
            OG (OverlapGraph): The "to" directed graph ("from" is us).
            fname (string): Path to which plot is saved, passed as is to
                :func:`draw`.
        """
        G = OverlapGraph()
        V1, E1 = self.V(), self.E()
        V2, E2 = OG.V(), OG.E()
        sE1, sE2 = set(self.nxG.edges()), set(OG.nxG.edges())
        for edge in sE1.union(sE2):
            G.nxG.add_node(edge[0], **V1[edge[0]])
            G.nxG.add_node(edge[1], **V1[edge[1]])
            kw = E1[edge] if edge in E1 else E2[edge]
            G.nxG.add_edge(edge[0], edge[1], **kw)

        V, E = G.V(), G.E()
        both, missing, added = sE1.intersection(sE2), sE1 - sE2, sE2 - sE1
        edge_colors = []
        for edge in G.nxG.edges():
            if edge in both:
                edge_colors += ['black']
            elif edge in missing:
                edge_colors += ['red']
            elif edge in added:
                edge_colors += ['green']

        pos = nx.fruchterman_reingold_layout(G.nxG, k=2.5)
        G.draw(fname, pos=pos, edge_colors=edge_colors, edge_width=2)

    def save(self, fname):
        """Saves the graph in GML format

        Args:
            fname (str): path to GML file.
        """
        nx.write_gml(self.nxG, fname)

    def break_cycles(self):
        """Breaks all directed cycles in the graph. No optimality guarantee
        is made, it seems to work in practice and time complexity is
        :math:`O(n\lg n)` where n is the number of cycles in the graph.
        """
        V, E = self.V(), self.E()
        # returns the edges of a cycle given the vertex list:
        fn_cycle_es = lambda x: set([(x[i-1], x[i]) for i in range(1, len(x))] + [(x[-1], x[0])])
        # returns the lightest edge in a list of edges
        fn_weakest = lambda C: sorted(C, key=lambda x: E[x]['weight'])[0]
        cycles = [c for c in nx.algorithms.cycles.simple_cycles(self.nxG)]
        if not cycles:
            return

        # process cycles in increasing order of weight:
        cycles = sorted(cycles, key=lambda c: sum([E[(c[i-1], c[i])]['weight'] for i in range(1, len(c))]))

        sys.stdout.write('Graph is not acyclic, breaking cycles: \n')
        # the intersection of all cycle edges seen so far:
        cands = fn_cycle_es(cycles[0])
        # the set of all edges to be eventually deleted:
        rm = set()
        for cycle in cycles:
            es = fn_cycle_es(cycle)
            if es.intersection(rm):
                # contains an edge that we have already decided to delete:
                continue
            # everytime the rolling intersection is about to become empty,
            # remove the lightest edge from the current set of candidates.
            new_cands = cands.intersection(es) if cands else es
            if not new_cands and cands:
                e = fn_weakest(cands)
                rm.update([e])
                sys.stdout.write('removed edge: %d --[%.2f]--> %d\n' % (e[0], E[e]['weight'], e[1]))
            cands = new_cands
        if cands:
            e = fn_weakest(cands)
            rm.update([e])
            sys.stdout.write('removed edge: %d --x[%.2f]x--> %d ' % (e[0], E[e]['weight'], e[1]))

        # actually remove the set of cycle-breaking edges:
        for e in rm:
            self.nxG.remove_edge(*e)
        sys.stdout.write('\n')
        if not nx.algorithms.dag.is_directed_acyclic_graph(self.nxG):
            sys.stdout.write('Err: failed to resolve all cycles of overlap graph.\n')

class OverlapBuilder(object):
    """Provided a :class:`align.tuples.Index` builds an overlap graph of all the
    sequences. All sequences must have already been indexed in the
    tuples database. For example::

        B = tuples.TuplesDB('path/to/file.db', alphabet=seq.Alphabet("ACGT"))
        I = tuples.Index(B, wordlen=10)
        C = align.AlignParams(
            ... # snip
        )
        G = assembly.OverlapBuilder(I, C).build()
        G.save(path)

    Args:
        index (tuples.Index): A tuples index that responds to
            :func:`align.tuples.Index.seeds`.
        align_params (pw.AlignParams): The alignment parameters for the
            rolling alignment.

    Keyword Args:
        drop_threshold (Optional[float]): What constitutes a drop in the
            score from one window to the next, default is 0. This means that
            if the overall score does not strictly increase (or the score of the
            new window is not positive) we drop the seed.
        window (Optional[int]): The size of the rolling window.
        max_succ_drops (Optional[int]): Maximum number of "drops" until the
            segment is dropped, default is 3.
        min_correct_seeds (Optional[int]): Minimum number of seeds extending
            successfully to a proper overlap alignment required to decide
            on an overlap edge. Default is 3.
    """

    def __init__(self, index, align_params, **kwargs):
        self.index, self.align_params = index, align_params
        self.window = kwargs.get('window', 20)
        self.drop_threshold = kwargs.get('drop_threshold', 0)
        self.max_succ_drops = kwargs.get('max_succ_drops', 3)
        self.min_correct_seeds = kwargs.get('min_correct_seeds')

    def build(self):
        """Builds a weighted, directed graph by using tuple methods. The process
        has 2 steps:

        * Find all seeds using :func:`align.tuples.Index.seeds`,
        * Extend all seeds to suffix-prefix segments using :func:`extend`.

        The resulting graph may not necessarily be acyclic. For further
        processing (e.g to find the layout) we need to ensure the overlap
        graph is acyclic. For this, see :func:`OverlapGraph.break_cycles`.

        Args:
            tuplesdb (tuples.TuplesDB): The tuples database.
            align_params (pw.AlignParams): Alignment parameters for overlap
                alignments.

        Returns:
            networkx.DiGraph
        """
        G = OverlapGraph()
        seqinfo = self.index.tuplesdb.seqinfo()
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
                G.nxG.add_node(S_id, name=S_name)
                G.nxG.add_node(T_id, name=T_name)

                # do they overlap?
                S = self.index.tuplesdb.loadseq(S_id)
                T = self.index.tuplesdb.loadseq(T_id)
                seeds = self.index.seeds(S_id, T_id)
                if not seeds:
                    continue
                overlap = self.extend(S, T, seeds)
                if not overlap:
                    continue
                #print set(['S->T' if x.tx.idx_S < x.tx.idx_T else 'T->S' for x in segments])
                #overlap.tx.pretty_print(tuplesdb.loadseq(S_id), tuplesdb.loadseq(T_id), sys.stdout)
                S_len = self._S_len(overlap.tx.opseq)
                T_len = self._T_len(overlap.tx.opseq)
                if abs(overlap.tx.idx_S - overlap.tx.idx_T) < self.window or \
                    abs(overlap.tx.idx_S + S_len - (overlap.tx.idx_T + T_len)) < self.window:
                    # end points are too close, ignore
                    continue
                if overlap.tx.idx_S == 0 and overlap.tx.idx_T == 0:
                    if S_len < T_len:
                        G.nxG.add_edge(S_id, T_id, weight=overlap.tx.score)
                    elif S_len > T_len:
                        G.nxG.add_edge(T_id, S_id, weight=overlap.tx.score)
                elif overlap.tx.idx_T == 0:
                    G.nxG.add_edge(S_id, T_id, weight=overlap.tx.score)
                elif overlap.tx.idx_S == 0:
                    G.nxG.add_edge(T_id, S_id, weight=overlap.tx.score)

        sys.stdout.write('\n')
        return G

    def _extend_fwd_once(self, S, T, segment, window):
        """Helper method for ``extend1d``."""
        S_len, T_len = self._S_len(segment.tx.opseq), self._T_len(segment.tx.opseq)
        align_problem_kw = {
            'S': S, 'T': T, 'params': self.align_params,
            'align_type': pw.GLOBAL,
        }
        align_problem_kw.update({
            'S_min_idx': segment.tx.idx_S + S_len,
            'T_min_idx': segment.tx.idx_T + T_len,
        })
        align_problem_kw.update({
            'S_max_idx': align_problem_kw['S_min_idx'] + window,
            'T_max_idx': align_problem_kw['T_min_idx'] + window
        })

        with pw.AlignProblem(**align_problem_kw) as P:
            score = P.solve()
            assert(score is not None)
            transcript = P.traceback()
            if transcript is None:
                return None

        tx = pw.Transcript(
            idx_S=segment.tx.idx_S,
            idx_T=segment.tx.idx_T,
            score=segment.tx.score + transcript.score,
            opseq=segment.tx.opseq + transcript.opseq
        )
        return tuples.Segment(S_id=segment.S_id, T_id=segment.T_id, tx=tx)

    def _extend_bwd_once(self, S, T, segment, window):
        """Helper method for ``extend1d``."""
        align_problem_kw = {
            'S': S, 'T': T, 'params': self.align_params,
            'align_type': pw.GLOBAL,
        }
        align_problem_kw.update({
            'S_max_idx': segment.tx.idx_S,
            'T_max_idx': segment.tx.idx_T,
        })
        align_problem_kw.update({
            'S_min_idx': align_problem_kw['S_max_idx'] - window,
            'T_min_idx': align_problem_kw['T_max_idx'] - window
        })

        with pw.AlignProblem(**align_problem_kw) as P:
            score = P.solve()
            assert(score is not None)
            transcript = P.traceback()
            if transcript is None:
                return None

        tx = pw.Transcript(
            idx_S=transcript.idx_S,
            idx_T=transcript.idx_T,
            score=transcript.score + segment.tx.score,
            opseq=transcript.opseq + segment.tx.opseq
        )
        return tuples.Segment(S_id=segment.S_id, T_id=segment.T_id, tx=tx)

    def _S_len(self, opseq):
        return sum([opseq.count(op) for op in 'DMS'])

    def _T_len(self, opseq):
        return sum([opseq.count(op) for op in 'IMS'])

    def _extend1d(self, S, T, segment, backwards=False):
        """Helper method for ``extend()``"""
        cur_seg = segment
        score_history = [segment.tx.score]
        while True:
            if backwards:
                w = min(self.window, min(cur_seg.tx.idx_S, cur_seg.tx.idx_T))
            else:
                S_wiggle = S.length - (cur_seg.tx.idx_S + self._S_len(cur_seg.tx.opseq))
                T_wiggle = T.length - (cur_seg.tx.idx_T + self._T_len(cur_seg.tx.opseq))
                w = min(self.window, min(S_wiggle, T_wiggle))

            if w == 0:
                # hit the end:
                return cur_seg

            if backwards:
                seg = self._extend_bwd_once(S, T, cur_seg, w)
            else:
                seg = self._extend_fwd_once(S, T, cur_seg, w)

            if seg is None:
                # no non-empty alignment found.
                break

            score_history += [seg.tx.score - segment.tx.score]
            if len(score_history) > self.max_succ_drops:
                score_history = score_history[-self.max_succ_drops:]

            if len(score_history) == self.max_succ_drops and \
                all([x <= self.drop_threshold for x in score_history]):
                break

            cur_seg = seg

        return None

    def extend(self, S, T, segments):
        """Given two sequences and a number of matching segments finds all
        extended, potentially gap-containing segments by repeatedly aligning a
        rolling frame along the two sequences and dropping a segment once it
        is observed to not be a high-scoring overlap alignment.

        Args:
            S (seq.Sequence): The "from" sequence.
            T (seq.Sequence): The "to" sequence.
            segments (List[tuples.Segment]): The starting segments. If called
                from :func:`build`, these are seeds but no assumption is made.
        """
        res = []
        for segment in segments:
            fwd = self._extend1d(S, T, segment)
            bwd = self._extend1d(S, T, segment, backwards=True)
            if fwd and bwd and min(fwd.tx.score, bwd.tx.score) > self.drop_threshold:
                assert(bwd.tx.idx_S == 0 or bwd.tx.idx_T == 0)
                opseq = bwd.tx.opseq[:-len(segment.tx.opseq)] + fwd.tx.opseq
                score = self.align_params.score(
                    S, T, opseq,
                    S_min_idx=bwd.tx.idx_S,
                    T_min_idx=bwd.tx.idx_T
                )
                tx = pw.Transcript(
                    idx_S=bwd.tx.idx_S,
                    idx_T=bwd.tx.idx_T,
                    score=score,
                    opseq=opseq
                )
                res += [tuples.Segment(S_id=segment.S_id, T_id=segment.T_id, tx=tx)]
                if len(res) >= self.min_correct_seeds:
                    break
        if not res:
            return None
        return sorted(res, key=lambda s: s.tx.score, reverse=True)[0]
