import re
import sys
import igraph
from termcolor import colored
import time
from . import tuples, pw, homopolymeric, ProgressIndicator


class OverlapGraph(object):
    """Wraps an :class:`igraph.Graph` object with additional methods to build
    and process an overlap graph.

    Attributes:
        iG (igraph.Graph): The graph object.

    Args:
        G (Optional[igraph.Graph]): The graph object to initialize with; no
            processing is done and if the object is ``None`` a new directed
            graph is instantiated.
    """
    def __init__(self, G=None):
        self.iG = G if G else igraph.Graph(directed=True)
        assert(isinstance(self.iG, igraph.Graph))
        self.v_highlight = '#b5ffb5'
        self.e_highlight = '#00b400'

    def _endpoint_names(self, eid):
        """Internal helper: igraph is not too happy when we munge vertex IDs
        since it keeps renumbering them according to its memory allocation
        scheme. Instead convert everything to "name"s which are the original
        sequence IDs.

        See: https://lists.nongnu.org/archive/html/igraph-help/2010-03/msg00078.html
        """
        if isinstance(eid, igraph.Edge):
            eid = eid.index
        uid, vid = self.iG.es[eid].tuple
        return self.iG.vs[uid]['name'], self.iG.vs[vid]['name']

    def _vids_to_names(self, vids):
        return [self.iG.vs[vid]['name'] for vid in vids]

    def _edge_weight(self, u, v):
        """Internal helper; given to vertex names returns the weight of the
        edge connecting them.
        """
        return self.iG.es['weight'][self.iG.get_eid(u, v)]

    def eid_to_str(self, eid, maxlen=50):
        """Prepares an edge for pretty printing. Truncates and paths the end
        point labels (``name`` is used as label) to ensure they both have
        length ``maxlen``.
        """
        u, v = self._endpoint_names(eid)
        u, v = u[:maxlen].ljust(maxlen), v[:maxlen].rjust(maxlen)
        w = self.iG.es[eid]['weight']
        w = ('+--[%.2f]-->' % w).ljust(20)
        return '%s %s %s\n' % (u, w, v)

    def break_cycles(self, method='ip'):
        """Removes a
        `feedback arc set <https://en.wikipedia.org/wiki/Feedback_arc_set>`__
        from the graph. Depending on the ``method`` the result may not be
        optimal.

        Keyword Args:
            method (str): The FAS discovery algorithm; passed to
                :func:`igraph.Graph.feedback_arc_set`. Default uses an
                integer programming formulation which is guaranteed to be
                optimal but is slow on large graphs. The alternative is
                ``eades`` which uses a suboptimal `heuristic
                <http://www.sciencedirect.com/science/article/pii/002001909390079O>`__.
        """
        if self.iG.is_dag():
            return
        rm = self.iG.feedback_arc_set(
            weights=self.iG.es['weight'], method=method
        )
        for e in rm:
            sys.stderr.write('removed edge: %s' % self.eid_to_str(e))
        self.iG.delete_edges(rm)

    def longest_path(self, exclude=[]):
        """Finds the longest path (i.e heaviest path) of the graph, excluding
        vertices with names specified in ``exclude``. This, naturally requires
        that the graph is acyclic. Assuming the graph is a DAG, we can find the
        longest path in two steps:
            - Find a topological ordering of the graph in :math:`O(|V|+|E|)`
              time,
            - Find a heaviest path using the sorting in :math:`O(|V|)` time.

        Keyword Arguments:
            exclude (Optional[List[str]]): A list of vertex names to be
                excluded from the graph when finding the longest path. This is
                only of use to :func:`all_longest_paths`.

        Returns:
            list[str]: A list of vertex names in order of appearance in the
                longest path.
        """
        sorting = self._vids_to_names(self.iG.topological_sorting())
        sorting = [v for v in sorting if v not in exclude]
        # longest paths ending at each vertex keyed by vertex. Each entry is a
        # tuple of (<weight, from>) where `from` is any of the predecessors
        # giving the maximum weight.
        longest_paths = {}
        for v in sorting:
            if v in exclude:
                continue
            incoming = self._vids_to_names(self.iG.predecessors(v))
            incoming = [x for x in incoming if x not in exclude]
            if not incoming:
                longest_paths[v] = (0, None)
            else:
                w = lambda x: longest_paths[x][0] + self._edge_weight(x, v)
                cands = [(w(u), u) for u in incoming]
                longest_paths[v] = sorted(
                    cands, key=lambda x: x[0], reverse=True
                )[0]

        if not longest_paths:
            return []

        # Find the terminal vertex of the longest path:
        end = sorted(
            longest_paths.items(), key=lambda x: x[1][0], reverse=True
        )[0][0]
        path = []
        # Trace back the entire path:
        while end and longest_paths:
            path = [end] + path
            end = longest_paths.pop(end)[1]

        # Don't report trivial paths:
        return path if len(path) > 1 else []

    def all_longest_paths(self):
        """Repeatedly finds the longest path in the graph while excluding
        vertices that are already included in a path. See :func:`longest_path`.

        Returns:
            List[List[str]]: A list of paths, each a list of vertex names in
                order of appearance in the path.
        """
        paths = []
        exclude = []
        while True:
            path = self.longest_path(exclude=exclude)
            if not path:
                break
            paths += [path]
            exclude += path
        return paths

    def layout(self, full=False):
        """Finds the heaviest path of the directed graph and creates a new
        :class:`OverlapGraph` containing only this layout path.

        Keyword Args:
            full (bool): If truthy, an effort is made to add other paths to
                cover all vertices of the graph.

        Returns:
            assembly.OverlapGraph: A linear subgraph (the heaviest path).

        Raises:
            AssertionError: If the graph is not acyclic.
        """
        assert(self.iG.is_dag())
        if full:
            paths = self.all_longest_paths()
        else:
            paths = [self.longest_path()]
        eids = []
        for path in paths:
            for idx in range(1, len(path)):
                eids += [self.iG.get_eid(path[idx-1], path[idx])]

        return OverlapGraph(self.iG.subgraph_edges(eids))

    # the paths are names not ids
    def draw(self, fname, **kw):
        """Draws the graph and potentially highlights provided paths.

        Keyword Arguments:
            highlight_paths ([List[List[str]]]): A list of paths to be
                highlighted. All edges of the path and the starting vertex
                are highlighted green.
            edge_color ([List|str]): Passed to :func:`igraph.Graph.plot`.
                Default is all black unless paths to be highlighted are
                specified. If provided, overrides path highlighting.
            vertex_color ([List|str]): Passed to :func:`igraph.Graph.plot`.
                Default is all white unless paths to be highlighted are
                specified in which case starting vertices are green.
            edge_width ([List|float]): Passed to :func:`igraph.Graph.plot`.
                Default is 10 for edges in highlighted path and 1 otherwise.
            edge_arrow_widge ([List|float]): Passed to
                :func:`igraph.Graph.plot`. Default is 3 for highlighted edges
                and 1 otherwise.
            edge_curvred (float): Passed to :func:`igraph.Graph.plot`. Default
                is 0.1.
        """
        highlight_paths = kw.get('highlight_paths', [])

        def e_in_path(eid):
            u, v = self._endpoint_names(eid)
            return any([
                u in p and v in p and
                p.index(u) == p.index(v) - 1 for p in highlight_paths
            ])

        v_start_path = lambda v: any([p[0] == v for p in highlight_paths])

        # Sugiyama works on non-DAG graphs as well
        n = len(self.iG.vs)
        layout_kw = {'maxiter': n * 20, 'weights': None}
        if 'weight' in self.iG.es.attributes():
            layout_kw['weights'] = 'weight'
        plot_kw = {
            'layout': self.iG.layout_sugiyama(**layout_kw),
            'bbox': (n*150, n*150),
            'vertex_size': 100,
            'vertex_label': [x.replace(' ', '\n') for x in self.iG.vs['name']],
            'vertex_label_size': 18,
            'vertex_color': kw.get('vertex_color', [self.v_highlight if v_start_path(v) else 'white' for v in self.iG.vs['name']]),
            'edge_width': kw.get('edge_width', [10 if e_in_path(e) else 1 for e in self.iG.es]),
            'edge_arrow_width': kw.get('edge_arrow_width', [3 if e_in_path(e) else 1 for e in self.iG.es]),
            'edge_color': kw.get('edge_color', [self.e_highlight if e_in_path(e) else 'black' for e in self.iG.es]),
            'edge_curved': kw.get('edge_curved', 0.1),
            'margin': 200,
        }
        igraph.plot(self.iG, fname, **plot_kw)

    def diff_text(self, OG, f=sys.stdout, summary_only=True):
        """Prints a diff-style comparison of our :attr:`iG` against another
        given :class:`OverlapGraph` and writes the output to the given file
        handle. Missing edges are printed in red with a leading '-' and added
        edges are printed in green with a leading '+'.

        Args:
            OG (OverlapGraph): The "to" directed graph ("from" is us).
            f (file): File handle to write output to.
        """
        sE1 = set([self._endpoint_names(e) for e in self.iG.es])
        sE2 = set([OG._endpoint_names(e) for e in OG.iG.es])

        def _edge_str(endpoints):
            if endpoints in sE1:
                return self.eid_to_str(self.iG.get_eid(*endpoints))
            elif endpoints in sE2:
                return OG.eid_to_str(OG.iG.get_eid(*endpoints))
            else:
                raise RuntimeError("This should not have happened")
        missing, added, both = sE1 - sE2, sE2 - sE1, sE1.intersection(sE2)
        missing_pctg = len(missing)*100.0/len(sE1)
        added_pctg = len(added)*100.0/len(sE1)
        f.write(
            'G1 (%d edges) --> G2 (%d edges): %%%.2f lost, %%%.2f added\n' %
            (len(sE1), len(sE2), missing_pctg, added_pctg)
        )
        if summary_only:
            return
        diff = [('-', edge) for edge in missing] + \
            [('+', edge) for edge in added] + [(None, edge) for edge in both]
        for edge in sorted(diff, cmp=lambda x, y: cmp(x[1], y[1])):
            color = None
            prefix = ' ' if edge[0] is None else edge[0]
            line = '%s %s' % (prefix, _edge_str(edge[1]))
            if edge[0] == '-':
                color = 'red'
            elif edge[0] == '+':
                color = 'green'

            if color:
                f.write(colored(line, color=color))
            else:
                f.write(line)

    def diff_draw(self, OG, fname, figsize=None):
        """Draws the difference between our :attr:`iG` against another
        given :class:`OverlapGraph`. Shared edges are in black, missing edges
        (from ours to ``OG``) are in red and added edges are in green.

        Args:
            OG (OverlapGraph): The "to" directed graph ("from" is us).
            fname (string): Path to which plot is saved, passed as is to
                :func:`draw`.
        """
        e_to_names = lambda G, e: (G.vs[e[0]]['name'], G.vs[e[1]]['name'])
        sE1 = set([self._endpoint_names(e) for e in self.iG.es])
        sE2 = set([OG._endpoint_names(e) for e in OG.iG.es])
        G = OverlapGraph()
        G.iG.add_vertices(list(
            set(self.iG.vs['name']).union(set(OG.iG.vs['name']))
        ))
        G.iG.add_edges(list(sE1.union(sE2)))
        both, missing, added = sE1.intersection(sE2), sE1 - sE2, sE2 - sE1
        edge_color = []
        for e in G.iG.es:
            e = G._endpoint_names(e)
            if e in both:
                edge_color += ['black']
            elif e in missing:
                edge_color += ['red']
            elif e in added:
                edge_color += ['green']

        vertex_color = ['white' if v.degree(mode=igraph.IN) else self.v_highlight for v in G.iG.vs]

        G.draw(fname, edge_color=edge_color, vertex_color=vertex_color,
               edge_width=5, edge_arrow_width=3, edge_curved=0.01)

    def save(self, fname):
        """Saves the graph in GML format

        Args:
            fname (str): path to GML file.
        """
        self.iG.write_gml(fname)


class OverlapBuilder(object):
    """Provided a :class:`align.tuples.Index` builds an overlap graph of all
    the sequences. All sequences must have already been indexed in the
    tuples database. For example::

        B = tuples.TuplesDB('path/to/file.db', alphabet=seq.Alphabet("ACGT"))
        I = tuples.Index(B, wordlen=10)
        C = align.AlignParams(
            ... # snip
        )
        G = assembly.OverlapBuilder(I, C).build()
        G.save(path)

    All arguments to the constructor become class attributes with the same
    name.

    Attributes:
        index (tuples.Index): A tuples index that responds to
            :func:`align.tuples.Index.seeds`.
        align_params (pw.AlignParams): The alignment parameters for the
            rolling alignment.
        hp_condenser (homopolymeric.HpCondenser): If specified,
            all alignments will be performed in condensed alphabet.
            Consequently, all other arguments are interpretted in the condensed
            alphabet.
        drop_threshold (float): What constitutes a drop in the
            score from one window to the next, default is 0. This means that
            if the overall score does not strictly increase (or the score of
            the new window is not positive) we drop the seed.
        window (int): The size of the rolling window.
        max_succ_drops (int): Maximum number of "drops" until the
            segment is dropped, default is 3.
    """

    def __init__(self, index, align_params, **kwargs):
        self.hp_condenser = kwargs.get('hp_condenser', None)
        if self.hp_condenser:
            assert(isinstance(self.hp_condenser, homopolymeric.HpCondenser))
        self.index, self.align_params = index, align_params
        self.window = kwargs.get('window', 20)
        self.drop_threshold = kwargs.get('drop_threshold', 0)
        self.max_succ_drops = kwargs.get('max_succ_drops', 3)

    def build(self, profile=False):
        """Builds a weighted, directed graph by using tuple methods. The
        process has 2 steps:

        * Find all seeds using :func:`align.tuples.Index.seeds`,
        * Extend all seeds to suffix-prefix segments using :func:`extend`.

        The resulting graph may not necessarily be acyclic. For further
        processing (e.g to find the layout) we need to ensure the overlap
        graph is acyclic. For this, see :func:`OverlapGraph.break_cycles`.

        Keyword Args:
            profile (Optional[bool]): If truthy, instead of reporting
                percentage progress time consumption is reported at *every*
                step (for every pair of sequences). This generates *a lot* of
                output.

        Returns:
            assembly.OverlapGraph: The overlap graph, potentially containing
                cycles.
        """
        vs = set()
        es, ws = [], []
        seqinfo = self.index.tuplesdb.seqinfo()
        seqids = seqinfo.keys()
        num_pairs = (len(seqids) * (len(seqids)-1)) / 2
        indicator = ProgressIndicator('discovering overlaps', num_pairs)
        progress_cnt = 0
        if not profile:
            indicator.start()
        for sid_idx in range(len(seqids)):
            for tid_idx in range(sid_idx + 1, len(seqids)):
                S_id, T_id = seqids[sid_idx], seqids[tid_idx]
                S_info, T_info = seqinfo[S_id], seqinfo[T_id]
                S_min_idx, T_min_idx = S_info['start'], T_info['start']
                S_max_idx = S_info['start'] + S_info['length']
                T_max_idx = T_info['start'] + T_info['length']
                S_name = '%s %d-%d #%d' \
                    % (S_info['name'], S_min_idx, S_max_idx, S_id)
                T_name = '%s %d-%d #%d' \
                    % (T_info['name'], T_min_idx, T_max_idx, T_id)

                if profile:
                    sys.stderr.write('"%s" and "%s": ' % (S_name, T_name))
                else:
                    indicator.progress()
                vs = vs.union([S_name, T_name])

                # do they have any seeds in common?
                _t_seeds = time.time()
                seeds = self.index.seeds(S_id, T_id)
                if profile:
                    _t_seeds = 1000 * (time.time() - _t_seeds)
                    sys.stderr.write(
                        'found %d seeds (%.0f ms)' % (len(seeds), _t_seeds)
                    )
                    if not seeds:
                        sys.stderr.write('.\n')

                if not seeds:
                    continue

                # are the seeds part of an overlap?
                S = self.index.tuplesdb.loadseq(S_id)
                T = self.index.tuplesdb.loadseq(T_id)
                if self.hp_condenser:
                    seeds = [self.hp_condenser.condense_seed(S, T, seed) for seed in seeds]
                    S = self.hp_condenser.condense_sequence(S)
                    T = self.hp_condenser.condense_sequence(T)

                _t_extend = time.time()
                overlap = self.extend(S, T, seeds)
                if profile:
                    _t_extend = 1000 * (time.time() - _t_extend)
                    sys.stderr.write(
                        ' overlaps (%.0f ms): %s\n'
                        % (_t_extend, '+' if overlap else '-')
                    )

                if not overlap:
                    continue

                S_len = self._S_len(overlap.tx.opseq)
                T_len = self._T_len(overlap.tx.opseq)
                if abs(overlap.tx.idx_S - overlap.tx.idx_T) < self.window or \
                    abs(overlap.tx.idx_S + S_len - (overlap.tx.idx_T + T_len)) < self.window:
                    # end points are too close, ignore
                    continue
                if overlap.tx.idx_S == 0 and overlap.tx.idx_T == 0:
                    if S_len < T_len:
                        es += [(S_name, T_name)]
                    elif S_len > T_len:
                        es += [(T_name, S_name)]
                elif overlap.tx.idx_T == 0:
                    es += [(S_name, T_name)]
                elif overlap.tx.idx_S == 0:
                    es += [(T_name, S_name)]
                else:
                    raise RuntimeError("This should not have happened")

                ws += [overlap.tx.score]

        if profile:
            sys.stderr.write('\n')
        else:
            indicator.finish()

        G = OverlapGraph()
        G.iG.add_vertices(list(vs))
        es = [(G.iG.vs.find(name=u), G.iG.vs.find(name=v)) for u, v in es]
        G.iG.add_edges(es)
        G.iG.es['weight'] = ws
        return G

    def _extend_fwd_once(self, S, T, segment, window):
        """Helper method for ``extend1d``."""
        S_len = self._S_len(segment.tx.opseq)
        T_len = self._T_len(segment.tx.opseq)
        align_problem_kw = {
            'S': S, 'T': T, 'params': self.align_params,
            'align_type': pw.START_ANCHORED_OVERLAP,
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
            'align_type': pw.END_ANCHORED_OVERLAP,
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
                S_end = cur_seg.tx.idx_S + self._S_len(cur_seg.tx.opseq)
                T_end = cur_seg.tx.idx_T + self._T_len(cur_seg.tx.opseq)
                S_wiggle = S.length - S_end
                T_wiggle = T.length - T_end
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

            if (len(score_history) == self.max_succ_drops and
                all([x <= self.drop_threshold for x in score_history])):
                break

            cur_seg = seg

        return None

    def extend(self, S, T, segments):
        """Given two sequences and a number of matching segments returns the
        first fully extended segment. A fully extended segment is one that
        corresponds to an suffix-prefix alignment of the sequences.

        Each seed is extended by a rolling frame of size :attr:`window` along
        the two sequences and is dropped once as many as :attr:`max_succ_drops`
        successive bad scores are observed. A bad score (over a single frame)
        is one that has a score of smaller than :attr:`drop_threshold`.

        Args:
            S (seq.Sequence): The "from" sequence.
            T (seq.Sequence): The "to" sequence.
            segments (List[tuples.Segment]): The starting segments. If called
                from :func:`build`, these are seeds but no assumption is made.

        Returns:
            tuples.Segment: A segment corresponding to an overlap alignment.
        """
        assert(S.alphabet.letters == T.alphabet.letters)
        if self.hp_condenser:
            assert S.alphabet.letters == self.hp_condenser.dst_alphabet.letters
        res = []
        for segment in segments:
            fwd = self._extend1d(S, T, segment)
            bwd = self._extend1d(S, T, segment, backwards=True)
            if (fwd and bwd and
                min(fwd.tx.score, bwd.tx.score) > self.drop_threshold):
                assert(bwd.tx.idx_S == 0 or bwd.tx.idx_T == 0)
                opseq = bwd.tx.opseq[:-len(segment.tx.opseq)] + fwd.tx.opseq
                score = self.align_params.score(
                    S, T, opseq,
                    S_min_idx=bwd.tx.idx_S, T_min_idx=bwd.tx.idx_T
                )
                tx = pw.Transcript(
                    idx_S=bwd.tx.idx_S, idx_T=bwd.tx.idx_T,
                    score=score, opseq=opseq
                )
                return tuples.Segment(
                    S_id=segment.S_id, T_id=segment.T_id, tx=tx
                )
