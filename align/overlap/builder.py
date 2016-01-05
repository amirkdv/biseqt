import re
import sys
import time
from .. import tuples, pw, homopolymeric, ProgressIndicator, lib, ffi
from . import OverlapGraph, SeedExtensionParams, extend_segments, most_signitifcant_shift

class OverlapBuilder(object):
    """Provided a :class:`align.tuples.Index` builds an overlap graph of all
    the sequences. All sequences must have already been indexed in the
    tuples database. For example::

        B = tuples.TuplesDB('path/to/file.db', alphabet=seq.Alphabet("ACGT"))
        I = tuples.Index(B, wordlen=10)
        C = align.AlignParams(
            ... # snip
        )
        G = overlap.OverlapBuilder(I, align_params=C, window=50,
            max_new_mins=10, min_overlap_score=400)
        G.build().save(path)

    All arguments to the constructor become class attributes with the same
    name.

    Overlaps are found as follows for any pair of potentially overlapping reads:

    * Try to decide if the reads overlap based on their shift distribution,
      see :func:`overlap_by_seed_shift_distribution`. The *shift* of a
      seed with coordinates :math:`(i_S,i_T)` is the integer :math:`i_S-i_T`.
    * If previous step was nonconclusive, try to extend seeds with shifts close
      to the shift mode, see :func:`overlap_by_seed_extension`.

    Attributes:
        index (tuples.Index): A tuples index that responds to
            :func:`seeds() <align.tuples.Index.seeds>`.
        align_params (pw.AlignParams): The alignment parameters for the
            rolling alignment.
        hp_condenser (homopolymeric.HpCondenser): If specified,
            all alignments will be performed in condensed alphabet.
            Consequently, all other arguments are interpretted in the condensed
            alphabet.
        window (int): The size of the rolling window.
        min_overlap_score (float): The minimum required score for an alignment
            to be reported.
        min_margin (int): The minimum margin required for the direction of an
            overlap to be reliable; default is :attr:`window`.
        shift_rolling_sum_width (int): The width of the rolling sum used to
            find the mode of the shifts distribution.
        rw_collect (Optional[bool]): Whether to ask libalign to dump score
            random walk data for all tried extensions to 'scores.txt';
            default is False.
    """
    def __init__(self, index, **kwargs):
        self.hp_condenser = kwargs.get('hp_condenser', None)
        self.index = index

        self.seed_ext_params = {k:kwargs[k] for k in ['window', 'min_overlap_score', 'max_new_mins', 'align_params']}
        self.seed_ext_params = SeedExtensionParams(**self.seed_ext_params)

        self.shift_rolling_sum_width = kwargs['shift_rolling_sum_width']
        self.lower_log_pvalue_cutoff = kwargs['lower_log_pvalue_cutoff']
        self.upper_log_pvalue_cutoff = kwargs['upper_log_pvalue_cutoff']
        self.min_margin = kwargs['min_margin']
        self.rw_collect = bool(kwargs.get('rw_collect', False))
        self.seqinfo = self.index.tuplesdb.seqinfo()

    def build(self):
        """Builds a weighted, directed graph by using tuple methods. The
        process has 3 steps:

        * Find all seeds using :func:`align.tuples.Index.seeds`,
        * Try to judge whether they are overlapping by judging their shift
          distribution; see :func:`overlap_by_seed_shift_distribution`.
        * If previous step was nonconclusive, extend select seeds (those with
          shift close to the shift mode found above) to suffix-prefix segments;
          see :func:`overlap_by_seed_extension`.

        The resulting graph may not necessarily be acyclic. For further
        processing (e.g to find the layout) we need to ensure the overlap
        graph is acyclic. For this, see :func:`align.overlap.OverlapGraph.break_cycles`.

        Returns:
            align.overlap.OverlapGraph: The overlap graph, potentially containing
                cycles.
        """
        process_time_spent = {
            'extension': {'t.p.': 0, 'f.p.': 0, 't.n.': 0, 'f.n.': 0},
            'p-values': 0,
            'seeds': 0
        }
        vs = set()
        es, ws = [], []
        seqids = self.seqinfo.keys()
        msg = 'Extending seeds on potentially homologous sequences'
        indicator = ProgressIndicator(msg,
            self.index.num_potential_homolog_pairs(), percentage=False)
        num_shift_decided = 0
        indicator.start()
        for S_id in seqids:
            for T_id in self.index.potential_homologs(S_id):
                S_info, T_info = self.seqinfo[S_id], self.seqinfo[T_id]
                S_min_idx, T_min_idx = S_info['start'], T_info['start']
                S_max_idx = S_info['start'] + S_info['length']
                T_max_idx = T_info['start'] + T_info['length']
                S_name = '%s %d-%d #%d' \
                    % (S_info['name'], S_min_idx, S_max_idx, S_id)
                T_name = '%s %d-%d #%d' \
                    % (T_info['name'], T_min_idx, T_max_idx, T_id)

                indicator.progress()
                vs = vs.union([S_name, T_name])

                # used to profile false/true positive/negatives.
                cheat_overlaps = not (S_max_idx < T_min_idx or T_max_idx < S_min_idx)
                # if not cheat_overlaps:
                #     continue

                # do they have any seeds in common?
                _t = time.clock()
                seeds = self.index.seeds(S_id, T_id)
                true_shift = T_min_idx - S_min_idx
                process_time_spent['seeds'] += 1000 * (time.clock() - _t)

                if not seeds:
                    continue

                # do the seeds obviously (statistically) indicate an overlap?
                _t = time.clock()
                overlap = self.overlap_by_seed_shift_distribution(seeds, S_id, T_id)
                process_time_spent['p-values'] += 1000 * (time.clock() - _t)
                if isinstance(overlap, pw.Segment):
                    # FIXME some of these are later discarded because of margins,
                    # we then get things like "10 out of 9 edges were decided by shift distribution"
                    num_shift_decided += 1
                if isinstance(overlap, list):
                    seeds = overlap
                    _t = time.clock()
                    seeds = tuples.Index.maximal_seeds(seeds, S_id, T_id)

                    process_time_spent['seeds'] += 1000 * (time.clock() - _t)
                    _t = time.clock()
                    overlap = self.overlap_by_seed_extension(seeds, S_id, T_id)
                    t_extension = 1000 * (time.clock() - _t)
                    if overlap:
                        if cheat_overlaps:
                            process_time_spent['extension']['t.p.'] += t_extension
                        else:
                            process_time_spent['extension']['f.p.'] += t_extension
                    else:
                        if cheat_overlaps:
                            process_time_spent['extension']['f.n.'] += t_extension
                        else:
                            process_time_spent['extension']['t.n.'] += t_extension


                if not overlap:
                    continue

                S_len = lib.tx_seq_len(overlap.tx.c_obj, 'S')
                T_len = lib.tx_seq_len(overlap.tx.c_obj, 'T')
                assert(overlap.tx.T_idx * overlap.tx.S_idx == 0)
                lmargin = abs(overlap.tx.S_idx - overlap.tx.T_idx)
                rmargin = abs(overlap.tx.S_idx + S_len - (overlap.tx.T_idx + T_len))
                if lmargin < self.min_margin or \
                    (lmargin == 0 and rmargin < self.min_margin):
                    # end points are too close, ignore
                    continue
                if overlap.tx.S_idx == 0 and overlap.tx.T_idx == 0:
                    if S_len < T_len:
                        es += [(S_name, T_name)]
                    elif S_len > T_len:
                        es += [(T_name, S_name)]
                elif overlap.tx.T_idx == 0:
                    es += [(S_name, T_name)]
                elif overlap.tx.S_idx == 0:
                    es += [(T_name, S_name)]
                else:
                    raise RuntimeError("This should not have happened")

                ws += [overlap.tx.score]

        indicator.finish()
        # Make the process_time_spent dict look like yaml output
        yamlify_depth1 = lambda v: '\n' + '\n'.join( '  %s: %.2f' % (p,q) for p,q in v.items())
        report = '\n'.join(
            '%s: %s' % (k, yamlify_depth1(v) if isinstance(v, dict) else '%.2f' % v) \
            for k,v in process_time_spent.items() \
        )
        sys.stderr.write('* decomposition of time spent (in milliseconds):\n')
        sys.stderr.write('  %s\n' % report.replace('\n', '\n  '))

        G = OverlapGraph()
        G.iG.add_vertices(list(vs))
        es = [(G.iG.vs.find(name=u), G.iG.vs.find(name=v)) for u, v in es]
        G.iG.add_edges(es)
        sys.stderr.write('* %d out of %d edges where chosen by shift distribution\n' % (num_shift_decided, G.iG.ecount()))
        G.iG.es['weight'] = ws
        return G

    # FIXME return a best shift so we can sort them in order (the order dies
    # when we do maximal seeds)
    # FIXME docs are out of date
    def overlap_by_seed_shift_distribution(self, seeds, S_id, T_id):
        """Decides whether the shift distribution of seeds for a given sequence
        pair is "indicative" enough of an overlap or lack thereof.

        Args:
            seeds (list[Segment]): Exactly matching seeds as returned by
                :attr:`index`.
            S_id (int): The database ID of the "from" sequence.
            T_id (int): The database ID of the "to" sequence.

        Returns:
            None|Segment|list[Segment]: Corresponding to scenarios described above.

        """
        S_len, T_len = self.seqinfo[S_id]['length'], self.seqinfo[T_id]['length']
        mode_shift, log_pvalue = most_signitifcant_shift(S_len, T_len,
            seeds, self.shift_rolling_sum_width)

        if log_pvalue > self.upper_log_pvalue_cutoff:
            # definitely not overlapping:
            return None
        overlap_length = min(S_len, T_len + mode_shift) - max(0, mode_shift) if mode_shift < min(S_len,T_len) else 0
        if log_pvalue < self.lower_log_pvalue_cutoff:
            # definitely overlapping; fake the transcripts:
            if mode_shift >= 0:
                tx = pw.Transcript(S_idx=mode_shift, T_idx=0, score=overlap_length, opseq='M'*overlap_length)
                return pw.Segment(S_id=S_id, T_id=T_id, tx=tx)
            else:
                tx = pw.Transcript(S_idx=0, T_idx=-mode_shift, score=overlap_length, opseq='M'*overlap_length)
                return pw.Segment(S_id=S_id, T_id=T_id, tx=tx)

        # Only return those seeds that have a shift close to the mode:
        # FIXME should we return only some of the seeds?
        # seeds = [seed for seed in seeds if abs(seed.tx.S_idx - seed.tx.T_idx - mode_shift) < 10 * self.shift_rolling_sum_width]
        return sorted(seeds, key=lambda x: abs(x.tx.S_idx - x.tx.T_idx - mode_shift))

    def overlap_by_seed_extension(self, seeds, S_id, T_id):
        """Tries to find a seed among given seeds that extends to a full overlap alignment by consecutive
        start/end-anchored overlap alignments in a moving window along the two
        reads:

          * If any such seed is found, the fully extended segment is returned.
          * If no such seeds are found, ``None`` is returned.

        Args:
            seeds (list[Segment]): Exactly matching seeds as returned by
                :attr:`index`.
            S_id (int): The database ID of the "from" sequence.
            T_id (int): The database ID of the "to" sequence.

        Returns:
            Segment|None: Depending on whether any seed successfully extends to boundaries.
        """
        S = self.index.tuplesdb.loadseq(S_id)
        T = self.index.tuplesdb.loadseq(T_id)
        # Calculate the score of each seed; FIXME do we need this?
        for seed in seeds:
            seed.tx.score = self.seed_ext_params.align_params.score(
                S, T, seed.tx.opseq,
                S_min_idx=seed.tx.S_idx, T_min_idx=seed.tx.T_idx
            )
        if self.hp_condenser:
            # condense the sequences and their seeds:
            S = self.hp_condenser.condense_sequence(S)
            T = self.hp_condenser.condense_sequence(T)
            seeds = (self.hp_condenser.condense_seed(S, T, s) for s in seeds)
            seeds = filter(lambda x: x, seeds)

        return extend_segments(S, T, seeds, self.seed_ext_params, rw_collect=self.rw_collect)
