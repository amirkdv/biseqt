import re
import sys
import time
from collections import namedtuple
from .. import tuples, pw, homopolymeric, ProgressIndicator, lib, ffi
from . import OverlapGraph, SeedExtensionParams, OverlapDiscoveryParams, extend_segments, most_signifcant_shift, discover_overlap

# FIXME this should be merged into graph.py as a stand alone function
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
    def __init__(self, **kwargs):
        self.index, self.min_margin = kwargs['index'], kwargs['min_margin']
        names = ['window', 'min_overlap_score', 'max_new_mins', 'align_params']
        seed_ext_params = SeedExtensionParams(**{k:kwargs[k] for k in names})
        hp_condenser = kwargs.get('hp_condenser', None)
        overlap_discover_args = {

        }
        self.params = OverlapDiscoveryParams(
            hp_condenser=hp_condenser,
            seed_ext_params=seed_ext_params,
            shift_rolling_sum_width=kwargs['shift_rolling_sum_width']
        )
        self.rw_collect = bool(kwargs.get('rw_collect', False))

    def build(self, true_overlaps=[]):
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
        vs = set()
        es, ws = [], []
        seqinfo = self.index.tuplesdb.seqinfo()
        seqids = seqinfo.keys()
        msg = 'Extending seeds on potentially homologous sequences'
        indicator = ProgressIndicator(msg,
            self.index.num_potential_homolog_pairs(), percentage=False)
        num_shift_decided = 0
        indicator.start()
        for S_id in seqids:
            for T_id in self.index.potential_homologs(S_id):
                indicator.progress()
                S_name = '%s %d' % (seqinfo[S_id]['name'], S_id)
                T_name = '%s %d' % (seqinfo[T_id]['name'], T_id)
                vs = vs.union([S_name, T_name])

                # FIXME OverlapDiscoveryParams does not exist here yet
                overlap = discover_overlap(S_id, T_id, index=self.index,
                    params=self.params, rw_collect=self.rw_collect)
                if not overlap:
                    continue

                S_tx_len = lib.tx_seq_len(overlap.tx.c_obj, 'S')
                T_tx_len = lib.tx_seq_len(overlap.tx.c_obj, 'T')

                lmargin = abs(overlap.tx.S_idx - overlap.tx.T_idx)
                rmargin = abs(overlap.tx.S_idx + S_len - (overlap.tx.T_idx + T_len))
                if lmargin < self.min_margin or \
                    (lmargin == 0 and rmargin < self.min_margin):
                    # end points are too close, ignore
                    continue
                if overlap.tx.S_idx == 0 and overlap.tx.T_idx == 0:
                    if S_tx_len < T_tx_len:
                        es += [(S_name, T_name)]
                    elif S_tx_len > T_tx_len:
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
