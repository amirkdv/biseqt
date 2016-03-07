from collections import namedtuple
from math import log
import os.path
import sys
from .. import pw, words, lib, ffi, ProgressIndicator
from . import OverlapGraph

# TODO make this a C struct so the C code cleans up.
# FIXME the documentation formatting is weird
class SeedExtensionParams(namedtuple('SeedExtensionParams',
    ['window', 'min_overlap_score', 'max_new_mins', 'align_params'])):
    """Represents the set of tuning parameters for seed extension.

    Attributes:
        window (int): The size of the rolling window.
        min_overlap_score (float): The minimum required score for an alignment
            to be reported
        max_new_mins (int): Maximum number of new minima observed in the score
            random walk until a segement is dropped.
        align_params (pw.AlignParams): The alignment parameters for the
            rolling alignment.
    """

# FIXME docs
OverlapDiscoveryParams = namedtuple('OverlapDiscoveryParams', [
    'seed_ext_params', 'shift_rolling_sum_width'
])

def extend_segments(S, T, segments, params, rw_collect=False):
    """Wraps :c:func:`extend()`: given two sequences and a number of
    matching segments returns the first fully extended segment.

    Args:
        S (seq.Sequence): The "from" sequence.
        T (seq.Sequence): The "to" sequence.
        segments (list[pw.Segment]): The starting segments. If called
            from :func:`build`, these are seeds but no assumption is made.
        params (SeedExtensionParams): The parameters for seed extension.

    Returns:
        pw.Segment|None: A segment corresponding to an overlap alignment or
            None if none found.
    """
    if not segments:
        return None

    segs = ffi.new('segment* []', [seg.c_obj for seg in segments])
    res = lib.extend(
        segs, len(segs),
        S.c_idxseq, T.c_idxseq, len(S), len(T), params.align_params.c_obj,
        params.window, params.max_new_mins, params.min_overlap_score, int(bool(rw_collect))
    )
    return pw.Segment(c_obj=res) if res != ffi.NULL else None

# FIXME docs
def rolling_sum(data, width):
    if len(data) < width:
        return
    cur = 0
    for idx in range(len(data)):
        if idx >= width:
            cur -= data[idx - width]
        if idx < len(data):
            cur += data[idx]
        # FIXME yield the corresponding indices so we can double check
        yield cur

# FIXME docs
class ShiftWindow(namedtuple('ShiftWindow', ['S_len', 'T_len', 'width', 'shift'])):
    # FIXME docs
    def significance(self, num):
        diag_len_squared = (self.S_len - abs(self.shift))**2 + \
            (self.T_len - abs(self.shift))**2
        log_pvalue = - log(self.S_len) - log(self.T_len) \
            + log(self.width) \
            + 0.5 * log(diag_len_squared)
        # we have num observations (each a seed) with the same probability
        # of being matched by the null hypothesis:
        log_pvalue *= num
        # we are testing S_len+T_len simultaneous hypotheses; apply a Bonferroni
        # correction:
        log_pvalue += log(self.S_len + self.T_len)
        return log_pvalue

# FIXME find the strip with max number of shifts *per unit area*.
# FIXME make some of these keyword arguments
def most_signifcant_shift(S_len, T_len, seeds, rolling_sum_width):
    """Builds a smoothed distribution (via a rolling sum of known width) of
    shifts for the provided seeds of a pair of sequences of known lengths and
    returns the shift with most number of seeds and its p-value.

    The shift window corresponding to shift :math:`s` is the interval
    :math:`[s-w, s]` of shift values where :math:`w` is the rolling sum width.

    The p-value for a given shift is calculated as follows: consider the
    distribution of seeds in the :math:`(i_S,i_T)` plane where :math:`i_S,i_T`
    are the starting coordinates of each seed. In this plane, all seeds reside
    within the rectangle :math:`M = [0, |S|] \\times [0, |T|]`. We assume that
    under the null hypothesis (the two sequences are not overlapping) seeds
    are uniformly distributed over :math:`M`. Further, in this plane a shift
    window :math:`[s-w, s]` corresponds to a diagonal strip at horizontal
    and vertical distance :math:`s` of the main quadrant diagonal. The p-value
    is then calculated as the probability that as many seeds would be observed
    in the corresponding strip:

        :math:`\\Pr(X_s\\ge n) \\simeq \\left(\\frac{A_s}{|S|\cdot|T|}\\right)^n`

    where :math:`A_s` is the area of the diagonal :math:`w`-strip corresponding
    to shift :math:`s`:

        :math:`A_s = w\\sqrt{(|S|-|s|)^2 + (|T|-|s|)^2}`

    The final formula for the returned pvalue logarithm is:

        :math:`\\log(|S|+|T|) + n \\left[\\log(A_s) - \\log(|S|) - \\log(|T|)\\right]`

    where the first term is a Bonferroni correction since the process involves
    testing the same hypothesis for :math:`|S|+|T|` different strips of shift
    windows.

    Args:
        seeds (list[Segment]): Exactly matching seeds.
        S_len (int): The length of the "from" sequence.
        T_len (int): The length of the "to" sequence.
        rolling_sum_width(int): The parameter :math:`w` in the above.

    Returns:
        (mode_shift, log_pvalue): A tuple of the form ``(int, float)``.

    """
    shift_range = range(-T_len, S_len)
    shift_coverage = {shift:0 for shift in shift_range}
    for seed in seeds:
        # all seeds are the same length at this stage (and they
        # are potentially overlapping):
        shift_coverage[seed.tx.S_idx - seed.tx.T_idx] += 1

    shift_coverage = [x[1] for x in sorted(shift_coverage.items())]
    shift_distrib = rolling_sum(shift_coverage, rolling_sum_width)
    mode_idx, mode = max(enumerate(shift_distrib), key=lambda x: x[1])
    most_dense_window = ShiftWindow(S_len=S_len, T_len=T_len,
        width=rolling_sum_width,
        shift=shift_range[min(mode_idx, len(shift_range)-1)]
    )

    return most_dense_window.shift, most_dense_window.significance(mode)

# FIXME docs
def discover_overlap(S_id, T_id, rw_collect=False, **kwargs):
    """
    Args:
        S_id (int):
        T_id (int):

    Keyword Args:
        rw_collect (Optional[bool]):
        index (words.Index):
        params (OverlapDiscoveryParams):
    """
    index, params = kwargs['index'], kwargs['params']
    seqinfo = index.seqdb.seqinfo()
    S_len, T_len = seqinfo[S_id]['length'], seqinfo[T_id]['length']
    seeds = index.seeds(S_id, T_id)

    if not seeds:
        return None

    # TODO is there any use to the significance value itself?
    best_shift, _ = most_signifcant_shift(S_len, T_len, seeds,
        params.shift_rolling_sum_width)

    seeds = words.maximal_seeds(seeds, S_id, T_id)
    seeds = sorted(seeds, key=lambda x: abs(x.tx.S_idx - x.tx.T_idx - best_shift))

    S = index.seqdb.loadseq(S_id)
    T = index.seqdb.loadseq(T_id)

    seg = extend_segments(S, T, seeds, params.seed_ext_params, rw_collect=rw_collect)
    if not seg:
        return None

    assert(seg.tx.T_idx * seg.tx.S_idx == 0)
    return seg

# FIXME docs, simplify?
def _satisfies_min_margin(seg, S_tx_len, T_tx_len, min_margin):
    assert(seg.tx.S_idx * seg.tx.T_idx == 0)
    lmargin = abs(seg.tx.S_idx - seg.tx.T_idx)
    rmargin = abs(seg.tx.S_idx + S_tx_len - (seg.tx.T_idx + T_tx_len))
    if lmargin < min_margin:
        return False
    if lmargin == 0 and rmargin < min_margin:
        return False
    return True

# FIXME docs
def _overlap_direction(seg, S_tx_len, T_tx_len):
    assert(seg.tx.S_idx * seg.tx.T_idx == 0)

    if seg.tx.S_idx == 0 and seg.tx.T_idx == 0:
        # Edge case: the two sequences align with no gap at (0,0):
        return '+' if S_tx_len < T_tx_len else '-'
        # TODO what to do with the case where S_tx_len = T_tx_len?
    else:
        # We know exactly one of `seg.tx.{S_idx,T_idx}` is zero:
        return '+' if seg.tx.T_idx == 0 else '-'

def overlap_graph(index, od_params, min_margin=10, rw_collect=False):
    vs = set()
    es, ws = [], []
    seqinfo = index.seqdb.seqinfo()
    seqids = seqinfo.keys()
    seqids = [1,2,3] # FIXME
    msg = 'Extending seeds on potentially homologous sequences'
    indicator = ProgressIndicator(msg,
        index.num_potential_homolog_pairs(), percentage=False)
    num_shift_decided = 0
    indicator.start()
    for S_id in seqids:
        for T_id in index.potential_homologs(S_id):
            if T_id not in [1,2,3]:
                continue
            indicator.progress()
            S_name = '%s #%d' % (seqinfo[S_id]['name'], S_id)
            T_name = '%s #%d' % (seqinfo[T_id]['name'], T_id)
            vs = vs.union([S_name, T_name])

            seg = discover_overlap(S_id, T_id, index=index,
                params=od_params, rw_collect=rw_collect)

            if not seg:
                continue

            S_tx_len = lib.tx_seq_len(seg.tx.c_obj, 'S')
            T_tx_len = lib.tx_seq_len(seg.tx.c_obj, 'T')

            if not _satisfies_min_margin(seg, S_tx_len, T_tx_len, min_margin):
                continue

            d = _overlap_direction(seg, S_tx_len, T_tx_len)
            es += [(S_name, T_name)] if d == '+' else [(T_name, S_name)]
            ws += [seg.tx.score]

    indicator.finish()

    G = OverlapGraph()
    G.iG.add_vertices(list(vs))
    es = [(G.iG.vs.find(name=u), G.iG.vs.find(name=v)) for u, v in es]
    G.iG.add_edges(es)
    G.iG.es['weight'] = ws
    return G