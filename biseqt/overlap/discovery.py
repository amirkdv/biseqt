from collections import namedtuple
from math import log
import os.path
import sys
from .. import pw, words, lib, ffi, ProgressIndicator, CffiObject

class SeedExtensionParams(CffiObject):
    """Wraps the C struct ``seedext_params``, see ``pwlib.h``.

    Attributes:
        scores (pw.AlignScores): Substitution and gap scores.
        window (int): The width of the rolling window of alignment.
        min_score (float): The minimum required score for an overlap to be
            reported.
        max_new_mins (int): Maximum number of times a new minimum score is
            tolerated before alignment is aborted.
    """
    def __init__(self, **kw):
        if 'c_obj' in kw:
            self.c_obj = kw['c_obj']
        else:
            self.c_obj = ffi.new('seedext_params*', {
                'window': kw['window'],
                'min_score': kw['min_score'],
                'max_new_mins': kw['max_new_mins'],
                'scores': kw['scores'].c_obj
            })

def most_significant_shift(S_id, T_id, index, min_overlap=-1):
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
    seeds = index.seeds(S_id, T_id)
    if not seeds:
        return None, None
    seqinfo = index.seqdb.seqinfo()
    S_len, T_len = seqinfo[S_id]['length'], seqinfo[T_id]['length']
    scores = {}
    for seed in seeds:
        shift = seed.tx.S_idx - seed.tx.T_idx
        r, s, s0 = index.seed_pvalue_contribution((S_len, T_len), shift)
        for d in range(shift - r, shift + r + 1):
            scores[d] = s0 if d not in scores else scores[d]
            scores[d] += s

    if min_overlap > 0:
        for d in range(- T_len, - T_len + min_overlap):
            scores.pop(d, 0)
        for d in range(S_len - min_overlap, S_len):
            scores.pop(d, 0)

    if not scores:
        return (None, None)

    return max(scores.items(), key=lambda x: scores[x[0]])

# FIXME docs
def discover_overlap(S_id, T_id, index, mode='banded alignment', **kwargs):
    """
    Args:
        S_id (int):
        T_id (int):

    Keyword Args:
        rw_collect (Optional[bool]):
        index (words.Index):
        params (OverlapDiscoveryParams):
    """
    assert(mode in ['seed extension', 'banded alignment'])
    min_sig = kwargs['min_shift_significance']
    seqinfo = index.seqdb.seqinfo()
    S_len, T_len = seqinfo[S_id]['length'], seqinfo[T_id]['length']
    seeds = index.seeds(S_id, T_id)

    if not seeds:
        return None

    shift, sig = most_signifcant_shift(S_len, T_len, index)
    if sig is None or sig < min_sig:
        return None

    S, T = index.seqdb.loadseq(S_id), index.seqdb.loadseq(T_id)
    tx = None
    if mode == 'banded_alignment':
        radius = index.band_radius((S_len, T_len), shift, kwargs['gap_prob'])
        F = pw.AlignFrame(S, T)
        assert(isinstance(kwargs['aln_scores'], AlignScores))
        with pw.AlignTable(F, kwargs['aln_scores'], alnmode=BANDED_MODE, alntype=B_OVERLAP) as T:
            score = T.solve()
            if score is not None and score > kwargs['min_alignment_score']:
                tx = T.traceback()
    elif mode == 'seed extension':
        seeds = words.maximal_seeds(seeds, S_id, T_id)
        seeds = sorted(seeds, key=lambda x: abs(x.tx.S_idx - x.tx.T_idx - best_shift))
        txs = ffi.new('transcript* []', [seed.tx.c_obj for seed in seeds])
        F = pw.AlignFrame(S, T)
        assert(isinstance(kwargs['seed_ext_params'], SeedExtensionParams))
        tx = lib.extend(txs, len(txs), F.c_obj, kwargs['seed_ext_params'].c_obj)
        tx = None if tx == ffi.NULL else tx

    # starting diagonal of the alignment must be far enough from the 0 diagonal:
    if tx and abs(tx.S_idx - tx.T_idx) < kwargs.get('min_margin', 1000):
        return None

    return tx

