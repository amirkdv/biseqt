from collections import namedtuple
from math import log
from .. import pw, lib, ffi

# TODO make this a C struct so the C code cleans up.
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
    pass

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

def analyze_shifts(seeds, S_len, T_len, rolling_sum_width):
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
    def window_log_pvalue(shift, num):
        log_pvalue = -log(S_len) - log(T_len) + log(rolling_sum_width) \
            + 0.5 * log( (S_len - abs(shift))**2 + (T_len - abs(shift))**2 )
        # 1- we have num observations (each a seed) with the same probability
        #    of being matched by the null hypothesis
        # 2- we are testing S_len+T_len simultaneous hypotheses;
        #    apply a Bonferroni correction:
        return log(S_len + T_len) + num * log_pvalue

    def rolling_sum(data):
        if len(data) < rolling_sum_width:
            return
        cur = 0
        for idx in range(0, len(data)):
            if idx >= rolling_sum_width:
                cur -= data[idx - rolling_sum_width]
            if idx < len(data):
                cur += data[idx]
            yield cur

    shift_range = range(-T_len, S_len)
    shift_coverage = {shift:0 for shift in shift_range}
    for seed in seeds:
        # all seeds are the same length at this stage (and they
        # are potentially overlapping):
        shift_coverage[seed.tx.S_idx - seed.tx.T_idx] += 1

    shift_distrib = [x for x in rolling_sum([x[1] for x in sorted(shift_coverage.items())])]
    mode_idx, mode = max(enumerate(shift_distrib), key=lambda x: x[1])
    mode_shift = shift_range[min(mode_idx, len(shift_range)-1)]
    log_pvalue = window_log_pvalue(mode_shift, mode)

    return (mode_shift, log_pvalue)
