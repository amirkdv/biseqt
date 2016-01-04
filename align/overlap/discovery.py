from .. import pw, lib, ffi
from collections import namedtuple

class SeedExtensionParams(namedtuple('SeedExtensionParams',
    ['window', 'min_overlap_score', 'max_new_mins', 'align_params'])):
    """Represents the set of tuning parameters for seed extension.

    Attrs:
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
