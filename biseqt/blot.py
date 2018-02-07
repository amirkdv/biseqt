# -*- coding: utf-8 -*-
# TODO
"""
Diagonals are numbered as follows in the dynamic programming table::

    0 -1 -2 -3 ... -len1
    1
    2
    3
    .
    .
    .
    +len0

"""
import numpy as np
from scipy.special import erfcinv


# A band (i-r, i+r) is considered of interest if xs[i] >threshold. This
# function returns a maximal (disjoint) set of bands of interest (i.e if two
# bands overlap they are reported as one bigger band).
def find_peaks(xs, rs, threshold):
    """Finds maximal (disjoint) peak regions in a sequence of real numbers.
    Each value that is at least as large as the threshold constitutes the
    center of a peak with radius according to its position. In the output all
    overlapping peaks are merged into maximal peaks.

    Args:
        xs: the 1D data sequence of interest
        rs: the radii for peaks defined at every point exceeding threshold,
            could be a 1D sequence or a fixed number,
        threshold: cutoff value to compare with ``xs`` values.

    Returns:
        list (tuple): A least of "peaks", each a tuple of ``(left, right)``
        coordinates in ``xs``. Returned peaks are guaranteed to be disjoint.
    """
    peaks = []
    cur_peak = None
    for idx, x in enumerate(xs):
        radius = rs[idx] if isinstance(rs, list) else rs
        assert isinstance(radius, int)
        if x < threshold:
            continue
        peak_l = max(0, idx - radius)
        peak_r = min(len(xs) - 1, idx + radius)
        if cur_peak is None:
            cur_peak = (peak_l, peak_r)
            continue
        if peak_l < cur_peak[1]:  # overlaps with cur_peak
            assert peak_r >= cur_peak[1]
            cur_peak = (cur_peak[0], peak_r)
        else:
            peaks.append(cur_peak)
            cur_peak = (peak_l, peak_r)
    peaks.append(cur_peak)
    return [(int(l), int(r)) for (l, r) in peaks]


def expected_overlap_len(len0, len1, diag, gap_prob):
    """Calculates the expected length of an overlap alignment given its starting
    coordinates:

    .. math::
        K = \\left(\\frac{2}{2 - g}\\right) L

    where :math:`L` is the maximum possible length of the alignment:

    .. math::
        L = \\min(l_0 - d, l_1) + \\min(d, 0)

    with :math:`l_0,l_1` being the length of the sequences (i.e ``len0`` and
    ``len1`` arguments) and :math:`d` the starting diagonal (i.e ``diag``
    argument).

    Args:
        len0 (int): Length of the 1st sequence.
        len1 (int): Length of the 2nd sequence.
        diag (int): Starting diagonal of alignments to consider.
        gap_prob (float): Probability of indels occuring at any position.
    Returns:
        int: Expected length of an overlap alignment.
    """
    max_len = min(len0 - diag, len1) + min(diag, 0)
    expected_len = (2. / (2 - gap_prob)) * max_len
    assert expected_len >= 0
    return int(np.ceil(expected_len))


# band radius for edit path of length K
def band_radius(expected_len, gap_prob, sensitivity):
    """Calculates the smallest band radius in the dynamic programming table
    such that an alignment of given expected length, with the given gap
    probability, stays entirely within the diagonal band with probability given
    as sensitivity. This is given by:

    .. math::
        r = 2\\sqrt{g(1-g)K}
            \\mathrm{erf}^{-1}\\left(1-\\epsilon\\right)

    where :math:`g` is the gap probability, :math:`1-\\epsilon` is the desired
    sensitivity, and :math:`K` is the given expected length.

    Args:
        expected_len (int):
            minimum expected length of similar region.
        gap_prob (float): Probability of indels occuring at any position.
        sensitivity (float): The probability that an alignment with given gap
            probability remains entirely within the band.
    Returns:
        int: The smallest band radius guaranteeing the required sensitivity.
    """
    assert 0 < gap_prob < 1 and 0 < sensitivity < 1
    epsilon = 1. - sensitivity
    C = 2 * erfcinv(epsilon) * np.sqrt(gap_prob * (1 - gap_prob))
    radius = C * np.sqrt(expected_len)
    return max(1, int(np.ceil(radius)))


def band_radii(expected_lens, gap_prob, sensitivity):
    assert 0 < gap_prob < 1 and 0 < sensitivity < 1
    epsilon = 1. - sensitivity
    C = 2 * erfcinv(epsilon) * np.sqrt(gap_prob * (1 - gap_prob))
    return np.array([max(1, int(np.ceil(C * np.sqrt(K))))
                     for K in expected_lens])


