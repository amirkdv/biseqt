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

from .seeds import SeedIndex


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
        peak_l, peak_r = max(0, idx - 1), min(len(xs) - 1, idx + 1)
        if cur_peak is None:
            cur_peak = (peak_l, peak_r)
            continue
        if peak_l < cur_peak[1] + radius:  # overlaps with cur_peak
            assert peak_r >= cur_peak[1]
            cur_peak = (cur_peak[0], peak_r)
        else:
            peaks.append(cur_peak)
            cur_peak = (peak_l, peak_r)
    if cur_peak is not None:
        peaks.append(cur_peak)
    return [(int(l), int(r)) for (l, r) in peaks]


def wall_to_wall_distance(len0, len1, diag):
    """Wall to wall distance for a diagonal position

    .. math::
        L = \\min(l_0 - d, l_1) + \\min(d, 0)

    with :math:`l_0,l_1` being the length of the sequences (i.e ``len0`` and
    ``len1`` arguments) and :math:`d` the starting diagonal (i.e ``diag``
    argument).
    """
    return min(len0 - diag, len1) + min(diag, 0)


def expected_overlap_len(len0, len1, diag, gap_prob):
    """Calculates the expected length of an overlap alignment given its starting
    coordinates:

    .. math::
        K = \\left(\\frac{2}{2 - g}\\right) L

    where :math:`L` is the :func:`wall_to_wall_distance`.

    Args:
        len0 (int): Length of the 1st sequence.
        len1 (int): Length of the 2nd sequence.
        diag (int): Starting diagonal of alignments to consider.
        gap_prob (float): Probability of indels occuring at any position.
    Returns:
        int: Expected length of an overlap alignment.
    """
    L = wall_to_wall_distance(len0, len1, diag)
    expected_len = (2. / (2 - gap_prob)) * L
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
    """Same as :func:`band_radius` but for bulk calculations (e.g. finding
    overlaps).

    Args:
        expected_lens (list):
            List of integers for which to calculate band radii.
        gap_prob (float):
            As in :func:`band_radius`.
        sensitivity (float):
            As in :func:`band_radius`.
    """
    assert 0 < gap_prob < 1 and 0 < sensitivity < 1
    epsilon = 1. - sensitivity
    C = 2 * erfcinv(epsilon) * np.sqrt(gap_prob * (1 - gap_prob))
    return np.array([max(1, int(np.ceil(C * np.sqrt(K))))
                     for K in expected_lens])


class HomologyFinder(SeedIndex):
    """A homology finder based on m-dependent CLT statistics.

    Attributes:
        gap_prob (float):
            Probability of indels at a given position in mutation model.
        subst_prob (float):
            Substitition (mismatch) probabilities in mutation model.
        sensitivity (float):
            Desired sensitivity of bands.
    """
    def __init__(self, S, T, gap_prob=None, subst_prob=None, sensitivity=None,
                 **kw):
        self.gap_prob = gap_prob
        self.subst_prob = subst_prob
        self.sensitivity = sensitivity
        super(HomologyFinder, self).__init__(S, T, **kw)
        self._n_by_d_ = self.seed_count_by_d_()

    def score_num_seeds(self, **kw):
        """Calculates our key central statistics based on m-dependent CLT. For
        a given observation of number of seeds in a region of interest (ROI),
        calculates the z-core against the H0 (unrelated) and H1 (related)
        models.

        Keyword Args:
            num_seeds (int):
                Number of observed seed in the ROI.
            area (int|float):
                Area of the ROI.
            seglen (int):
                Homologous segment length in H1 model.

        Returns:
            tuple (float): z-scores in H0 and H1 models
        """
        num_seeds = kw['num_seeds']
        A = kw['area']
        K = kw['seglen']

        if A == 0:
            return float('-inf'), float('-inf')

        p_neg = 1. / len(self.alphabet)
        pw_neg = p_neg ** self.wordlen
        p_pos = (1 - self.gap_prob) * (1 - self.subst_prob)
        pw_pos = p_pos ** self.wordlen

        mu_neg = A * pw_neg
        mu_pos = mu_neg + K * pw_pos
        sd_neg = np.sqrt(A * (
            (1 - pw_neg) * (pw_neg + 2 * p_neg * pw_neg / (1 - p_neg)) -
            2 * self.wordlen * pw_neg ** 2
        ))
        sd_pos = np.sqrt(sd_neg**2 + K * (
            (1 - pw_pos) * (pw_pos + 2 * p_pos * pw_pos / (1 - p_pos)) -
            2 * self.wordlen * pw_pos ** 2
        ))
        s0 = (num_seeds - mu_neg) / sd_neg  # score under H0
        s1 = (num_seeds - mu_pos) / sd_pos  # score under H1
        return s0, s1

    def band_radius(self, K):
        """Wraps :func:`band_radius` with our mutation parameters and sequence
        lengths.

        Args:
            K (int): expected alignment length of interest.

        Returns:
            int: radius of band for desired :attr:`sensitivity`.
        """
        return band_radius(K, self.gap_prob, self.sensitivity)

    def band_radii(self, Ks):
        """Wraps :func:`band_radii` with our mutation parameters and sequence
        lengths.

        Args:
            K (int): expected alignment lengths of interest.

        Returns:
            int: radius of band for desired :attr:`sensitivity`.
        """
        return band_radii(Ks, self.gap_prob, self.sensitivity)

    def score_diagonal_bands(self, seglens):
        """Scores all diagonal bands looking for specified segment lengths
        against both H0 and H1.

        Args:
            seglens (list|int):
                (minimum) segment lengths of interest for each diagonal
                position or constant number.
            d_center (int):
                center of diagonal band.
            d_radius (int):
                radius of diagonal band.

        Returns:
            np.array: H0 and H1 scores of each diagonal band.
        """
        try:
            iter(seglens)
            radii = self.band_radii(seglens)
        except TypeError:
            radius = self.band_radius(seglens)
            seglens = np.ones(len(self._n_by_d_)) * seglens
            radii = np.ones(len(self._n_by_d_)) * radius

        assert len(seglens) == len(radii) == len(self._n_by_d_)
        n_by_d_cum = np.cumsum(self._n_by_d_)
        # Each band gets two scores, one for H0 and one for H1
        scores_by_d_ = np.zeros((len(n_by_d_cum), 2))
        for d_ in range(len(self._n_by_d_)):
            seglen = seglens[d_]
            d = d_ - self.d0
            # don't score things too close to the edges, the small area
            # leads to unnecessarily high scores.
            if d_ < seglen or d > len(self.S) - seglen:
                scores_by_d_[d_][0] = float('-inf')
                scores_by_d_[d_][1] = float('-inf')
                continue

            radius = radii[d_]
            m = int(max(0, d_ - radius))
            M = int(min(len(self._n_by_d_) - 1, d_ + radius))
            n_in_band = n_by_d_cum[M] - n_by_d_cum[m]

            L = wall_to_wall_distance(len(self.S), len(self.T), d)
            area = L * radius * 2
            s0, s1 = self.score_num_seeds(num_seeds=n_in_band, area=area,
                                          seglen=seglen)
            scores_by_d_[d_][0], scores_by_d_[d_][1] = s0, s1
        return scores_by_d_

    def score_anti_bands(self, seglen, d_center=None, d_radius=None):
        """Scores antidiagonal bands within a given diagonal band for both H0
        and H1.

        Args:
            seglen (int):
                (minimum) segment length of interest.
            d_center (int):
                center of diagonal band.
            d_radius (int):
                radius of diagonal band.

        Returns:
            np.array: scores of each antidiagonal position in the band.
        """
        d_min, d_max = d_center - d_radius, d_center + d_radius
        area = (d_max - d_min) * seglen

        # X[a] = number of seeds in (d_min, d_max) at anti position a
        n_by_a = self.seed_count_by_a(d_center, d_radius)
        n_by_a_cum = np.cumsum(n_by_a)

        # X[a] = (s0, s1) the scores for ROI
        scores_by_a = -1e4 * np.ones((len(n_by_a_cum), 2))
        for anti in range(len(n_by_a_cum)):
            a_radius = int(np.ceil(seglen/2))
            m = max(0, anti - a_radius)
            M = min(len(n_by_a_cum) - 1, anti + a_radius)
            n_in_band = n_by_a_cum[M] - n_by_a_cum[m]
            s0, s1 = self.score_num_seeds(num_seeds=n_in_band, area=area,
                                          seglen=seglen)
            scores_by_a[anti][0], scores_by_a[anti][1] = s0, s1
        return scores_by_a

    def similar_segments(self, min_seglen, mode='H1'):
        """A sub-quadratic Local similarity search algorithm that finds
        diagonal regions of similarity between given sequences.

        Args:
            min_seglen (int):
                minimum length of similar segments to look for.
            mode (str):
                either 'H0' or 'H1' specifying the null hypothesis.

        Yields:
            tuple: coordinates of similar region in diagonal coordinates
            ``((d_min, d_max), (a_min, a_max))``.
        """
        self.log('finding local homologies between %s and %s' %
                 (self.S.content_id[:8], self.T.content_id[:8]))
        assert min_seglen > 0
        threshold = {'H0': 2, 'H1': -1.5}[mode]
        key = {'H0': 0, 'H1': 1}[mode]
        K = min_seglen
        d_radius = int(np.ceil(self.band_radius(K)))
        a_radius = int(np.ceil(K / 2))

        scores_by_d_ = self.score_diagonal_bands(K)
        d_peaks = find_peaks(scores_by_d_[:, key], d_radius, threshold)
        for d_min_, d_max_ in d_peaks:
            d_min, d_max = d_min_ - self.d0, d_max_ - self.d0
            center, radius = (d_max + d_min) / 2, (d_max - d_min) / 2
            scores_by_a = self.score_anti_bands(min_seglen, d_center=center,
                                                d_radius=radius)
            a_peaks = find_peaks(scores_by_a[:, key], a_radius, threshold)

            for (a_min, a_max) in a_peaks:
                yield (d_min, d_max), (a_min, a_max)
