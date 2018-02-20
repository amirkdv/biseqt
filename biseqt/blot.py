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
    C = 2 * erfcinv(epsilon) * np.sqrt(gap_prob)
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
    C = erfcinv(epsilon) * np.sqrt(2 * gap_prob)
    return np.array([max(1, int(np.ceil(C * np.sqrt(K))))
                     for K in expected_lens])


def H0_moments(alphabet_len, wordlen, area):
    p_H0 = 1. / alphabet_len
    pw_H0 = p_H0 ** wordlen

    # HACK for bio data, wordlen = 8
    # pw_H0 *= 2.5

    mu_H0 = area * pw_H0
    sd_H0 = np.sqrt(area * (
        (1 - pw_H0) * (pw_H0 + 2 * p_H0 * pw_H0 / (1 - p_H0)) -
        2 * wordlen * pw_H0 ** 2
    ))
    return mu_H0, sd_H0


def H1_moments(alphabet_len, wordlen, area, seglen, p_match):
    mu_H0, sd_H0 = H0_moments(alphabet_len, wordlen, area)

    p_H1 = p_match
    pw_H1 = p_H1 ** wordlen

    # HACK
    # pw_H1 *= 2.5

    mu_H1 = mu_H0 + seglen * pw_H1
    sd_H1 = np.sqrt(sd_H0 ** 2 + seglen * (
        (1 - pw_H1) * (pw_H1 + 2 * p_H1 * pw_H1 / (1 - p_H1)) -
        2 * wordlen * pw_H1 ** 2
    ))
    return mu_H1, sd_H1


class HomologyFinder(SeedIndex):
    """A homology finder based on m-dependent CLT statistics.

    Attributes:
        g_max (float):
            Upper bound for indel probabilities in mutation model.
        sensitivity (float):
            Desired sensitivity of bands.
    """
    def __init__(self, S, T, g_max=None, sensitivity=None, **kw):
        assert 0 < g_max < 1 and 0 < sensitivity < 1
        self.g_max = g_max
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
            p_match (float):
                Expected match probability at each position.

        Returns:
            tuple (float): z-scores in H0 and H1 models
        """
        # FIXME calculate all this stuff once, keep in object
        num_seeds = kw['num_seeds']
        area = kw['area']

        if area == 0:
            return float('-inf'), float('-inf')

        mu_H0, sd_H0 = H0_moments(len(self.alphabet), self.wordlen, area)
        mu_H1, sd_H1 = H1_moments(len(self.alphabet), self.wordlen, area,
                                  kw['seglen'], kw['p_match'])

        z_H0 = (num_seeds - mu_H0) / sd_H0  # score under H0
        z_H1 = (num_seeds - mu_H1) / sd_H1  # score under H1
        return z_H0, z_H1

    def band_radius(self, K):
        """Wraps :func:`band_radius` with our mutation parameters and sequence
        lengths.

        Args:
            K (int): expected alignment length of interest.

        Returns:
            int: radius of band for desired :attr:`sensitivity`.
        """
        return band_radius(K, self.g_max, self.sensitivity)

    def band_radii(self, Ks):
        """Wraps :func:`band_radii` with our mutation parameters and sequence
        lengths.

        Args:
            K (int): expected alignment lengths of interest.

        Returns:
            int: radius of band for desired :attr:`sensitivity`.
        """
        return band_radii(Ks, self.g_max, self.sensitivity)

    def score_diagonal_bands(self, Ks, p_match, edge_margin=100):
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
            edge_margin(int):
                The width of the diagonal margin ignored on either extreme.
                This is used to avoid division by small areas in calculating
                scores.

        Returns:
            np.array: H0 and H1 scores of each diagonal band.
        """
        try:
            iter(Ks)
            radii = self.band_radii(Ks)
        except TypeError:
            radius = self.band_radius(Ks)
            Ks = np.ones(len(self._n_by_d_)) * Ks
            radii = np.ones(len(self._n_by_d_)) * radius

        assert len(Ks) == len(radii) == len(self._n_by_d_)
        n_by_d_cum = np.cumsum(self._n_by_d_)
        # Each band gets two scores, one for H0 and one for H1
        scores_by_d_ = np.zeros((len(n_by_d_cum), 2))
        for d_ in range(len(self._n_by_d_)):
            K = Ks[d_]
            d = d_ - self.d0
            # don't score things too close to the edges, the small area
            # leads to unnecessarily high scores.
            if d_ < edge_margin or d > len(self.S) - edge_margin:
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
                                          seglen=K, p_match=p_match)
            scores_by_d_[d_][0], scores_by_d_[d_][1] = s0, s1
        return scores_by_d_

    def highest_scoring_overlap_band(self, p_min):
        """Finds the highest scoring diagonal band with respect to both H0 and
        H1 statistics with segment length calculated such that only overlap
        similarities are considered.

        Returns:
            tuple: containing the the highest scoring band and their score for
            each of H0 and H1 models.  ``((d0_min, d0_max), s0), ((d1_min,
            d1_max), s1)``.
        """
        self.log('scoring possible overlaps between %s and %s' %
                 (self.S.content_id[:8], self.T.content_id[:8]))

        # FIXME embed this directly in for loop below (only user), the name is
        # confusing (wrt indexing)
        ds = range(- len(self.T) + 1, len(self.S))
        # NOTE we are using 1 - p as an (over)estimate of g (g_max may be too
        # wild).
        Ks = [expected_overlap_len(len(self.S), len(self.T), d, 1 - p_min)
              for d in ds]
        scores_by_d_ = self.score_diagonal_bands(Ks, p_min)

        d_H0_ = np.argmax(scores_by_d_[:, 0])
        d_H1_ = np.argmax(scores_by_d_[:, 1])
        K_H0, K_H1 = Ks[d_H0_], Ks[d_H1_]
        s_H0, s_H1 = scores_by_d_[d_H0_, 0], scores_by_d_[d_H1_, 1]
        r_H0, r_H1 = self.band_radius(K_H0), self.band_radius(K_H1)

        d_H0, d_H1 = d_H0_ - self.d0, d_H1_ - self.d0
        seg_H0 = (d_H0 - r_H0, d_H0 + r_H0)
        seg_H1 = (d_H1 - r_H1, d_H1 + r_H1)
        return (seg_H0, s_H0), (seg_H1, s_H1)

    # returns area and alignment length
    def segment_dims(self, d_band=None, a_band=None):
        a_min, a_max = a_band
        d_min, d_max = d_band
        L = (a_max - a_min)
        # FIXME what to do with this?
        # K = (2. / (2 - self.gap_prob)) * L
        K = L
        area = K * (d_max - d_min)
        return K, area

    def estimate_match_probability(self, num_seeds, d_band=None, a_band=None):
        K, area = self.segment_dims(d_band=d_band, a_band=a_band)
        word_p_null = (1./len(self.alphabet)) ** self.wordlen
        word_p = (num_seeds - area * word_p_null) / K
        match_p = np.exp(np.log(word_p) / self.wordlen)
        return min(match_p, 1)

    def similar_segments(self, K_min, p_min, mode='H1'):
        """A sub-quadratic local similarity search algorithm that finds
        diagonal regions of similarity between given sequences.

        Args:
            K_min (int):
                minimum required length of homology.
            p_min (float):
                Minimum required match probability at each position.
            mode (str):
                either 'H0' or 'H1' specifying the null hypothesis.

        Yields:
            tuple: coordinates of similar region in diagonal coordinates, with
            a z-score, and estimated match probability:
            ``((d_min, d_max), (a_min, a_max)), z_score, match_prob``.
        """
        self.log('finding local homologies between %s and %s' %
                 (self.S.content_id[:8], self.T.content_id[:8]))
        assert K_min > 0
        threshold = {'H0': 10, 'H1': -1.5}[mode]
        key = {'H0': 0, 'H1': 1}[mode]
        d_radius = int(np.ceil(self.band_radius(K_min)))
        a_radius = int(np.ceil(K_min / 2))

        scores_by_d_ = self.score_diagonal_bands(K_min, p_min)
        d_peaks = find_peaks(scores_by_d_[:, key], d_radius, threshold)
        for d_min_, d_max_ in d_peaks:
            d_min, d_max = d_min_ - self.d0, d_max_ - self.d0
            area = (d_max - d_min) * K_min
            # n_by_a[a] = number of seeds in (d_min, d_max) at anti position a
            n_by_a = self.seed_count_by_a(d_min, d_max)
            n_by_a_cum = np.cumsum(n_by_a)
            scores_by_a = -1e4 * np.ones((len(n_by_a_cum), 2))
            for anti in range(len(n_by_a_cum)):
                m = max(0, anti - a_radius)
                M = min(len(n_by_a_cum) - 1, anti + a_radius)
                n_in_band = n_by_a_cum[M] - n_by_a_cum[m]
                # NOTE the segment in question at this point has
                # a_max - a_min = K
                s0, s1 = self.score_num_seeds(num_seeds=n_in_band, area=area,
                                              seglen=K_min, p_match=p_min)
                scores_by_a[anti][0], scores_by_a[anti][1] = s0, s1

            a_peaks = find_peaks(scores_by_a[:, key], a_radius, threshold)

            for (a_min, a_max) in a_peaks:
                segment = ((d_min, d_max), (a_min, a_max))
                num_seeds_in_segment = n_by_a_cum[a_max] - n_by_a_cum[a_min]
                match_p = self.estimate_match_probability(
                    num_seeds_in_segment, d_band=segment[0], a_band=segment[1])
                # NOTE we if we re-score the segment based on observed seglen
                # (i.e K_hat) we get scores that are eventually, indicative of
                # quality only and do not depend on segment length.
                #
                # z_score = np.max(scores_by_a[a_min:a_max, key])
                K_hat, area = self.segment_dims(d_band=(d_min, d_max),
                                                a_band=(a_min, a_max))
                z_score = self.score_num_seeds(num_seeds=num_seeds_in_segment,
                                               area=area, seglen=K_hat,
                                               p_match=p_min)[key]
                yield segment, z_score, match_p
