# -*- coding: utf-8 -*-
"""
.. wikisection:: overview
    :title: (6) Alignment-free local similarity detection with Word-Blot

    The :mod:`biseqt.blot` module implements the WordBlot algorithm.

    >>> from biseqt.blot import WordBlot
    >>> from biseqt.sequence import Sequence, Alphabet
    >>> from biseqt.stochastics import rand_seq, MutationProcess
    >>> A = Alphabet('ACGT')
    >>> S = rand_seq(A, 100)
    >>> M = MutationProcess(A, go_prob=.1, ge_prob=.1, subst_probs=.2)
    >>> T, _ = M.mutate(S)
    >>> WB = WordBlot(S, T, wordlen=5, alphabet=A, g_max=.5, \
                            sensitivity=.99)
    >>> list(WB.similar_segments(50, .7))
    [(((-28, 7), (0, 99)), 1.1988215486845972, 0.7630903108462143)]
    >>>
"""
import sys
import warnings
import numpy as np
import logging
from scipy.special import erfcinv
from scipy.spatial import cKDTree
from .seeds import SeedIndex, SeedIndexMultiple
from .kmers import as_kmer_seq
from .util import Logger


# so we can catch numpy warnings
warnings.filterwarnings('error')


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
        r = \\mathrm{erf}^{-1}\\left(1-\\epsilon\\right)
            \\sqrt{2gK}

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
    C = erfcinv(epsilon) * np.sqrt(2 * gap_prob)
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
    """The mean and standrad deviation of the limiting normal distribution
    under the :math:`H_0` (unrelated) model given by:

    .. math::
        \\begin{aligned}
            \mu_0 & = Ap^w \\\\
            \sigma_0^2 & =
                A(1 - p^w)\left(p^w + \\frac{2p^{w+1}}{1 - p}\\right)
                - 2wp^{2w}
        \\end{aligned}

    where :math:`w` is the word length, :math:`A` is the area of the ROI, and
    :math:`p = \\frac{1}{|\Sigma|}` with :math:`|\Sigma|` being the alphabet
    length.
    """
    p_H0 = 1. / alphabet_len
    pw_H0 = p_H0 ** wordlen

    mu_H0 = area * pw_H0
    sd_H0 = np.sqrt(area * (
        (1 - pw_H0) * (pw_H0 + 2 * p_H0 * pw_H0 / (1 - p_H0)) -
        2 * wordlen * pw_H0 ** 2
    ))
    return mu_H0, sd_H0


def H1_moments(alphabet_len, wordlen, area, seglen, p_match):
    """The mean and standrad deviation of the limiting normal distribution under
    the :math:`H_1` (related) model given by:

    .. math::
        \\begin{aligned}
            \mu_1 & = \mu_0 + Kp^w \\\\
            \sigma_1^2 & = \sigma_0^2
                + K(1 - p^w)\left(p^w + \\frac{2p^{w+1}}{1 - p}\\right)
                - 2wp^{2w}
        \\end{aligned}

    where :math:`w` is the word length, :math:`K` is the similarity length, and
    :math:`p` is the match probability.
    """
    mu_H0, sd_H0 = H0_moments(alphabet_len, wordlen, area)

    p_H1 = p_match
    if p_H1 == 1.:
        # we can't let p_H1 == 1 because we get division by zero below
        p_H1 = 1 - np.finfo(float).eps
    pw_H1 = p_H1 ** wordlen

    mu_H1 = mu_H0 + seglen * pw_H1
    sd_H1 = np.sqrt(sd_H0 ** 2 + seglen * (
        (1 - pw_H1) * (pw_H1 + 2 * p_H1 * pw_H1 / (1 - p_H1)) -
        2 * wordlen * pw_H1 ** 2
    ))
    return mu_H1, sd_H1


# FIXME the fact that we have organized our data as self.S and self.T is the
# main blocker for merging pairwise and multiple sequence implementations.
# Also involved: SeedIndexMultiple
class WordBlot(SeedIndex):
    """A similarity finder based on m-dependent CLT statistics.

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
        super(WordBlot, self).__init__(S, T, **kw)

    def score_num_seeds(self, **kw):
        """Calculates our key central statistics based on m-dependent CLT. For
        a given observation of number of seeds in a region of interest (ROI),
        calculates the z-core against the :math:`H_0` (unrelated) and
        :math:`H_1` (related) models. In either case, the distribution of
        number of seeds, being a sum of m-dependent identically distributed
        random variables, is asymptotically normal with parameters given by
        :func:`H0_moments`, :func:`H1_moments`.

        Keyword Args:
            num_seeds (int):
                Number of observed seed in the ROI.
            area (int|float):
                Area of the ROI.
            seglen (int):
                Similar segment length in H1 model.
            p_match (float):
                Expected match probability at any given position.

        Returns:
            tuple (float): z-scores in H0 and H1 models
        """
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

    def segment_dims(self, d_band=None, a_band=None):
        """Calculate the edit path length :math:`K` and the area of the given
        diagonal/antiodiagonal segment.

        Keyword Args:
            d_band (tuple): lower and upper diagonals limiting the segment.
            a_band (tuple): lower and upper antidiagonal positions limiting the
                segment.

        Returns:
            tuple: edit path length and segment area.
        """
        a_min, a_max = a_band
        d_min, d_max = d_band
        K = (a_max - a_min) / 2
        A = (d_max - d_min) * K
        return K, A

    def estimate_match_probability(self, num_seeds, d_band=None, a_band=None):
        """Estimate the edit path match probability given the provided observed
        number of seeds :math:`n` in given sigment:

        .. math::
            \hat{p} = \\left(\\frac{n - Ap_\circ^w}{K}\\right)^{\\frac{1}{w}}

        where :math:`n, A, p_\circ, K, w` are the number of seeds, segment
        area, null probability of matching nucleotides, segment length, and the
        word length, respectively.

        Args:
            num_seeds (int): number of seeds observed in segment.

        Keywords Args:
            d_band (tuple):
                lower and upper diagonals limiting the segment.
            a_band (tuple):
                lower and upper antidiagonal positions limiting the segment.

        Returns:
            float: estimated match probability
        """
        K, area = self.segment_dims(d_band=d_band, a_band=a_band)
        word_p_null = (1./len(self.alphabet)) ** self.wordlen
        # NOTE K effectively has become the projected alignment length not the
        # full length. This is REALLY important in interpretation but I think
        # must things are currently consistent (TODO full review needed).
        word_p = (num_seeds - area * word_p_null) / K
        try:
            match_p = np.exp(np.log(word_p) / self.wordlen)
        except Warning:
            # presumably this happened because word_p was too small for log
            match_p = 0
        return min(match_p, 1)

    @classmethod
    def find_all_neighbors(cls, seeds, d_radius, a_radius):
        """For each seed finds all seeds in its neighborhood defined by:

        .. math::

            U_{(d, a)} = \\{(d', a'): |d - d'| < r_d,  |a - a'| < r_a \\}

        This is done using a Quad-Tree in :math:`O(m lg m)` time where m is the
        number of seeds.

        Returns:
            list: tuples ``((d, a), neighs)`` where ``neighs`` is a list of
                  neighbor indices.
        """
        # normalize the two diameters so we can use a standard L∞ neighborhood.
        # typically a_diam is larger, so scale up d values proportionally
        d_coeff = 1. * a_radius / d_radius
        radius = a_radius

        all_seeds = list(cls.to_diagonal_coordinates(i, j) for i, j in seeds)
        if not all_seeds:
            return []
        all_seeds_scaled = np.array([(d * d_coeff, a) for d, a in all_seeds])
        quad_tree = cKDTree(all_seeds_scaled)
        all_neighs = quad_tree.query_ball_tree(quad_tree, radius,
                                               p=float('inf'))
        # all_neighs[i] is the indices of the neighbors of all_seeds[i]; this
        # always contains the seed itself (i.e always: i in neighs[i])
        for idx, _ in enumerate(all_neighs):
            all_neighs[idx].remove(idx)
        return zip(all_seeds, all_neighs)

    def score_seeds(self, K):
        """Find the neighbors of each seed in the sense of
        :func:`find_all_neighbors` and estimates the match probability of the
        segment centered at the coordinates of each seed.

        Args:
            K (int): the similarity legnth of interest that dictates
                neighborhood shapes and estimated match probabilities.

        Returns:
            list(dict): List of dictionaries with keys: ``seed`` (coordinates
            of exactly matching kmer in diagonal coordinates), ``neighs`` (list
            of indices of neighbors of this seed in the appropriate diagonal
            strip), ``p`` the estimated match probability
            of a segment centered at the seed.
        """
        d_radius = int(np.ceil(self.band_radius(K)))
        a_radius = K

        seeds_with_neighs = self.find_all_neighbors(
            self.seeds(exclude_trivial=True), d_radius, a_radius
        )

        def _p(d, a, n):
            d_band = (d - d_radius, d + d_radius)
            a_band = (a - a_radius, a + a_radius)
            return self.estimate_match_probability(n, d_band=d_band,
                                                   a_band=a_band)

        return [{'seed': (d, a), 'neighs': neighs, 'p': _p(d, a, len(neighs))}
                for (d, a), neighs in seeds_with_neighs]

    def similar_segments(self, K_min, p_min, score=True, at_least_one=False):
        """Find all maximal local similarities of given minium length and match
        probability. Additionally for each segment, the match probability is
        estimated and H0/H1 scores are calculated.

        Args:
            K_min (int):
                minimum required length of similarity.
            p_min (float):
                Minimum required match probability at each position.
            score (bool):
                Whether to score the segment against H1 by counting seeds
                inside.

        Yields:
            dict: dictionary with keys: ``segment`` (coordinates of similar
            region in diagonal coordinates ``((d_min, d_max), (a_min,
            a_max))``)), ``p`` the estimated match probability, and ``score``
            the H1 z-score if keyword argument ``score`` is true.
        """
        self.log('finding local similarities between %s and %s' %
                 (self.S.content_id[:8], self.T.content_id[:8]))
        d_radius = int(np.ceil(self.band_radius(K_min)))
        a_radius = K_min
        scored_seeds = self.score_seeds(K_min)

        def _update_seg(seg, seed):
            d, a = seed
            if seg is None:
                d_min, d_max = d - d_radius, d + d_radius
                a_min, a_max = a - a_radius, a + a_radius
            else:
                (d_min, d_max), (a_min, a_max) = seg
                d_min = min(d - d_radius, d_min)
                d_max = max(d + d_radius, d_max)
                a_min = min(a - a_radius, a_min)
                a_max = max(a + a_radius, a_max)
            return (d_min, d_max), (a_min, a_max)

        avail = [rec['p'] >= p_min for rec in scored_seeds]
        if not any(avail) and at_least_one:
            # we're obliged to return something, let the highest probability
            # seed go through.
            assert len(scored_seeds), 'no seeds found while at_least_one=True'
            avail[np.argmax([rec['p'] for rec in scored_seeds])] = True
        while True:
            try:
                seed_idx = avail.index(True)
            except ValueError:
                break
            stack = [seed_idx]
            avail[seed_idx] = False
            ps_in_seg = [scored_seeds[seed_idx]['p']]
            seg = None
            while stack:
                idx = stack.pop()
                ps_in_seg.append(scored_seeds[idx]['p'])
                seg = _update_seg(seg, scored_seeds[idx]['seed'])
                for neigh in scored_seeds[idx]['neighs']:
                    if avail[neigh]:
                        stack.append(neigh)
                        avail[neigh] = False
            if seg is None:
                break
            else:
                (d_min, d_max), (a_min, a_max) = seg
                d_min = min(len(self.S), max(d_min, -len(self.T)))
                d_max = min(len(self.S), max(d_max, -len(self.T)))
                a_min = max(a_min + a_radius, 0)
                a_max = min(a_max - a_radius, len(self.S) + len(self.T))
                seg = (d_min, d_max), (a_min, a_max)
            # NOTE the following is more justifiable but it matches the
            # average. TODO turn this in into an experiment to justify
            # p_hat = self.estimate_match_probability(
            #   n, d_band=seg[0], a_band=seg[1])
            p_hat = sum(ps_in_seg) / len(ps_in_seg)
            res = {'segment': seg, 'p': p_hat}
            if score:
                n = self.seed_count(d_band=seg[0], a_band=seg[1])
                K_hat, area_hat = self.segment_dims(d_band=seg[0],
                                                    a_band=seg[1])
                scores = self.score_num_seeds(num_seeds=n, area=area_hat,
                                              seglen=K_hat, p_match=p_hat)
                res['scores'] = scores
            yield res


class WordBlotOverlap(WordBlot):
    """A specialized version of WordBlot for detecting overlap
    (suffix-prefix) similarities between sequences (e.g. in a sequencing
    context)."""
    # FIXME get rid of p_min argument it's only affecting our alignment length
    # estimates. If too worried about g_max being wild, can use wall to wall
    # instead of expected_overlap_len!
    def score_seeds(self, p_min):
        """For each seed finds all seeds in its neighborhood defined by:

        .. math::

            U_{(d, a)} = \\{(d', a'): |d - d'| < r_d(d),  |a - a'| < r_a(d) \\}

        in such a way that each seed's neighborhood is the entire diagonal band
        containing it. Each seed recieves an estimated match probability for a
        similarity on its diagonal band and a z-score with respect to the H1
        model.

        Returns:
            list: dicts with keys ``seed`` (diagonal coordinates of seed),
                  ``p`` (estimated match probability of overlap alignment),
                  ``L`` (the length of an alignment through the seed's band
                  according to sequence lengths), ``score`` (the z-score with
                  respect to H1), and ``r`` (band radius at the seed's
                  coordinates).
        """
        # use 1 - p_min as an overestimate of gap prob
        gap_prob = 1 - p_min

        def _len(d):
            return expected_overlap_len(len(self.S), len(self.T), d, gap_prob)

        def _rad(d):
            return np.ceil(self.band_radius(_len(d)))

        all_seeds = list(self.to_diagonal_coordinates(i, j)
                         for i, j in self.seeds(exclude_trivial=True))
        if not all_seeds:
            return []
        all_seeds_scaled = np.array([(d / _rad(d), a / _len(d))
                                     for d, a in all_seeds])
        quad_tree = cKDTree(all_seeds_scaled)
        all_neighs = quad_tree.query_ball_tree(quad_tree, 1, p=float('inf'))
        # all_neighs[i] is the indices of the neighbors of all_seeds[i]; this
        # always contains the seed itself (i.e always: i in neighs[i])
        for idx, _ in enumerate(all_neighs):
            all_neighs[idx].remove(idx)
        seeds_with_neighs = zip(all_seeds, all_neighs)

        def _p(d, n):
            L = _len(d)
            d_radius = int(np.ceil(self.band_radius(L)))
            area = 2 * d_radius * L
            word_p_null = (1./len(self.alphabet)) ** self.wordlen
            word_p = (n - area * word_p_null) / L
            try:
                match_p = np.exp(np.log(word_p) / self.wordlen)
            except Warning:
                # presumably this happened because word_p was too small for log
                match_p = 0
            return min(match_p, 1)

        return [{'seed': (d, a), 'r': _rad(d), 'L': _len(d),
                 'p': _p(d, len(neighs))}
                for (d, a), neighs in seeds_with_neighs]

    def highest_scoring_overlap_band(self, p_min, score=True):
        """Finds the highest scoring diagonal band according to probabiliy
        estimations of :func:`score_seeds`.

        Returns:
            dict: with same keys ``p, score, d_band`` consistent with the
                  output of :func:`score_seeds` for the highest scoring
                  diagonal band.
        """
        scored_seeds = self.score_seeds(p_min)
        if not scored_seeds:
            return None
        idx = max(range(len(scored_seeds)), key=lambda i: scored_seeds[i]['p'])
        seed, rad = scored_seeds[idx]['seed'], scored_seeds[idx]['r']
        p_hat, overlap_len = scored_seeds[idx]['p'], scored_seeds[idx]['L']
        d_band = seed[0] - rad, seed[0] + rad
        res = {'d_band': d_band, 'p': p_hat, 'len': overlap_len}
        if score:
            area = 2 * rad * overlap_len
            mu_H1, sd_H1 = H1_moments(len(self.alphabet), self.wordlen, area,
                                      overlap_len, p_hat)
            num_seeds = self.seed_count(d_band=d_band)
            z_H1 = (num_seeds - mu_H1) / sd_H1
            res['score'] = z_H1
        return res


class WordBlotOverlapRef(WordBlotOverlap):
    """An in-memory, SQL-free version of :class:`WordBlotOverlap` for faster
    comparisons.  Due to implementation details the word length is constrained
    above by the available memory.

    Attributes:
        allowed_memory (int|float): allocatable memory in GB for kmers index.
    """
    def __init__(self, ref, allowed_memory=1, **kw):
        self.wordlen = kw['wordlen']
        self.alphabet = kw['alphabet']
        self.g_max = kw['g_max']
        self.sensitivity = kw['sensitivity']
        self.log_level = kw.get('log_level', logging.INFO)
        self.S = ref
        num_kmers = len(self.alphabet) ** self.wordlen
        mem_needed = sys.getsizeof(num_kmers) * num_kmers
        mem_needed_gb = np.power(2, np.log2(mem_needed) - 30)
        assert allowed_memory > 0, 'allowed memory must be positive'
        self.allowed_memory = allowed_memory
        if mem_needed_gb > self.allowed_memory:
            msg = 'not enough memory (max = %.2f GB) ' % allowed_memory
            msg += 'to store %d-mers ' % self.wordlen
            msg += '(%.2f GB needed)' % mem_needed_gb
            raise MemoryError(msg)
        self.kmer_hits = [[] for _ in range(num_kmers)]
        for pos, kmer in enumerate(as_kmer_seq(ref, self.wordlen)):
            self.kmer_hits[kmer].append(pos)
        self.T = None
        relpath = 'python-object'
        log_header = '%d-mer cache (%s)' % (self.wordlen, relpath)
        self._logger = Logger(log_level=self.log_level, header=log_header)

    def seeds(self, exclude_trivial=True):
        assert self.T is not None
        for pos, kmer in enumerate(as_kmer_seq(self.T, self.wordlen)):
            for pos_ref in self.kmer_hits[kmer]:
                if self.S == self.T and exclude_trivial and pos == pos_ref:
                    continue
                yield pos_ref, pos

    def score_seeds_(self, seq, p_min):
        self.T = seq
        return super(WordBlotOverlapRef, self).score_seeds(p_min)

    def highest_scoring_overlap_band(self, seq, p_min):
        self.T = seq
        return super(WordBlotOverlapRef, self).highest_scoring_overlap_band(
            p_min, score=False
        )


class WordBlotLocalRef(WordBlot):
    """An in-memory, SQL-free version of :class:`WordBlot` for faster
    comparisons. Due to implementation details the word length is constrained
    above by the available memory.

    Attributes:
        allowed_memory (int|float): allocatable memory in GB for kmers index.
    """
    def __init__(self, ref, allowed_memory=1, **kw):
        self.wordlen = kw['wordlen']
        self.alphabet = kw['alphabet']
        self.g_max = kw['g_max']
        self.sensitivity = kw['sensitivity']
        self.log_level = kw.get('log_level', logging.INFO)
        self.S = ref
        num_kmers = len(self.alphabet) ** self.wordlen
        mem_needed = sys.getsizeof(num_kmers) * num_kmers
        mem_needed_gb = np.power(2, np.log2(mem_needed) - 30)
        assert allowed_memory > 0, 'allowed memory must be positive'
        self.allowed_memory = allowed_memory
        if mem_needed_gb > self.allowed_memory:
            msg = 'not enough memory (max = %.2f GB) ' % allowed_memory
            msg += 'to store %d-mers ' % self.wordlen
            msg += '(%.2f GB needed)' % mem_needed_gb
            raise MemoryError(msg)
        self.kmer_hits = [[] for _ in range(num_kmers)]
        for pos, kmer in enumerate(as_kmer_seq(ref, self.wordlen)):
            self.kmer_hits[kmer].append(pos)
        self.T = None
        relpath = 'python-object'
        log_header = '%d-mer cache (%s)' % (self.wordlen, relpath)
        self._logger = Logger(log_level=self.log_level, header=log_header)

    def seeds(self, exclude_trivial=True):
        assert self.T is not None
        for pos, kmer in enumerate(as_kmer_seq(self.T, self.wordlen)):
            for pos_ref in self.kmer_hits[kmer]:
                if self.S == self.T and exclude_trivial and pos == pos_ref:
                    continue
                yield pos_ref, pos

    def score_seeds_(self, seq, K):
        self.T = seq
        return super(WordBlotLocalRef, self).score_seeds(K)

    def similar_segments(self, seq, K_min, p_min, at_least_one=False):
        self.T = seq
        kw = {'score': False, 'at_least_one': at_least_one}
        for res in super(WordBlotLocalRef, self).similar_segments(K_min, p_min,
                                                                  **kw):
            yield res


# FIXME lots of code duplication here, can it be cleaned up? (cf. note above
# WordBlot; the main issue is self.S and self.T being used everywhere there,
# and presumably lots of hidden pairwise assumptions). The good news is:
# every function implemented below should in principle work exactly as is for
# pairwise.
class WordBlotMultiple(SeedIndexMultiple):
    """A multiple sequence similarity finder based on m-dependent CLT
    statistics.

    Attributes:
        g_max (float):
            Upper bound for indel probabilities in mutation model.
        sensitivity (float):
            Desired sensitivity of bands.
    """
    def __init__(self, *seqs, **kw):
        g_max, sensitivity = kw.pop('g_max'), kw.pop('sensitivity')
        assert 0 < g_max < 1 and 0 < sensitivity < 1
        self.g_max = g_max
        self.sensitivity = sensitivity
        super(WordBlotMultiple, self).__init__(*seqs, **kw)

    def band_radius(self, K):
        """Wraps :func:`band_radius` with our mutation parameters and sequence
        lengths.

        Args:
            K (int): expected alignment length of interest.

        Returns:
            int: radius of band for desired :attr:`sensitivity`.
        """
        return band_radius(K, self.g_max, self.sensitivity)

    def wall_to_wall_distance(self, ds):
        """Wall to wall distance :math:`L` for a diagonal position :math:`d_1,
        \ldots, d_{n-1}`.

        The distance :math:`L` is the largest possible value of the
        antidiagonal coordinate once all diagonal coordinates are fixed.
        """
        raise NotImplementedError

    def estimate_match_probability(self, num_seeds, K, volume):
        """Estimate the edit path match probability given the provided observed
        number of seeds in given sigment.

        .. math::
            \hat{p} = \\left(
                            \\frac{n - Vp_\circ^{w(N-1)}}{K}
                      \\right)^{\\frac{1}{w(N-1)}}

        where :math:`n, N, V, p_\circ, w` are the number of seeds, number of
        sequences, segment n-d volume, null probability of matching
        nucleotides, and word length, respectively.

        Args:
            num_seeds (int|float): number of seeds observed in segment.
            K (int|float): the expected length of the alignment.
            volume (int|float): the n-d volume of the region of interest.

        Returns:
            float: estimated match probability
        """
        power = self.wordlen * (len(self.seqs) - 1)
        word_p_null = (1. / len(self.alphabet)) ** power

        if num_seeds > 0:
            word_p = (num_seeds - volume * word_p_null) / K
            try:
                match_p = np.exp(np.log(word_p) / self.wordlen)
            except Warning:
                # presumably this happened because word_p was too small for log
                match_p = 0
        else:
            match_p = 0
        return min(match_p, 1)

    def score_seeds(self, K):
        """Find the neighbors of each seed in the sense of
        :func:`find_all_neighbors` and estimates the match probability of the
        segment centered at the coordinates of each seed.

        Args:
            K (int): the similarity legnth of interest that dictates
                neighborhood shapes and estimated match probabilities.

        Returns:
            list: list of tuples ``((d, a), neighs, p)`` where ``(d, a)``
                  is the diagonal coordinates of the seed, ``neighs`` is a list
                  of integer indices of seeds in its diagonal/antidiagonal
                  neighborhood, and ``p`` is the estimated match probability of
                  a segment of length ``K`` centered at the seed.
        """
        d_radius = int(np.ceil(self.band_radius(K)))
        a_radius = int(np.ceil(len(self.seqs) * K / 2.))
        seeds_with_neighs = self.find_all_neighbors(d_radius, a_radius)
        volume = (2 * d_radius) ** (len(self.seqs) - 1) * 2 * a_radius

        def _p(n):
            return self.estimate_match_probability(n, K, volume)

        return [{'seed': (ds, a), 'neighs': neighs, 'p': _p(len(neighs))}
                for (ds, a), neighs in seeds_with_neighs]

    def find_all_neighbors(self, d_radius, a_radius):
        """For each seed finds all seeds in its neighborhood defined by:

        .. math::

            U_{(d_1,\ldots,d_{n-1}, a)} =
                \\{(d_1', \ldots, d_{n-1}', a'):
                    |d_k - d_k'| < r_d,  |a - a'| < r_a \\}

        This is done using a kD-Tree in :math:`O(nm lg m)` time where m is the
        number of seeds and n is the number of sequences.

        Returns:
            list: tuples ``((ds, a), neighs)`` where ``neighs`` is a list of
                  neighbor indices.
        """
        self.log('finding all neighbors using a kD-tree')
        # normalize the two diameters so we can use a standard L∞ neighborhood.
        # typically a_diam is larger, so scale up d values proportionally
        d_coeff = 1. * a_radius / d_radius
        radius = a_radius

        all_seeds = list(self.seeds())
        if not all_seeds:
            return []
        all_seeds_scaled = np.array([[d * d_coeff for d in ds] + [a]
                                     for ds, a in all_seeds])
        quad_tree = cKDTree(all_seeds_scaled)
        all_neighs = quad_tree.query_ball_tree(quad_tree, radius,
                                               p=float('inf'))
        # all_neighs[i] is the indices of the neighbors of all_seeds[i]; this
        # always contains the seed itself (i.e always: i in neighs[i])
        for idx, _ in enumerate(all_neighs):
            all_neighs[idx].remove(idx)
        return zip(all_seeds, all_neighs)

    def score_num_seeds(self, num_seeds, **kw):
        """The exrtension of :func:`WordBlot.score_num_seeds` to multiple
        sequence comparisons.

        Keyword Args:
            num_seeds (int):
                Number of observed seed in the ROI.
            volume (int|float):
                Area of the ROI.
            seglen (int):
                Similar segment length in H1 model.
            p_match (float):
                Expected match probability at any given position.

        Returns:
            tuple (float): z-scores in H0 and H1 models
        """
        p_match = kw['p_match']
        seglen = kw['seglen']
        volume = kw['volume']
        if volume == 0:
            return float('-inf'), float('-inf')
        if p_match == 1.:
            # we can't let p_H1 == 1 because we get division by zero below
            p_match = 1 - np.finfo(float).eps

        p_H0 = (1. / len(self.alphabet)) ** (len(self.seqs) - 1)
        pw_H0 = p_H0 ** self.wordlen

        mu_H0 = volume * pw_H0
        sd_H0 = np.sqrt(volume * (
            (1 - pw_H0) * (pw_H0 + 2 * p_H0 * pw_H0 / (1 - p_H0)) -
            2 * self.wordlen * pw_H0 ** 2
        ))

        p_H1 = p_match
        pw_H1 = p_H1 ** self.wordlen

        mu_H1 = mu_H0 + seglen * pw_H1
        sd_H1 = np.sqrt(sd_H0 ** 2 + seglen * (
            (1 - pw_H1) * (pw_H1 + 2 * p_H1 * pw_H1 / (1 - p_H1)) -
            2 * self.wordlen * pw_H1 ** 2
        ))

        z_H0 = (num_seeds - mu_H0) / sd_H0  # score under H0
        z_H1 = (num_seeds - mu_H1) / sd_H1  # score under H1
        return z_H0, z_H1

    def similar_segments(self, K_min, p_min, score=True, at_least_one=False):
        """Find all maximal local similarities of given minium length and match
        probability. Additionally for each segment, the match probability is
        estimated and H0/H1 scores are calculated.

        Args:
            K_min (int):
                minimum required length of similarity.
            p_min (float):
                Minimum required match probability at each position.
            score (bool):
                Whether to score the segment against H1 by counting seeds
                inside.

        Yields:
            dict: dictionary with keys: ``segment`` (coordinates of similar
            region in diagonal coordinates, ``p`` the estimated match
            probability, and ``score`` the H1 z-score if keyword argument
            ``score`` is true.
        """
        d_radius = int(np.ceil(self.band_radius(K_min)))
        a_radius = int(np.ceil(len(self.seqs) * K_min / 2.))
        scored_seeds = self.score_seeds(K_min)
        self.log('finding local similarities between %d sequences' %
                 len(self.seqs))

        def _update_seg(seg, seed):
            ds, a = seed
            if seg is None:
                d_ranges = [None] * (len(self.seqs) - 1)
                for i in range(len(self.seqs) - 1):
                    d_ranges[i] = ds[i] - d_radius, ds[i] + d_radius
                a_range = a - a_radius, a + a_radius
            else:
                d_ranges, a_range = seg
                assert all(len(r) == 2 for r in d_ranges)  # pairs of min, max
                for i in range(len(self.seqs) - 1):
                    d_min, d_max = d_ranges[i]
                    d_ranges[i] = (min(ds[i], d_min),
                                   max(ds[i], d_max))
                a_min, a_max = seg[-1]
                a_range = (min(a - a_radius, a_min), max(a + a_radius, a_max))
            return d_ranges, a_range

        avail = [rec['p'] >= p_min for rec in scored_seeds]
        if not any(avail) and at_least_one:
            # we're obliged to return something, let the highest probability
            # seed go through.
            assert len(scored_seeds), 'no seeds found while at_least_one=True'
            avail[np.argmax([rec['p'] for rec in scored_seeds])] = True
        while True:
            try:
                seed_idx = avail.index(True)
            except ValueError:
                break
            stack = [seed_idx]
            avail[seed_idx] = False
            ps_in_seg = [scored_seeds[seed_idx]['p']]
            seg = None
            while stack:
                idx = stack.pop()
                ps_in_seg.append(scored_seeds[idx]['p'])
                seg = _update_seg(seg, scored_seeds[idx]['seed'])
                for neigh in scored_seeds[idx]['neighs']:
                    if avail[neigh]:
                        stack.append(neigh)
                        avail[neigh] = False
            if seg is None:
                break
            else:
                # clip overflowing values so translating back to standard
                # coordinates gives meaningful numbers.
                ds_range, a_range = seg
                for idx in range(len(self.seqs) - 1):
                    d_min, d_max = ds_range[idx]
                    d_min = min(
                        len(self.seqs[0]),
                        max(d_min, -len(self.seqs[idx + 1]))
                    )
                    d_max = min(
                        len(self.seqs[0]),
                        max(d_max, -len(self.seqs[idx + 1]))
                    )
                    ds_range[idx] = d_min, d_max
                a_min, a_max = a_range
                a_min = max(a_min, 0)
                a_max = min(a_max, sum(len(seq) for seq in self.seqs))
                seg = ds_range, (a_min, a_max)
            # NOTE the following is more justifiable but it matches the
            # average. TODO turn this in into an experiment to justify
            # p_hat = self.estimate_match_probability(
            #   n, d_band=seg[0], a_band=seg[1])
            p_hat = sum(ps_in_seg) / len(ps_in_seg)
            res = {'segment': seg, 'p': p_hat}
            if score:
                ds_band, a_band = seg
                n = self.seed_count(ds_band=ds_band, a_band=a_band)
                a_radius = int(np.ceil(len(self.seqs) * K_min / 2.))
                K_hat = np.ceil((a_band[1] - a_band[0]) / len(self.seqs))
                volume = a_band[1] - a_band[0]
                for d_band in ds_band:
                    volume *= d_band[1] - d_band[0]
                scores = self.score_num_seeds(num_seeds=n, volume=volume,
                                              seglen=K_hat, p_match=p_hat)
                res['scores'] = scores
            yield res
