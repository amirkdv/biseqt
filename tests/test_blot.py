# -*- coding: utf-8 -*-
import pytest
import numpy as np

from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess
from biseqt.blot import find_peaks
from biseqt.blot import band_radius
from biseqt.blot import expected_overlap_len
from biseqt.blot import band_radii
from biseqt.blot import HomologyFinder


def test_find_peaks():
    threshold = 100
    radius = 3
    xs = [0, 1, 2, 3, 100, 5, 6, 7, 8]
    #              (   ^   )
    assert find_peaks(xs, radius, threshold) == [(3, 5)], \
        'single peak should be correctly detected'

    xs = [0, 1, 2, 3, 100, 5, 6, 7, 8, 9, 10, 11, 100, 13, 14, 15, 16]
    #              (   ^   )                  (    ^   )
    assert find_peaks(xs, radius, threshold) == [(3, 5), (11, 13)], \
        'two disjoint peaks should be correctly detected'

    xs = [0, 1, 2, 3, 100, 5, 6, 7, 100, 9, 10, 11, 12, 13]
    #              (   ^             ^   )
    assert find_peaks(xs, radius, threshold) == [(3, 9)], \
        'two overlapping peaks should be correctly merged'


def test_expected_overlap_len():
    n = 50  # sequence lengths
    gap = .1
    ds = range(-n, n + 1)
    aln_lens = np.array([expected_overlap_len(n, n, d, gap) for d in ds])
    assert np.all(np.diff(aln_lens[0: n]) >= 0), \
        'expected alignment length should increase as d goes from -|T| to 0'

    assert aln_lens[n] >= n, \
        'expected alignment length at diagonal 0 should be at least n'

    assert np.all(np.diff(aln_lens[n:]) <= 0), \
        'expected alignment length should increase as d goes from 0 to |S|'

    d = 0
    gaps = [i * .05 for i in range(6)]
    aln_lens = np.array([expected_overlap_len(n, n, d, g) for g in gaps])
    assert np.all(np.diff(aln_lens) >= 0), \
        'expected alignment length should increase with larger gap probability'

    d = 0
    N = [i * 100 for i in range(1, 5)]
    aln_lens = np.array([expected_overlap_len(n, N, d, g) for g in gaps])
    assert np.all(aln_lens) < 70, \
        'expected overlap length is controlled by the smaller sequence length'


def test_band_radius():
    gap_prob = .1
    sensitivity = 1 - 1e-3
    Ks = [i * 200 for i in range(1, 10)]
    Rs = [band_radius(K, gap_prob, sensitivity) for K in Ks]
    ratios = [Rs[i] / np.sqrt(Ks[i]) for i in range(len(Ks))]
    ratios_0_ctrd = np.array(ratios) - ratios[0] * np.ones(len(ratios))
    assert np.allclose(ratios_0_ctrd, np.zeros(len(ratios)), atol=1e-1), \
        'band radius must increase proportional to square root of alignment'

    n = 50
    sensitivity = .99
    gaps = [i * .05 for i in range(1, 7)]
    radii = np.array([band_radius(n, g, sensitivity) for g in gaps])
    assert np.all(np.diff(radii) >= 0), \
        'band radius should increase with increased gap probability'

    sensitivities = [1 - i * .05 for i in range(1, 7)]
    radii = np.array([band_radius(n, gap_prob, s)
                      for s in sensitivities])
    assert np.all(np.diff(radii) <= 0), \
        'band radius should decrease with decreased sensitivity'


def test_band_radii():
    n = 50  # sequence lengths
    g = .1
    sensitivity = .99
    radii = band_radii(range(n), gap_prob=g, sensitivity=sensitivity)

    assert len(radii) == n, 'correct number of radii returned'
    assert np.all(np.diff(radii) >= 0), \
        'radius must increase with expected length'


def test_overlap_band_radius():
    n = 50  # sequence lengths
    g = .1
    sensitivity = .99
    lens = [expected_overlap_len(n, n, diag, g) for diag in range(-n, n)]
    radii = band_radii(lens, g, sensitivity)

    assert len(radii) == 2 * n, 'correct number of diagonals'
    assert np.all(np.diff(radii[0: n]) >= 0), \
        'overlap band radius should increase as d goes from -|T| to 0'
    assert np.all(np.diff(radii[n:]) <= 0), \
        'overlap band radius should increase as d goes from 0 to |S|'


@pytest.mark.parametrize('wordlen', [8, 15],
                         ids=['k=8', 'k=15'])
@pytest.mark.parametrize('K', [500, 1000],
                         ids=['K=500', 'K=1000'])
@pytest.mark.parametrize('n', [2000, 5000],
                         ids=['n=1000', 'n=5000'])
def test_local_similarity(wordlen, K, n):
    gap, subst = .05, .05
    A = Alphabet('ACGT')
    M = MutationProcess(A, subst_probs=subst, ge_prob=gap, go_prob=gap)
    HF_kw = {'g_max': .2, 'sensitivity': .99, 'alphabet': A,
             'wordlen': wordlen, 'path': ':memory:'}

    hom = rand_seq(A, K)
    S = A.parse(str(hom) + str(rand_seq(A, n - K)))
    T = A.parse(str(M.mutate(hom)[0]) + str(rand_seq(A, n-K)))
    HF = HomologyFinder(S, T, **HF_kw)

    p_match = (1 - gap) * (1 - subst) * .9
    diagonal_scores = HF.score_diagonal_bands(K, p_match)
    assert abs(np.nanargmax(diagonal_scores[:, 1]) - HF.d0) < K, \
        'The highest H1 scoring diagonal band must be roughly around 0'

    found_homs = list(HF.similar_segments(K, p_match))
    assert len(found_homs) == 1, 'Only one similar segment should be found'
    ((d_min, d_max), (a_min, a_max)), score, match_p = found_homs[0]
    assert d_min < 10 and d_max > -10 and a_min < K, \
        'The coordinates of the homologous segment must be correct'
    assert 0 <= match_p <= 1, 'estimated match probability should be in [0, 1]'


@pytest.mark.parametrize('wordlen', [8, 15],
                         ids=['k=8', 'k=15'])
@pytest.mark.parametrize('K', [500, 1000],
                         ids=['K=500', 'K=1000'])
@pytest.mark.parametrize('n', [2000, 5000],
                         ids=['n=1000', 'n=5000'])
def test_overlap_detection(wordlen, K, n):
    gap, subst = .05, .05
    A = Alphabet('ACGT')
    M = MutationProcess(A, subst_probs=subst, ge_prob=gap, go_prob=gap)
    HF_kw = {'g_max': .2, 'sensitivity': .99, 'alphabet': A,
             'wordlen': wordlen, 'path': ':memory:'}

    p_match = (1 - gap) * (1 - subst)
    overlap = rand_seq(A, K)
    S = rand_seq(A, n - K) + overlap
    T = M.mutate(overlap)[0] + rand_seq(A, n - K)
    HF = HomologyFinder(S, T, **HF_kw)
    (seg0, s0), (seg1, s1) = HF.highest_scoring_overlap_band(p_match)
    assert seg0[0] < n - K < seg0[1], 'H0 must detect the correct overlap'
    assert seg1[0] < n - K < seg1[1], 'H1 must detect the correct overlap'
