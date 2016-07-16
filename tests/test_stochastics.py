# -*- coding: utf-8 -*-
import pytest
import mock
from itertools import combinations
from scipy.stats import norm, binom
from math import log

from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess, rand_read
from biseqt.stochastics import binomial_to_normal, normal_neg_log_pvalue
from biseqt.stochastics import band_radius
from biseqt.stochastics import np  # to mock


def test_rand_seq():
    _bak = np.random.choice
    np.random.choice = mock.Mock(return_value=[0, 0, 0])
    A = Alphabet('ACGT')
    assert rand_seq(A, 10) == A.parse('AAA')
    np.random.choice = _bak


def test_lossless_reads():
    A = Alphabet('ACGT')
    S = rand_seq(A, 100)
    with pytest.raises(AssertionError):
        next(rand_read(S, len_mean=200, num=1))  # len_mean must be < len(S)
    with pytest.raises(AssertionError):
        # at most one of num or expected_coverage given
        next(rand_read(S, len_mean=50, num=1, expected_coverage=1))

    assert sum(1 for _ in rand_read(S, len_mean=50, num=10)) == 10, \
        'The number of sampled reads should be controllable'
    assert sum(1 for _ in rand_read(S, len_mean=50)) == 1, \
        'If neither num or expected coverage is given only one sample is read'

    # there should be no noise added
    read, pos = next(rand_read(S, len_mean=40, num=1))
    assert S[pos:pos+len(read)] == read

    S = A.parse('ACT' * 100)
    reads = [x for x in rand_read(S, len_mean=100, len_sd=0.01, num=100)]
    assert set(len(read) for read, _ in reads) > 1, \
        'Read lengths should be randomly chosen'
    len_mean = sum(len(read) for read, _ in reads) / 100.
    assert len_mean > 50 and len_mean < 150, \
        'Normal distribution of read lengths works'

    # index edge cases
    A = Alphabet(['00', '01'])
    S = A.parse('01' * 10)
    _bak = np.random.normal
    np.random.normal = mock.Mock(return_value=[1])
    assert next(rand_read(S, len_mean=1, num=1))[0] == A.parse('01'), \
        'sequences in alphabets with > 1 long letters can be sampled too'
    np.random.normal = _bak


def test_lossy_reads():
    A = Alphabet('ACGT')
    S = A.parse('ACT' * 100)
    gap_kw = {'go_prob': 0.2, 'ge_prob': 0.3}
    M = MutationProcess(A, subst_probs=0.3, **gap_kw)
    read, pos, tx = next(M.noisy_read(S, len_mean=50, num=1))
    assert tx.count('S') > 0 and tx.count('I') + tx.count('D') > 0, \
        'Random mutations should be performed to get lossy reads'


def test_binomial_to_normal():
    n, p = 10000, 0.1
    mu, sd = binomial_to_normal(n, p)
    binomial = [binom.cdf(x, n, p) for x in range(0, 10000, 500)]
    normal = [norm(mu, sd).cdf(x) for x in range(0, 10000, 500)]
    assert np.allclose(binomial, normal, atol=0.01)


def test_normal_neg_log_pvalue():
    sds = range(1, 10)
    neg_log_pvalue = [normal_neg_log_pvalue(0, sd, 0) for sd in sds]
    assert np.allclose(neg_log_pvalue, [-log(0.5)] * len(sds))


def test_band_radius():
    len0 = len1 = 10
    diag = 0
    assert band_radius(len0, len1, diag, gap_prob=1e-10, sensitivity=0.9) \
        == 1, 'band radius can be shrunk down to 1'

    args = (len0, len1, diag)

    for sensitivity in [i/10. for i in range(1, 10)]:
        gap_probs = [1e-10, 1e-8, 1e-4, 1e-2, 1e-1]
        radii = [band_radius(*args, gap_prob=g, sensitivity=sensitivity)
                 for g in gap_probs]
        assert all(x <= y for x, y in zip(radii, radii[1:])), \
            'band radius must be increasing with increasing gap probability'

    for gap_prob in [i/100. for i in range(1, 10)]:
        sensitivities = [.5 + .05 * i for i in range(10)]
        radii = [band_radius(*args, gap_prob=gap_prob, sensitivity=s)
                 for s in sensitivities]
        assert all(x <= y for x, y in zip(radii, radii[1:])), \
            'band radius must be increasing with increasing sensitivity'

    kw = {'gap_prob': 0.2, 'sensitivity': 0.9}
    radii = [band_radius(len0, len1, d, **kw) for d in range(-len1, 0)]
    assert all(x <= y for x, y in zip(radii, radii[1:])), \
        'band radius must be increasing with increasing diag if diag < 0'

    kw = {'gap_prob': 0.2, 'sensitivity': 0.9}
    radii = [band_radius(len0, len1, d, **kw) for d in range(0, len0)]
    assert all(x >= y for x, y in zip(radii, radii[1:])), \
        'band radius must be decreasing with increasing diag if diag > 0'


def test_mutation_process():
    A = Alphabet('ACGT')
    S = A.parse('ACT' * 100)
    gap_kw = {'go_prob': 0, 'ge_prob': 0}
    T, tx = MutationProcess(A, subst_probs=0, **gap_kw).mutate(S)
    assert T == S and tx == 'MMM' * 100, \
        'all mutation probabilities can be set to zero'

    T, tx = MutationProcess(A, subst_probs=0.1, **gap_kw).mutate(S)
    assert all(op in 'MS' for op in tx) and 'S' in tx, \
        'there can be mutation processes with only substitutions'

    T, tx = MutationProcess(A, subst_probs=0.01, **gap_kw).mutate(S)
    assert tx.count('S') < 0.1 * len(S), 'substitution probabilities work'

    with pytest.raises(AssertionError):
        MutationProcess(A, go_prob=0.2, ge_prob=0.1)  # go_prob <= ge_prob

    gap_kw = {'go_prob': 0.05, 'ge_prob': 0.1}
    T, tx = MutationProcess(A, subst_probs=0, **gap_kw).mutate(S)
    indels = sum(1 for op in tx if op in 'ID')
    assert indels > 0 and indels < 0.5 * len(S), 'gap probabilities work'


def test_log_odds_scores():
    A = Alphabet('ACGT')
    # linear gap model
    P = MutationProcess(A, subst_probs=.1, ge_prob=.1, go_prob=.1)
    subst_scores, (go_score, ge_score) = P.log_odds_scores()
    assert go_score == 0. and ge_score < 0
    match_pos = [(i, i) for i in range(len(A))]
    mismatch_pos = [(i, j) for i, j in combinations(range(len(A)), 2)]
    assert all(subst_scores[i][j] < 0 for i, j in mismatch_pos)
    assert all(subst_scores[i][j] > 0 for i, j in match_pos)

    # affine gap model
    P = MutationProcess(A, subst_probs=.1, ge_prob=.2, go_prob=.1)
    subst_scores, (go_score, ge_score) = P.log_odds_scores()
    assert ge_score < 0

    # do mismatch scores go down if subst probs are decreased?
    P = MutationProcess(A, subst_probs=.01, ge_prob=.2, go_prob=.1)
    new_subst_scores, _ = P.log_odds_scores()
    assert new_subst_scores[0][1] < subst_scores[0][1], \
        'mismatch scores become more negative with lower mismatch probs'
