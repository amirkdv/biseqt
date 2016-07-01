# -*- coding: utf-8 -*-
import pytest
import mock

from biseqt.sequence import Alphabet
from biseqt.random import rand_seq, MutationProcess, rand_read
from biseqt.random import np  # to mock


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
