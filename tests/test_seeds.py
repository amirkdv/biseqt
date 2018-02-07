# -*- coding: utf-8 -*-
from biseqt.sequence import Alphabet
from biseqt.seeds import SeedIndex


def test_coordinate_change():
    pairs = {  # (i,j) => (d, a)
        (0, 0):     (0, 0),
        (10, 0):    (10, 0),
        (0, 10):    (-10, 0),
        (10, 10):   (0, 10)
    }
    for (i, j), (d, a) in pairs.items():
        assert SeedIndex.to_diagonal_coordinates(i, j) == (d, a), \
            'conversion to diagonal coordinates should be correct'
        assert SeedIndex.to_ij_coordinates(d, a) == (i, j), \
            'conversion from diagonal coordinates should be correct'


def test_index_seeds():
    A = Alphabet('ACGT')
    wordlen = 5
    kw = {'alphabet': A, 'wordlen': wordlen, 'path': ':memory:'}

    S = A.parse('G' * wordlen)
    T = A.parse('TC' + str(S))
    assert list(SeedIndex(S, T, **kw).seeds()) == [(0, 2)] and \
        list(SeedIndex(T, S, **kw).seeds()) == [(2, 0)], \
        'seed coordinates should be correctly calculated'

    S = A.parse('A' * 10)
    T = A.parse('A' * 20)
    n_seeds = (len(S) - wordlen + 1) * (len(T) - wordlen + 1)
    assert len(list(SeedIndex(S, T, **kw).seeds())) == n_seeds, \
        'correct number of seeds should be found'
    n_kmers = (len(S) - wordlen + 1)
    n_seeds = n_kmers * n_kmers - n_kmers  # don't count (i, i) twice
    assert len(list(SeedIndex(S, S, **kw).seeds())) == n_seeds, \
        'self comparison should not get confused'


def test_count_seeds_on_diagonals():
    A = Alphabet('ACGT')
    wordlen = 5
    kw = {'alphabet': A, 'wordlen': wordlen, 'path': ':memory:'}
    S = A.parse('G' * wordlen)
    T = A.parse('TT' + str(S))
    d0 = len(T) - 1

    count_by_d_ = SeedIndex(S, T, **kw).seed_count_by_d_()
    expected_count = [1 if i - d0 == -2 else 0
                      for i in range(len(S) + len(T))]
    assert list(count_by_d_) == expected_count, \
        'diagonal counts should be correctly calculated'
