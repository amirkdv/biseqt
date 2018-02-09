# -*- coding: utf-8 -*-
import pytest

from tempfile import NamedTemporaryFile
from biseqt.stochastics import rand_seq
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


@pytest.mark.parametrize('in_memory', [True, False],
                         ids=['in memory', 'on disk'])
@pytest.mark.parametrize('wordlen', [3, 13],
                         ids=['k=3', 'k=13'])
def test_index_integrity(in_memory, wordlen):
    def _tests(path):
        A = Alphabet('ACGT')
        kw = {'alphabet': A, 'wordlen': wordlen, 'path': path}

        S = rand_seq(A, 2 * wordlen)
        T = rand_seq(A, 2 * wordlen)
        index1 = SeedIndex(S, T, **kw)
        index2 = SeedIndex(S, T, **kw)
        assert set(index1.seeds()) == set(index2.seeds()), \
            'Multiple index objects for the same pair should not conflict'

        # S = {ACG}*, T = TTTT [S] TTTT
        # NOTE S must not have non-trivial seeds with itself
        S = A.parse('AAACCCGGGCAAGCC')
        T = A.parse(('T' * 2 * wordlen) + str(S) + ('T' * 2 * wordlen))
        index3 = SeedIndex(S, T, **kw)
        assert len(list(index3.seeds())) == len(S) - wordlen + 1,\
            'Multiple comparisons should work on the same datbase path'

    if in_memory:
        _tests(':memory:')
    else:
        with NamedTemporaryFile() as f:
            _tests(f.name)


@pytest.mark.parametrize('in_memory', [True, False],
                         ids=['in memory', 'on disk'])
@pytest.mark.parametrize('wordlen', [3, 13],
                         ids=['k=3', 'k=13'])
def test_index_seeds(in_memory, wordlen):
    def _tests(path):
        A = Alphabet('ACGT')
        kw = {'alphabet': A, 'wordlen': wordlen, 'path': path}

        S = A.parse('G' * wordlen)
        T = A.parse('TC' + str(S))
        assert list(SeedIndex(S, T, **kw).seeds()) == [(0, 2)] and \
            list(SeedIndex(T, S, **kw).seeds()) == [(2, 0)], \
            'seed coordinates should be correctly calculated'

        S = A.parse('A' * 5 * wordlen)
        T = A.parse('A' * 10 * wordlen)
        n_seeds = (len(S) - wordlen + 1) * (len(T) - wordlen + 1)
        assert len(list(SeedIndex(S, T, **kw).seeds())) == n_seeds, \
            'correct number of seeds should be found'
        n_kmers = (len(S) - wordlen + 1)
        n_seeds = n_kmers * n_kmers
        assert len(list(SeedIndex(S, S, **kw).seeds())) == n_seeds, \
            'self comparison should not get confused'
        assert len(list(SeedIndex(S, S, **kw).seeds())) == n_seeds, \
            'repeated self comparison should not get confused'

    if in_memory:
        _tests(':memory:')
    else:
        with NamedTemporaryFile() as f:
            _tests(f.name)


@pytest.mark.parametrize('in_memory', [True, False],
                         ids=['in memory', 'on disk'])
@pytest.mark.parametrize('wordlen', [3, 13],
                         ids=['k=3', 'k=13'])
def test_seed_counts(in_memory, wordlen):
    def _tests(path):
        A = Alphabet('ACGT')
        kw = {'alphabet': A, 'wordlen': wordlen, 'path': path}

        S = rand_seq(A, 10 * wordlen)
        T = rand_seq(A, 10 * wordlen)
        seed_index = SeedIndex(S, T, **kw)
        assert len(list(seed_index.seeds())) == seed_index.seed_count(), \
            'seeds() and seed_count() should agree'
        for d in range(-wordlen, wordlen):
            n1 = len(list(seed_index.seeds(d_center=d, d_radius=wordlen)))
            n2 = seed_index.seed_count(d_center=d, d_radius=wordlen)
            assert n1 == n2, 'seeds() and seed_count() should agree on bands'

        S = A.parse('G' * wordlen)
        T = A.parse('TT' + str(S))
        d0 = len(T) - 1

        seed_index = SeedIndex(S, T, **kw)
        count_by_d_ = seed_index.seed_count_by_d_()
        expected_count = [1 if i - d0 == -2 else 0
                          for i in range(len(S) + len(T) - 1)]
        assert list(count_by_d_) == expected_count, \
            'diagonal counts should be correctly calculated'

        S = A.parse('T' * wordlen + 'G' * wordlen)
        T = A.parse('G' * wordlen + 'T' * wordlen)
        seed_index = SeedIndex(S, T, **kw)
        assert seed_index.seed_count() == 2 and \
            seed_index.seed_count(d_center=-wordlen, d_radius=1) == 1 and \
            seed_index.seed_count(d_center=wordlen, d_radius=1) == 1, \
            'diagonal counts should be correctly calculated'

    if in_memory:
        _tests(':memory:')
    else:
        with NamedTemporaryFile() as f:
            _tests(f.name)
