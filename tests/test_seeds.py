# -*- coding: utf-8 -*-
import pytest

from tempfile import NamedTemporaryFile
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
def test_index_integrity(in_memory):
    def _tests(path):
        A = Alphabet('ACGT')
        wordlen = 5
        kw = {'alphabet': A, 'wordlen': wordlen, 'path': path}

        S = A.parse('ACGGGCTTTTCG')
        T = A.parse('GTTTCTGGGAGC')
        index1 = SeedIndex(S, T, **kw)
        index2 = SeedIndex(S, T, **kw)
        assert set(index1.seeds()) == set(index2.seeds()), \
            'Multiple index objects for the same pair should not conflict'

        T = A.parse('GGG' + str(S))
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
def test_index_seeds(in_memory):
    def _tests(path):
        A = Alphabet('ACGT')
        wordlen = 5
        kw = {'alphabet': A, 'wordlen': wordlen, 'path': path}

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
def test_count_seeds_on_diagonals(in_memory):
    def _tests(path):
        A = Alphabet('ACGT')
        wordlen = 5
        kw = {'alphabet': A, 'wordlen': wordlen, 'path': path}

        S = A.parse('G' * wordlen)
        T = A.parse('TT' + str(S))
        d0 = len(T) - 1

        count_by_d_ = SeedIndex(S, T, **kw).seed_count_by_d_()
        expected_count = [1 if i - d0 == -2 else 0
                          for i in range(len(S) + len(T))]
        assert list(count_by_d_) == expected_count, \
            'diagonal counts should be correctly calculated'

    if in_memory:
        _tests(':memory:')
    else:
        with NamedTemporaryFile() as f:
            _tests(f.name)
