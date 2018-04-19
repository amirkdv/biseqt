# -*- coding: utf-8 -*-
import pytest

from tempfile import NamedTemporaryFile
from biseqt.stochastics import rand_seq
from biseqt.sequence import Alphabet
from biseqt.seeds import SeedIndex, SeedIndexMultiple


def test_coordinate_change():
    pairs = {  # (i,j) => (d, a)
        (0, 0): (0, 0),
        (0, 1): (-1, 1),
        (1, 0): (1, 1),
        (1, 1): (0, 2)
    }
    for (i, j), (d, a) in pairs.items():
        assert SeedIndex.to_diagonal_coordinates(i, j) == (d, a), \
            '(i, j) = (%d, %d) is the same as (d, a) = (%d, %d)' % \
            (i, j, d, a)
        assert SeedIndex.to_ij_coordinates(d, a) == (i, j), \
            '(d, a) = (%d, %d) is the same as (i, j) = (%d, %d)' % \
            (d, a, i, j)

    # (d_min, d_max), (a_min, a_max)
    seg = (-2, 2), (0, 2)
    (i_start, i_end), (j_start, j_end) = SeedIndex.to_ij_coordinates_seg(seg)
    assert i_start == j_start == 0, \
        'segment %s in diagonal coordinates contains (0, 0)' % str(seg)
    assert i_end == j_end == 2, \
        'segment %s in diagonal coordinates contains (2, 2)' % str(seg)


def test_coordinate_change_multiple():
    pairs = {  # (i,j, k) => (d1, d2, a)
        (0, 0, 0): ((0, 0), 0),
        (0, 0, 1): ((0, -1), 1),
        (0, 1, 0): ((-1, 0), 1),
        (0, 1, 1): ((-1, -1), 2),
        (1, 0, 0): ((1, 1), 1),
        (1, 0, 1): ((1, 0), 2),
        (1, 1, 0): ((0, 1), 2),
        (1, 1, 1): ((0, 0), 3),
    }
    for std, diag in pairs.items():
        print std, diag
        assert SeedIndexMultiple.to_diagonal_coordinates(*std) == diag, \
            '(%s) standard == %s diagonal' % (str(std), str(diag))
        assert SeedIndexMultiple.to_ij_coordinates(*diag) == std, \
            '(%s) diag == %s standard' % (str(diag), str(std))

    # (d_min, d_max), (a_min, a_max)
    seg = ((-2, 2), (-2, 2)), (0, 2)
    std_ranges = SeedIndexMultiple.to_ij_coordinates_seg(seg)
    for std_range in std_ranges:
        assert std_range == (0, 2), \
            'segment %s in diagonal coordinates projects to (0, 2)' % str(seg)


@pytest.mark.parametrize('in_memory', [True, False],
                         ids=['in memory', 'on disk'])
@pytest.mark.parametrize('wordlen', [5, 15], ids=['k=5', 'k=15'])
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
@pytest.mark.parametrize('wordlen', [5, 15], ids=['k=5', 'k=15'])
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


@pytest.mark.parametrize('wordlen', [5, 15], ids=['k=5', 'k=15'])
def test_seed_counts(wordlen):
    A = Alphabet('ACGT')
    kw = {'alphabet': A, 'wordlen': wordlen, 'path': ':memory:'}

    S = rand_seq(A, 5 * wordlen)
    T = rand_seq(A, 5 * wordlen)
    seed_index = SeedIndex(S, T, **kw)
    assert len(list(seed_index.seeds())) == seed_index.seed_count(), \
        'seeds() and seed_count() should agree'
    for d in range(-wordlen, wordlen):
        n1 = len(list(seed_index.seeds(d_band=(d - wordlen, d + wordlen))))
        n2 = seed_index.seed_count(d_band=(d - wordlen, d + wordlen))
        assert n1 == n2, \
            'seeds() and seed_count() should agree on diagonal bands'

    T = A.parse(str(S))
    seed_index = SeedIndex(S, T, **kw)
    band = (0, 0)
    n1 = seed_index.seed_count(d_band=band)
    assert n1 == len(S) - wordlen + 1, \
        'self comparison should find all seeds on diagonal band'
    n2 = len(list(seed_index.seeds(d_band=band, exclude_trivial=True)))
    assert n2 == 0, \
        'excluding trivial seeds from self comparison should work'

    S = A.parse('T' * wordlen + 'G' * wordlen)
    T = A.parse('G' * wordlen + 'T' * wordlen)
    seed_index = SeedIndex(S, T, **kw)
    assert seed_index.seed_count() == 2, \
        '%s and %s have 2 seeds of length %d' % (S, T, wordlen)
    d_band = (-wordlen - 1, -wordlen + 1)
    assert seed_index.seed_count(d_band=d_band) == 1, \
        '%s and %s have 1 seed of length %d in diagonal band %s' % \
        (S, T, wordlen, str(d_band))
    d_band = (wordlen - 1, wordlen + 1)
    assert seed_index.seed_count(d_band=d_band) == 1, \
        '%s and %s have 1 seed of length %d in diagonal band %s' % \
        (S, T, wordlen, str(d_band))
    a_band = (wordlen, wordlen)
    assert seed_index.seed_count(a_band=a_band) == 2, \
        '%s and %s have 1 seed of length %d in antidiagonal band %s' % \
        (S, T, wordlen, str(a_band))

    T = A.parse('C' * wordlen + 'A' * wordlen)
    seed_index = SeedIndex(S, T, **kw)
    assert seed_index.seed_count() == 0, \
        '%s and %s have no seeds' % (S, T)


@pytest.mark.parametrize('n_seqs', [5, 15], ids=['n=5', 'n=15'])
@pytest.mark.parametrize('wordlen', [5, 15], ids=['k=5', 'k=15'])
def test_seed_counts_multiple(n_seqs, wordlen):
    A = Alphabet('ACGT')
    kw = {'alphabet': A, 'wordlen': wordlen, 'path': ':memory:'}

    S = rand_seq(A, 5 * wordlen)
    seqs = [S + rand_seq(A, 5 * wordlen) for _ in range(n_seqs)]
    seed_index = SeedIndexMultiple(*seqs, **kw)
    a_band = (0, n_seqs * (len(S) - wordlen))
    ds_band = [(0, 0) for _ in range(n_seqs - 1)]
    n_seeds = seed_index.seed_count(a_band=a_band, ds_band=ds_band)
    assert n_seeds == len(S) - wordlen + 1, \
        'number of seeds for multiple sequences should be correct'
