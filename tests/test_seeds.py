# -*- coding: utf-8 -*-
import pytest

from biseqt.sequence import Alphabet
from biseqt.seeds import SeedIndex


@pytest.fixture
def seed_index():
    """Creates a database, a kmer index, and a seed index with word length 5
    stored in memory and returns the seed index. The database is populated with
    3 random sequences of length 100 and all kmers and seeds are indexed."""
    A = Alphabet('ACGT')
    wordlen = 5

    seed_index = SeedIndex(path=':memory:', wordlen=wordlen, alphabet=A)
    return seed_index


def test_index_seeds(seed_index):
    A = seed_index.alphabet
    wordlen = seed_index.wordlen
    S = A.parse('G' * wordlen)
    T = A.parse('TC' + str(S))
    assert list(seed_index.seeds(S, T)) == [(0, 2)] and \
        list(seed_index.seeds(T, S)) == [(2, 0)], \
        'seed coordinates should be correctly calculated'

    S = seed_index.alphabet.parse('A' * 10)
    T = seed_index.alphabet.parse('A' * 20)
    n_seeds = (len(S) - wordlen + 1) * (len(T) - wordlen + 1)
    assert len(list(seed_index.seeds(S, T))) == n_seeds, \
        'correct number of seeds should be found'
    n_kmers = (len(S) - wordlen + 1)
    n_seeds = n_kmers * n_kmers - n_kmers  # don't count (i, i) twice
    assert len(list(seed_index.seeds(S, S))) == n_seeds, \
        'self comparison should not get confused'


def test_count_seeds_on_diagonals(seed_index):
    A = seed_index.alphabet
    wordlen = seed_index.wordlen
    S = A.parse('G' * wordlen)
    T = A.parse('TT' + str(S))
    d0 = len(T) - 1

    count_by_d_ = seed_index.seed_count_by_d_(S, T)
    expected_count = [1 if i - d0 == -2 else 0
                      for i in range(len(S) + len(T))]
    assert list(count_by_d_) == expected_count, \
        'diagonal counts should be correctly calculated'
