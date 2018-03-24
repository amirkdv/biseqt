# -*- coding: utf-8 -*-
import pytest
from itertools import product
from random import choice

from biseqt.stochastics import rand_seq
from biseqt.sequence import Alphabet, Sequence
from biseqt.kmers import kmer_as_int, as_kmer_seq, KmerIndex, KmerCache


def test_kmer_as_int_limitations():
    with pytest.raises(AssertionError):
        # word length too large
        A = Alphabet('ACGT')
        KmerIndex(':memory:', wordlen=32, alphabet=A)

    with pytest.raises(AssertionError):
        # alphabet too large
        A = Alphabet(str(i) for i in range(37))
        KmerIndex(':memory:', wordlen=5, alphabet=A)


@pytest.mark.parametrize('alphabet',
                         [Alphabet('ACGT'), Alphabet(['00', '01', '11'])],
                         ids=['one-letter alphabet', 'two-letter alphabet'])
@pytest.mark.parametrize('wordlen', [3, 6, 9],
                         ids=['w = %d' % i for i in [3, 6, 9]])
def test_kmer_as_int(alphabet, wordlen):
    kmers = product(range(len(alphabet)), repeat=wordlen)
    as_ints = [kmer_as_int(kmer, alphabet) for kmer in kmers]
    assert len(set(as_ints)) == len(alphabet) ** wordlen, \
        'for a fixed k the integer representation of kmers must be unique'


@pytest.mark.parametrize('alphabet',
                         [Alphabet('ACGT'), Alphabet(['00', '01', '11'])],
                         ids=['one-letter alphabet', 'two-letter alphabet'])
@pytest.mark.parametrize('wordlen', [3, 6, 9],
                         ids=['w = %d' % i for i in [3, 6, 9]])
def test_kmer_as_int_masked(alphabet, wordlen):
    S = Sequence(alphabet, contents=tuple([0] * 10))
    mask = [set([0])]
    assert all(i is None for i in as_kmer_seq(S, wordlen, mask=mask)), \
        'masking %s-only kmers should leave no kmers in %s' % (alphabet[0], S)
    S = Sequence(alphabet, contents=tuple([0] * 10 + [1]))
    n = sum(i for i in as_kmer_seq(S, wordlen, mask=mask) if i is not None)
    assert n == 1, 'masking %s-only kmers should leave 1 kmer in %s' % S

    mask = [set([1]), set([2]), set([1, 2])]  # e.g. CG mask
    S = Sequence(alphabet,
                 contents=tuple([choice([1, 2]) for _ in range(10)] + [0]))
    n = sum(int(i is not None) for i in as_kmer_seq(S, wordlen, mask=mask))
    assert n == 1, 'masking %s-only kmers should leave one %d-mer in %s' % \
                   ((alphabet[1] + alphabet[2]), wordlen, S)


@pytest.mark.parametrize('alphabet',
                         [Alphabet('ACGT'), Alphabet(['00', '01', '11'])],
                         ids=['one-letter alphabet', 'two-letter alphabet'])
@pytest.mark.parametrize('wordlen', [3, 13, 23],
                         ids=['k=%d' % i for i in [3, 6, 9]])
def test_as_kmer_seq(alphabet, wordlen):
    S = rand_seq(alphabet, 50)
    kmer_seq = as_kmer_seq(S, wordlen)
    assert all(isinstance(i, int) for i in kmer_seq), \
        'kmer representations should be integers'
    assert all(i < len(alphabet) ** wordlen for i in kmer_seq), \
        'kmer representations should be valid integers in base |alphabet|'
    assert len(kmer_seq) == len(S) - wordlen + 1, \
        'correct number of kmers should be scanned'


@pytest.fixture(ids=['wordlen 3', 'wordlen 13'], params=[3, 13])
def dna_kmer_index(request):
    """Returns a kmer index created on top of a sequence database (i.e
    ``biseqt.database.DB``) with DNA alphabet, parameterized word length,
    stored in memory."""
    A = Alphabet('ACGT')
    return KmerIndex(path=':memory:', alphabet=A, wordlen=request.param)


def test_index_kmers(dna_kmer_index):
    S = rand_seq(dna_kmer_index.alphabet, 50)
    seqid = dna_kmer_index.index_kmers(S)

    kmers = dna_kmer_index.kmers()
    assert set(kmers) == set(as_kmer_seq(S, dna_kmer_index.wordlen)), \
        'correct set of unique kmers should be observed'
    total_hits = sum(len(dna_kmer_index.hits(kmer)) for kmer in kmers)
    assert total_hits == len(S) - dna_kmer_index.wordlen + 1, \
        'correct number of hits should be oberved'

    assert seqid == dna_kmer_index.index_kmers(S), \
        'same sequence should not be indexed twice'
    T = rand_seq(dna_kmer_index.alphabet, 10)
    assert seqid != dna_kmer_index.index_kmers(T), \
        'different sequences should have different seqids'


def test_kmer_cache():
    A = Alphabet('ACGT')
    S = rand_seq(A, 50)
    wordlen = 5
    cache = KmerCache(path=':memory:', wordlen=wordlen, alphabet=A)
    assert cache.as_kmer_seq(S) == as_kmer_seq(S, wordlen), \
        'kmer cache must produce same results as as_kmer_seq()'
    assert len(cache.cached_seqs()) == 1, \
        'cache must be populated'
