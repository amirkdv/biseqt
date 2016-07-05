# -*- coding: utf-8 -*-
import pytest
from tempfile import NamedTemporaryFile
from itertools import product
from scipy.stats import norm, binom
from math import log
import numpy as np

from biseqt.random import rand_seq
from biseqt.sequence import Alphabet
from biseqt.database import DB
from biseqt.kmers import binomial_to_normal, normal_neg_log_pvalue
from biseqt.kmers import KmerIndex


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


def test_kmer_as_int_limitations():
    A = Alphabet('ACGT')
    with pytest.raises(AssertionError):
        # word length too large
        KmerIndex(DB(A, 'example.db'), 32)

    with pytest.raises(AssertionError):
        # alphabet too large
        KmerIndex(DB(Alphabet(str(i) for i in range(37))), 5)


@pytest.fixture(ids=['one-letter alphabet', 'two-letter alphabet'],
                params=[Alphabet('ACGT'), Alphabet(['00', '01', '11'])])
def db_gen(request):
    """Returns a function that generates a sequence database (i.e
    ``biseqt.database.DB``) with parametrized alphabet (single letter and
    double letter) stored in a temporary file."""
    A = request.param

    def f():
        with NamedTemporaryFile() as tmp:
            db = DB(tmp.name, A)
            yield db
    return f


@pytest.mark.parametrize('wordlen', [3, 6, 9],
                         ids=['k = %d' % i for i in [3, 6, 9]])
def test_kmer_as_int(db_gen, wordlen):
    db = next(db_gen())
    A = db.alphabet
    kmer_index = KmerIndex(db, wordlen)
    kmers = product(range(len(A)), repeat=wordlen)
    as_ints = [kmer_index.kmer_as_int(kmer) for kmer in kmers]
    assert len(set(as_ints)) == len(A) ** wordlen, \
        'for a fixed k the integer representation of kmers must be unique'


@pytest.mark.parametrize('wordlen', [3, 6, 9],
                         ids=['k = %d' % i for i in [3, 6, 9]])
def test_scan_kmers(db_gen, wordlen):
    db = next(db_gen())
    A = db.alphabet
    kmer_index = KmerIndex(db, wordlen)
    S = rand_seq(A, 50)
    assert sum(1 for _ in kmer_index.scan_kmers(S)) == len(S) - wordlen + 1, \
        'correct number of kmers should be scanned'


def test_index_kmers():
    A = Alphabet('ACGT')
    wordlen = 3
    with NamedTemporaryFile() as tmp:
        db = DB(tmp.name, A)
        kmer_index = KmerIndex(db, wordlen)
        db.initialize()

        S = A.parse('ATGCA', name='foo')
        S_id = db.insert(S).id
        assert kmer_index.num_kmers() == 3

        T = A.parse('ATGCC', name='bar')
        T_id = db.insert(T).id
        assert kmer_index.num_kmers() == 4
        assert kmer_index.total_length_indexed() == len(S) + len(T)

        # find the occurences of 'ATG'
        atg_int = kmer_index.kmer_as_int((0, 3, 2))
        atg_hits = [(hits, score) for kmer, hits, score in kmer_index.kmers()
                    if kmer == atg_int]
        assert len(atg_hits) == 1 and atg_hits[0][0] == [(S_id, 0), (T_id, 0)]

        # calculate scores
        assert atg_hits[0][1] is None
        kmer_index.score_kmers()
        atg_hits = [(hits, score) for kmer, hits, score in kmer_index.kmers()
                    if kmer == atg_int]
        atg_score = atg_hits[0][1]
        assert atg_score is not None and score > 0

        # ATG is the most common kmer, it most have the highest score:
        assert all(atg_score >= score for _, _, score in kmer_index.kmers())

        # shouldn't do anything
        kmer_index.score_kmers(only_missing=True)
