# -*- coding: utf-8 -*-
import pytest
from itertools import product
from StringIO import StringIO

from biseqt.random import rand_seq
from biseqt.sequence import Alphabet
from biseqt.database import DB
from biseqt.kmers import KmerIndex
from biseqt.io import write_fasta


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
def seq_db(request):
    """Returns a a sequence database (i.e ``biseqt.database.DB``) with
    parametrized alphabet (single letter and double letter) stored in
    memory."""
    return DB(':memory:', request.param)


@pytest.mark.parametrize('wordlen', [3, 6, 9],
                         ids=['k = %d' % i for i in [3, 6, 9]])
def test_kmer_as_int(seq_db, wordlen):
    db = seq_db
    A = db.alphabet
    kmer_index = KmerIndex(db, wordlen)
    kmers = product(range(len(A)), repeat=wordlen)
    as_ints = [kmer_index.kmer_as_int(kmer) for kmer in kmers]
    assert len(set(as_ints)) == len(A) ** wordlen, \
        'for a fixed k the integer representation of kmers must be unique'


@pytest.mark.parametrize('wordlen', [3, 6, 9],
                         ids=['k = %d' % i for i in [3, 6, 9]])
def test_scan_kmers(seq_db, wordlen):
    db = seq_db
    A = db.alphabet
    kmer_index = KmerIndex(db, wordlen)
    S = rand_seq(A, 50)
    assert sum(1 for _ in kmer_index.scan_kmers(S)) == len(S) - wordlen + 1, \
        'correct number of kmers should be scanned'


@pytest.fixture(ids=['wordlen 3'], params=[3])
def dna_kmer_index(request):
    """Returns a kmer index created on top of a sequence database (i.e
    ``biseqt.database.DB``) with DNA alphabet, parameterized word length,
    stored in memory."""
    return KmerIndex(DB(':memory:', Alphabet('ACGT')), request.param)


def test_index_kmers(dna_kmer_index):
    db = dna_kmer_index.db
    A = db.alphabet

    db.initialize()
    S = A.parse('ATGCAGGGC', name='foo')  # has no repeating kmers
    rec = db.insert(S)
    assert next(dna_kmer_index.scanned_sequences()) == (rec.id, len(S))
    assert dna_kmer_index.num_kmers() == len(S) - dna_kmer_index.wordlen + 1

    num = dna_kmer_index.total_length_indexed()
    dna_kmer_index.index_kmers(db.connection(), S, rec)
    assert dna_kmer_index.total_length_indexed() == num, \
        'calling inder_kmers on already-indexed sequences has no consequences'


def test_kmer_hits(dna_kmer_index):
    db, wordlen = dna_kmer_index.db, dna_kmer_index.wordlen
    A = db.alphabet

    db.initialize()
    S = A.parse('A' * 10, name='foo')
    db.insert(S)

    kmer_contents = S[:wordlen].contents
    num_hist = sum(1 for _ in dna_kmer_index.hits(kmer_contents))
    assert num_hist == len(S) - wordlen + 1, \
        'hits() should work with tuple representation of kmer'

    kmer_as_int = dna_kmer_index.kmer_as_int(S[:wordlen])
    num_hist = sum(1 for _ in dna_kmer_index.hits(kmer_as_int))
    assert num_hist == len(S) - wordlen + 1, \
        'hits() should work with integer representation of kmer'


def test_score_kmers(dna_kmer_index):
    A = dna_kmer_index.db.alphabet
    dna_kmer_index.db.initialize()

    dna_kmer_index.db.insert(A.parse('ATGCA', name='foo'))

    assert dna_kmer_index.kmers().next()[2] is None, \
        'score should be left None until explicitly calculated'

    dna_kmer_index.score_kmers()
    assert isinstance(dna_kmer_index.kmers().next()[2], float), \
        'score should be calculated when requested'

    # test only_missing: set the score of all observed kmers to 1, then
    # add a new sequence to DB, score kmers for only missing and make sure
    # the old kmers' scores remains 1. The kmer 'AAA' is 0.
    dna_kmer_index.db.insert(A.parse('AAA', name='bar'))
    dna_kmer_index.db.connection().cursor().execute("""
        UPDATE %s SET score = 1 WHERE kmer <> 0
    """ % dna_kmer_index.scores_table)
    dna_kmer_index.score_kmers(only_missing=True)
    scores = [score for kmer, h, score in dna_kmer_index.kmers() if kmer != 0]
    assert all(score == 1. for score in scores), \
        'score_kmers(only_missing=True) should work'

    dna_kmer_index.score_kmers(only_missing=False)
    scores = [score for kmer, h, score in dna_kmer_index.kmers() if kmer != 0]
    assert all(score != 1. for score in scores), \
        'score_kmers(only_missing=False) should work'

    with pytest.raises(StopIteration):
        # max_score should work
        dna_kmer_index.kmers(max_score=float('-inf')).next()


def test_bulk_index_kmers(dna_kmer_index):
    A = dna_kmer_index.db.alphabet
    dna_kmer_index.db.initialize()

    S = A.parse('ATGCA', name='foo')
    T = A.parse('ATGCC', name='bar')
    fasta = StringIO()
    write_fasta(fasta, [S, T])
    fasta.seek(0)

    S_rec, T_rec = dna_kmer_index.db.load_fasta(fasta)
    assert S_rec is not None and T_rec is not None
    assert dna_kmer_index.num_kmers() == 4
    assert dna_kmer_index.total_length_indexed() == len(S) + len(T)

    # find the occurences of 'ATG'
    S_id, T_id = S_rec.id, T_rec.id
    atg_int = dna_kmer_index.kmer_as_int((0, 3, 2))
    atg_hits = [(hits, score) for kmer, hits, score in dna_kmer_index.kmers()
                if kmer == atg_int]
    assert len(atg_hits) == 1 and atg_hits[0][0] == [(S_id, 0), (T_id, 0)]

    dna_kmer_index.score_kmers()
    atg_hits = [(hits, score) for kmer, hits, score in dna_kmer_index.kmers()
                if kmer == atg_int]
    atg_score = atg_hits[0][1]
    assert atg_score is not None and score > 0

    # ATG is the most common kmer, it most have the highest score:
    assert all(atg_score >= score for _, _, score in dna_kmer_index.kmers())

    # shouldn't do anything
    dna_kmer_index.score_kmers(only_missing=True)


def test_inf_score_kmers(dna_kmer_index):
    A = dna_kmer_index.db.alphabet
    dna_kmer_index.db.initialize()

    dna_kmer_index.db.insert(A.parse('A' * 1000, name='foo'))
    dna_kmer_index.score_kmers()
    _, _, score = dna_kmer_index.kmers().next()
    assert score == float('+inf')
