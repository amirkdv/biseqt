# -*- coding: utf-8 -*-
import pytest
from tempfile import NamedTemporaryFile
from itertools import product, combinations

from biseqt.sequence import Alphabet
from biseqt.random import rand_seq, MutationProcess
from biseqt.database import DB
from biseqt.kmers import KmerIndex
from biseqt.seeds import band_radius, SeedIndex, Seed


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


@pytest.fixture
def seed_index_gen():
    """Returns a context manager that creates a database, a kmer index, and a
    seed index with word length 5 stored in a temporary file. The database is
    populated with 3 sequences of length 100 and all kmers and seeds are
    indexed.
    """
    A = Alphabet('ACGT')
    num_seqs = 3
    seq_len = 100
    wordlen = 5

    class context(object):
        def __enter__(self):
            self.tmp = NamedTemporaryFile()
            db = DB(self.tmp.name, A)
            seed_index = SeedIndex(KmerIndex(db, wordlen))
            seed_index.db.initialize()
            for i in range(num_seqs):
                seq = rand_seq(A, seq_len).to_named('#%d' % i)
                db.insert(seq)
            seed_index.index_seeds()
            return seed_index

        def __exit__(self, *args):
            self.tmp.close()

    return context


def test_index_seeds(seed_index_gen):
    with seed_index_gen() as seed_index:
        scanned_sequences = seed_index.kmer_index.scanned_sequences()
        for (id0, len0), (id1, len1) in combinations(scanned_sequences, 2):
            for seed, score in seed_index.seeds(id0, id1):
                assert seed.id0 == id0 and seed.id1 == id1
                diag = seed.pos0 - seed.pos1
                assert diag <= len0 - seed_index.wordlen and \
                    diag >= -len1 + seed_index.wordlen


def test_seed_count(seed_index_gen):
    with seed_index_gen() as seed_index:
        scanned_sequences = seed_index.kmer_index.scanned_sequences()
        for (id0, len0), (id1, len1) in combinations(scanned_sequences, 2):
            count, diags = seed_index.cum_seed_count(id0, id1, len0, len1)
            assert all(x <= y for x, y in zip(count, count[1:])), \
                'cumulative seed count must be non-decreasing'
            jumps = [i - len1 for i in range(1, len(count))
                     if count[i] > count[i-1]]
            assert jumps == diags, \
                'for every diagonal with seeds there must be a jump in count'


# Parameters for fixture parameterization
_gaps = {
    'gap prob 1e-3': 1e-3,
    'gap prob 1e-2': 1e-2,
}
_substs = {
    'subst prob 1e-3': 1e-3,
    'subst prob 1e-2': 1e-2,
}
_words = {
    'k = 5': 5,
    'k = 8': 8,
}
_ids = product(_gaps.keys(), _substs.keys(), _words.keys())
_params = product(_gaps.values(), _substs.values(), _words.values())


@pytest.fixture(ids=[x for x in _ids], params=[x for x in _params])
def sequencing_sample(request):
    """Creates a context manager that creates a random sequence, generates
    reads, with parameterized mutation probabilities, of equal length starting
    at whole multiples of half of read length. It is expected that successive
    reads have an overlap starting at their halfway position.

    Upon invocation a tuple is returned containing: full genome sequence, a
    list of lossy reads, inserted records, gap probability, and the seed index.
    """
    A = Alphabet('ACGT')
    gap_prob, subst_prob, wordlen = request.param
    seq_len, read_len = 2000, 500
    seq = rand_seq(A, seq_len).to_named('genome')
    mutation_process = MutationProcess(A, subst_probs=subst_prob,
                                       go_prob=gap_prob, ge_prob=gap_prob)
    reads = []
    for i in range(0, seq_len - read_len, int(read_len/2)):
        read, _ = mutation_process.mutate(seq[i: i + read_len])
        reads += [read.to_named('read#%d' % i)]

    class context(object):
        def __enter__(self):
            self.tmp = NamedTemporaryFile()
            db = DB(self.tmp.name, A)
            kmer_index = KmerIndex(db, wordlen)
            seed_index = SeedIndex(kmer_index)
            seed_index.db.initialize()
            records = [db.insert(read) for read in reads]
            seed_index.index_seeds()
            seed_index.score_seeds(sensitivity=1-1e-4, gap_prob=gap_prob)
            return seq, reads, records, gap_prob, seed_index

        def __exit__(self, *args):
            self.tmp.close()

    return context


def _test_score_seeds(sequencing_sample):
    with sequencing_sample() as sample:
        seq, reads, records, gap_prob, seed_index = sample
        num_reads = len(reads)
        scanned_sequences = seed_index.kmer_index.scanned_sequences()
        assert sum(1 for _ in scanned_sequences) == num_reads

        with seed_index.db.connect() as conn:
            cursor = conn.cursor()
            query = 'SELECT COUNT(*) FROM seeds_%d WHERE score IS NULL' % \
                seed_index.wordlen
            cursor.execute(query)
            assert next(cursor)[0] == 0


def test_highest_scoring_seeds(sequencing_sample):
    with sequencing_sample() as sample:
        seq, reads, records, gap_prob, seed_index = sample
        num_reads = len(reads)
        avg_read_len = sum(len(read) for read in reads)/num_reads
        for i in range(num_reads - 1):
            id0, id1 = records[i].id, records[i + 1].id
            best_diag = seed_index.highest_scoring_band(id0, id1)
            assert best_diag is not None, \
                'there should always be a most dense diagonal'
            min_diag, max_diag = best_diag
            expected = avg_read_len/2
            err = avg_read_len/5
            assert expected >= min_diag - err and expected <= max_diag + err, \
                'Should find the most dense diagonal band upto %d' % err

            assert isinstance(next(seed_index.seeds(id0, id1))[0], Seed)
