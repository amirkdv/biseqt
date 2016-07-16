# -*- coding: utf-8 -*-
import pytest
from itertools import product, combinations
from StringIO import StringIO

from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess
from biseqt.database import DB
from biseqt.kmers import KmerIndex
from biseqt.seeds import band_radius, SeedIndex
from biseqt.io import write_fasta


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
def seed_index():
    """Creates a database, a kmer index, and a seed index with word length 5
    stored in memory and returns the seed index. The database is populated with
    3 random sequences of length 100 and all kmers and seeds are indexed."""
    A = Alphabet('ACGT')
    num_seqs = 3
    seq_len = 100
    wordlen = 5

    db = DB(':memory:', A)
    seed_index = SeedIndex(KmerIndex(db, wordlen))
    seed_index.db.initialize()

    fasta = StringIO()
    seqs = (rand_seq(A, seq_len).to_named('#%d' % i) for i in range(num_seqs))
    write_fasta(fasta, seqs)
    fasta.seek(0)

    db.load_fasta(fasta)
    seed_index.index_seeds()
    return seed_index


def test_index_seeds(seed_index):
    seqpairs = combinations(seed_index.kmer_index.scanned_sequences(), 2)
    for (id0, len0), (id1, len1) in seqpairs:
        seeds = list(seed_index.seeds(id0, id1))
        assert len(seeds) > 0, 'all seeds should be indexed'
        for seed in seeds:
            assert seed.id0 == id0 and seed.id1 == id1
            diag = seed.pos0 - seed.pos1
            assert diag <= len0 - seed_index.wordlen and \
                diag >= -len1 + seed_index.wordlen


def test_count_seeds_on_diagonals(seed_index):
    seed_index.count_seeds_on_diagonals()
    with seed_index.db.connection() as conn:
        cursor = conn.cursor()
        cursor.execute("""
            SELECT COUNT(*) FROM %s WHERE count IS NOT NULL
        """ % seed_index.diagonals_table)
        assert cursor.next()[0] > 0, 'seeds on diagonals should be counted'
        cursor.execute("""
            SELECT COUNT(*) FROM %s WHERE count IS NOT NULL
        """ % seed_index.diagonals_table)


def test_cum_seed_count(seed_index):
    seqpairs = combinations(seed_index.kmer_index.scanned_sequences(), 2)
    for (id0, len0), (id1, len1) in seqpairs:
        count, diags = seed_index.cum_seed_count(id0, id1, len0, len1)
        assert all(x <= y for x, y in zip(count, count[1:])), \
            'cumulative seed count must be non-decreasing'
        jumps = [i - len1 for i in range(1, len(count))
                 if count[i] > count[i-1]]
        assert jumps == diags, \
            'for every diagonal with seeds there must be a jump in count'


# Parameters for fixture parameterization
_gaps = {
    'gap prob 1e-2': 1e-2,
}
_substs = {
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
    """Creates a random sequence, generates reads, with parameterized mutation
    probabilities, of equal length starting at whole multiples of half of read
    length. It is expected that successive reads have an overlap starting at
    their halfway position.

    Returns:
        tuple:
            A tuple containing the full genome, a list of reads, the gap
            probability and the seed index.
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

    db = DB(':memory:', A)
    kmer_index = KmerIndex(db, wordlen)
    seed_index = SeedIndex(kmer_index)
    seed_index.db.initialize()
    records = [db.insert(r) for r in reads]
    seed_index.score_diagonals()
    return seq, reads, records, gap_prob, seed_index


def test_score_diagonals(sequencing_sample):
    seq, reads, records, gap_prob, seed_index = sequencing_sample
    seed_index.index_seeds()
    seed_index.count_seeds_on_diagonals()
    seed_index.calculate_band_radii(sensitivity=1-1e-4, gap_prob=gap_prob)
    seed_index.score_diagonals()
    num_reads = len(reads)
    avg_read_len = sum(len(read) for read in reads)/num_reads
    for i in range(num_reads - 1):
        id0, id1 = records[i].id, records[i + 1].id
        impossible_diag, score = seed_index.highest_scoring_band(
            id0, id1, min_band_score=float('+inf')
        )
        assert score is None or score == float('+inf'), \
            'min_band_score should work'

        best_diag, _ = seed_index.highest_scoring_band(id0, id1)

        assert best_diag is not None, \
            'there should always be a most dense diagonal'
        min_diag, max_diag = best_diag
        expected = avg_read_len/2
        err = avg_read_len/5
        assert expected >= min_diag - err and expected <= max_diag + err, \
            'should find the most dense diagonal band upto %d' % err

        seeds = seed_index.seeds(id0, id1)
        assert all(seed.id0 == id0 and seed.id1 == id1 for seed in seeds), \
            'seeds should be properly retrieved'

        seeds = seed_index.seeds(id0, id1, diag_range=best_diag)
        diags = (seed.pos0 - seed.pos1 for seed in seeds)
        assert all(diag <= max_diag and diag >= min_diag for diag in diags), \
            'seeds within a diagonal range should be properly retrieved'
