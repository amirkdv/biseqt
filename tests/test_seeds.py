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


def test_index_seeds():
    A = Alphabet('ACGT')
    with NamedTemporaryFile() as tmp:
        db = DB(tmp.name, A)
        kmer_index = KmerIndex(db, 5)
        seed_index = SeedIndex(kmer_index)

        kmer_index.db.initialize()
        records = []
        for i in range(3):
            seq = rand_seq(A, 100).to_named('#%d' % i)
            records.append(kmer_index.db.insert(seq))

        seed_index.index_seeds()
        scanned_sequences = kmer_index.scanned_sequences()
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
    'gap prob 1e-3': 1e-2,
    'gap prob 1e-1': 1e-1,
}
_substs = {
    'subst prob 1e-3': 1e-2,
    'subst prob 1e-1': 1e-1,
}
_words = {
    'k = 5': 5,
    'k = 8': 8,
}
_ids = product(_gaps.keys(), _substs.keys(), _words.keys())
_params = product(_gaps.values(), _substs.values(), _words.values())


@pytest.fixture(ids=[x for x in _ids], params=[x for x in _params])
def sequencing_sample(request):
    """Creates a random sequence of length 2000, and generates reads, with
    parameterized mutation probabilities, of length 500 starting at whole
    multiples of 250, i.e 0-500, 250-750, ..., 1500-2000. It is expected that
    successive reads have an overlap starting at position 250, upto 50
    positions of error should be tolerated in tests.

    Returns:
        tuple:
            full genome sequence, list of reads, gap probability, kmer index.
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

    with NamedTemporaryFile() as tmp:
        db = DB(tmp.name, A)
        kmer_index = KmerIndex(db, wordlen)
        return seq, reads, gap_prob, kmer_index


def test_score_seeds(sequencing_sample):
    seq, reads, gap_prob, kmer_index = sequencing_sample
    seed_index = SeedIndex(kmer_index)
    seed_index.db.initialize()

    num_reads = len(reads)
    avg_read_len = sum(len(read) for read in reads)/num_reads
    records = [seed_index.db.insert(read) for read in reads]
    assert sum(1 for _ in kmer_index.scanned_sequences()) == num_reads

    seed_index.index_seeds()
    seed_index.score_seeds(sensitivity=1-1e-4, gap_prob=gap_prob)

    for i in range(num_reads - 1):
        id0, id1 = records[i].id, records[i + 1].id
        best_diag = seed_index.highest_scoring_band(id0, id1)
        assert best_diag is not None, 'there should be a most dense diagonal'
        min_diag, max_diag = best_diag
        expected = avg_read_len/2
        err = avg_read_len/5
        assert expected >= min_diag - err and expected <= max_diag + err, \
            'Should find the most dense diagonal band upto %d' % err

        assert isinstance(next(seed_index.seeds(id0, id1))[0], Seed)
