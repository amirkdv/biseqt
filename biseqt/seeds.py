# -*- coding: utf-8 -*-
"""
.. wikisection:: overview
    :title: Statistical Seed Analysis

    The :mod:`biseqt.seeds` module provides tools for storing and analyzing
    matching segment pairs (aka seeds) between large numbers of sequences.

    >>> from biseqt.database import DB
    >>> from biseqt.sequence import Alphabet
    >>> from biseqt.kmers import KmerIndex
    >>> from biseqt.seeds import SeedIndex
    >>> A = Alphabet('ACGT')
    >>> db = DB('example.db', A)
    >>> seed_index = SeedIndex(KmerIndex(db))
    >>> db.initialize()
    >>> with open('example.fa') as f:
    ...     db.load_fasta(f)
    >>> seed_index.index_seeds()
    >>> seed_index.score_seeds(max_kmer_score=10)
    >>> seed_index.seeds(1, 2)  # yields all the seeds for sequences 1 and 2
    >>> diag_range = seed_index.highest_scoring_band(1, 2, min_band_score=10)
    >>> seed_index.seeds(1, 2, diag_range)  # only yields seeds in best band
"""
import numpy as np
from itertools import chain, combinations

from .kmers import KmerIndex, KmerDBWrapper
from .util import logging


class SeedIndex(KmerDBWrapper):
    """An index for seeds. Usage involves indexing seeds via
    :func:`index_seeds` and then FIXME

    Attributes:
        kmer_index (KmerIndex): The kmer index to operate on.
    """
    def __init__(self, S, T, path=':memory:', alphabet=None,
                 wordlen=None, log_level=logging.INFO):
        super(SeedIndex, self).__init__(path=path, wordlen=wordlen,
                                        alphabet=alphabet, log_level=log_level)
        self.self_comp = S == T
        self.S, self.T = S, T
        self.table = 'seeds_%s_%s' % (S.content_id[:8], T.content_id[:8])
        self.d0 = len(self.T) - 1
        if self._table_exists():
            self.log('seeds for %s and %s already indexed, skipping' %
                     (S.content_id[:8], T.content_id[:8]))
        else:
            self.log('Indexing seeds for %s and %s.' %
                     (S.content_id[:8], T.content_id[:8]))
            self._index_seeds()

    @classmethod
    def to_diagonal_coordinates(cls, i, j):
        d = i - j
        a = min(i, j)
        return d, a

    @classmethod
    def to_ij_coordinates(cls, d, a):
        i = a + max(d, 0)
        j = a - min(d, 0)
        return (i, j)

    def _table_exists(self):
        with self.connection() as conn:
            q = """
                SELECT name FROM sqlite_master
                WHERE type='table' AND name='%s';
            """ % self.table
            cursor = conn.cursor()
            cursor.execute(q)
            for name in cursor:
                return True
        return False

    # idempotent operation
    def _index_seeds(self):
        """Indexes all seeds and their diagonal positions. For each kmer with
        :math:`n` hits from *distinct* sequences :math:`n\\choose2` seeds
        are created. If a kmer occurs multiple times in a sequence seeds
        between different positions of the same sequence are not considered.
        If a minimum kmer score is required it is assumed that
        :func:`score_kmers <biseqt.kmers.KmerIndex.score_kmers>` has already
        been called.

        Keyword Args:
            max_kmer_score: The maximum score beyond which kmers are not
                considered for seeds. High scoring words are more likely to
                belong to repetitive regions (cf. :func:`score_kmers
                <biseqt.kmers.KmerIndex.score_kmers>`). Default is None in
                which case all kmers are considered.
        """
        with self.connection() as conn:
            conn.cursor().execute("""
                CREATE TABLE %s (
                  'd_' INTEGER,     -- zero-adjusted diagonal position
                  'a'  INTEGER      -- antidiagonal position
                );
            """ % self.table)

        kmer_index = KmerIndex(self.path, wordlen=self.wordlen,
                               alphabet=self.alphabet)
        kmer_index.index_kmers(self.S)
        if not self.self_comp:
            kmer_index.index_kmers(self.T)

        def _records():
            kmers = kmer_index.kmers()
            for kmer in kmers:
                hits = kmer_index.hits(kmer)
                if self.self_comp:
                    pairs = chain(combinations(hits, 2),
                                  [(x, x) for x in hits])
                else:
                    pairs = (((id0, pos0), (id1, pos1))
                             for (id0, pos0), (id1, pos1)
                             in combinations(hits, 2)
                             if id0 != id1)
                for (id0, pos0), (id1, pos1) in pairs:
                    d, a = self.to_diagonal_coordinates(pos0, pos1)
                    yield d + self.d0, a

        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.executemany(
                'INSERT INTO %s (d_, a) VALUES (?, ?)' % self.table,
                _records()
            )
            self.log('Creating SQL index for seeds table.')
            cursor.execute('CREATE INDEX %s_diagonal ON %s(d_);' %
                           (self.table, self.table))

        kmer_index.drop_data()

    def seed_count_by_d_(self):
        q = 'SELECT COUNT(a), d_ FROM %s GROUP BY d_' % self.table
        count_by_d_ = np.zeros(len(self.S) + len(self.T))
        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(q)
            for count, d_ in cursor:
                count_by_d_[d_] = count
        return count_by_d_

    def seed_in_band_count_by_a(self, center, radius):
        q = """
            SELECT COUNT(d_), a FROM %s
            WHERE d_ - ? BETWEEN ? AND ?
            GROUP BY a
        """ % self.table
        count_by_a = np.zeros(min(len(self.S), len(self.T)))
        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(q, (self.d0, center - radius, center + radius))
            for count, a in cursor:
                count_by_a[a] = count

        return count_by_a

    def seeds(self, d_center=None, d_radius=None):
        """Yields the :class:`seeds <Seed>` and their respective scores for a
        given pair of sequences and a given diagonal range (cf.
        :func:`highest_scoring_band`). Only seeds that are processed by
        :func:`score_seeds` are considered.

        Args:
            id0 (int): The identifier for the first (vertical) sequence.
            id1 (int): The idenfitier for the second (horizontal) sequence.

        Keyword Args:
            diag_range (tuple): The inclusive upper and lower bounds for the
                diagonals to consider. Default is None in which case all
                diagonals are considered.

        Yields:
            tuple:
                The :class:`Seed` object and the score of its diagonal as a
                ``float`` in descending order of score.
        """
        query = 'SELECT d_, a FROM %s' % self.table
        if d_center and d_radius:
            query += ' WHERE d_ - %d BETWEEN %d AND %d ' % \
                (self.d0, d_center - d_radius, d_center + d_radius)

        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(query)
            for d_, a in cursor:
                i, j = self.to_ij_coordinates(d_ - self.d0, a)
                yield (i, j)
                if self.self_comp and i != j:
                    yield (j, i)
