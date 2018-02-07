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

from itertools import combinations

from .kmers import KmerIndex, KmerDBWrapper
from .util import logging


class SeedIndex(KmerDBWrapper):
    """An index for seeds. Usage involves indexing seeds via
    :func:`index_seeds` and then FIXME

    Attributes:
        kmer_index (KmerIndex): The kmer index to operate on.
    """
    def __init__(self, path=':memory:', alphabet=None,
                 wordlen=None, log_level=logging.INFO):
        super(SeedIndex, self).__init__(path=path, wordlen=wordlen,
                                        alphabet=alphabet, log_level=log_level)

    def table_name(self, S, T):
        return 'seeds_%s_%s' % (S.content_id[:8], T.content_id[:8])

    def _table_exists(self, table):
        with self.connection() as conn:
            q = """
                SELECT name FROM sqlite_master
                WHERE type='table' AND name='%s';
            """ % table
            cursor = conn.cursor()
            cursor.execute(q)
            for name in cursor:
                return True
        return False

    # idempotent operation
    def index_seeds(self, S, T):
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
        table = self.table_name(S, T)
        if self._table_exists(table):
            self.log('seeds for %s and %s already indexed, skipping' %
                     (S.content_id[:8], T.content_id[:8]))
            return

        with self.connection() as conn:
            conn.cursor().execute("""
                CREATE TABLE %s (
                  'd_' INTEGER,     -- zero-adjusted diagonal position
                  'a'  INTEGER      -- antidiagonal position
                );
            """ % table)

        # FIXME the kmers table shouldn't be reused, use a separate file but
        # stay in memory if we are :memory:
        kmer_index = KmerIndex(self.path, wordlen=self.wordlen,
                               alphabet=self.alphabet)
        self_comp = S == T
        kmer_index.index_kmers(S)
        if not self_comp:
            kmer_index.index_kmers(T)

        def _records():
            kmers = kmer_index.kmers()
            for kmer in kmers:
                hits = kmer_index.hits(kmer)
                # FIXME check integrity of id0 and id1? we should behave
                # differently if we are expected to work under S vs S
                # comparison
                d0 = len(T) - 1
                for (id0, pos0), (id1, pos1) in combinations(hits, 2):
                    if not self_comp and id0 == id1:
                        continue
                    d = pos0 - pos1
                    d_ = d + d0
                    a = min(pos0, pos1)
                    yield d_, a

        self.log('Indexing seeds.')
        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.executemany(
                'INSERT INTO %s (d_, a) VALUES (?, ?)' % table,
                _records()
            )
            self.log('Creating SQL index for seeds table.')
            cursor.execute('CREATE INDEX %s_diagonal ON %s(d_);' %
                           (table, table))

    def seed_count_by_d_(self, S, T):
        self.index_seeds(S, T)
        q = 'SELECT COUNT(a), d_ FROM %s GROUP BY d_' % self.table_name(S, T)
        count_by_d_ = np.zeros(len(S) + len(T))
        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(q)
            for count, d_ in cursor:
                count_by_d_[d_] = count
        return count_by_d_

    def seed_in_band_count_by_a(self, S, T, center, radius):
        self.index_seeds(self, S, T)
        q = """
            SELECT COUNT(d_), a FROM %s
            WHERE d_ - ? BETWEEN ? AND ?
            GROUP BY a
        """ % self.table_name(S, T)
        d0 = len(T) - 1
        count_by_a = np.zeros(min(len(S), len(T)))
        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(q, (d0, center - radius, center + radius))
            for count, a in cursor:
                count_by_a[a] = count

        return count_by_a

    def seeds(self, S, T, d_center=None, d_radius=None):
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
        self.index_seeds(S, T)
        d0 = len(T) - 1
        query = 'SELECT d_, a FROM %s' % self.table_name(S, T)
        if d_center and d_radius:
            query += ' WHERE d_ - %d BETWEEN %d AND %d ' % \
                (d0, d_center - d_radius, d_center + d_radius)

        self_comp = S == T
        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(query)
            for d_, a in cursor:
                d = d_ - d0
                i, j = (a + max(d, 0), a - min(d, 0))
                yield (i, j)
                if self_comp and i != j:
                    yield (j, i)
