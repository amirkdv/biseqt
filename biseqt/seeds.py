# -*- coding: utf-8 -*-
"""
.. wikisection:: overview
    :title: Statistical Seed Analysis

    The :mod:`biseqt.seeds` module provides tools for storing and analyzing
    matching segment pairs (aka seeds) between large numbers of sequences.

    FIXME
"""
import numpy as np
from itertools import chain, combinations

from .kmers import KmerIndex, KmerDBWrapper


class SeedIndex(KmerDBWrapper):
    """An index for seeds in diagonal coordinates.

    Attributes:
        S (biseqt.sequence.Sequence): The 1st sequence.
        T (biseqt.sequence.Sequence): The 2nd sequence.
    """
    def __init__(self, S, T, kmer_cache=None, **kw):
        name = '%s_%s' % (S.content_id[:8], T.content_id[:8])
        super(SeedIndex, self).__init__(name=name, **kw)
        self.kmer_cache = kmer_cache
        self.self_comp = S == T
        self.S, self.T = S, T
        self.d0 = len(self.T) - 1
        if self._table_exists():
            self.log('seeds for %s and %s already indexed, skipping' %
                     (S.content_id[:8], T.content_id[:8]))
        else:
            self.log('Indexing seeds for %s and %s.' %
                     (S.content_id[:8], T.content_id[:8]))
            self._index_seeds()

    @property
    def seeds_table(self):
        return 'seeds_' + self.name

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
            """ % self.seeds_table
            cursor = conn.cursor()
            cursor.execute(q)
            for name in cursor:
                return True
        return False

    # idempotent operation
    def _index_seeds(self):
        """TODO
        """
        with self.connection() as conn:
            conn.cursor().execute("""
                CREATE TABLE %s (
                  'd_' INTEGER,     -- zero-adjusted diagonal position
                  'a'  INTEGER      -- antidiagonal position
                );
            """ % self.seeds_table)

        kmer_index_name = '%d_%s' % (self.wordlen, self.name)
        kmer_index = KmerIndex(path=self.path, name=kmer_index_name,
                               wordlen=self.wordlen, alphabet=self.alphabet,
                               log_level=self.log_level, kmer_cache=self.kmer_cache)
        kmer_index.index_kmers(self.S)
        if not self.self_comp:
            kmer_index.index_kmers(self.T)

        kmers = kmer_index.kmers()

        def _records():
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
                'INSERT INTO %s (d_, a) VALUES (?, ?)' % self.seeds_table,
                _records()
            )
            self.log('Creating SQL index for table %s.' % self.seeds_table)
            cursor.execute('CREATE INDEX %s_diagonal ON %s(d_);' %
                           (self.seeds_table, self.seeds_table))

    def seed_count_by_d_(self):
        q = 'SELECT COUNT(a), d_ FROM %s GROUP BY d_' % self.seeds_table
        count_by_d_ = np.zeros(len(self.S) + len(self.T) - 1)
        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(q)
            for count, d_ in cursor:
                count_by_d_[d_] = count
        return count_by_d_

    def seed_count_by_a(self, d_min, d_max):
        q = """
            SELECT COUNT(d_), a FROM %s
            WHERE d_ - ? BETWEEN ? AND ?
            GROUP BY a
        """ % self.seeds_table
        count_by_a = np.zeros(min(len(self.S), len(self.T)))
        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(q, (self.d0, d_min, d_max))
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
        query = 'SELECT d_, a FROM %s' % self.seeds_table
        if d_center is not None and d_radius is not None:
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

    def seed_count(self, d_center=None, d_radius=None):
        """Counts the number of seeds either in the whole table or in the
        specified diagonal band.

        Args:
            d_center (int|None):
                If specified, the diagonal number of the center of band.
            d_radius (int|None):
                If specified, the radius of the band.

        Returns:
            int: Number of seeds found in the entire table or in the specified
            diagonal band.
        """
        query = 'SELECT COUNT(*) FROM %s' % self.seeds_table
        if d_center is not None and d_radius is not None:
            query += ' WHERE d_ - %d BETWEEN %d AND %d ' % \
                (self.d0, d_center - d_radius, d_center + d_radius)

        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(query)
            for row in cursor:
                return row[0]
            return 0
