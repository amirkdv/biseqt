# -*- coding: utf-8 -*-
"""
.. wikisection:: overview
    :title: (5) Seeds (exactly matching kmers)

    The :mod:`biseqt.seeds` module provides tools for storing and analyzing
    exactly matching kmers, aka seeds.

    >>> from biseqt.seeds import SeedIndex
    >>> from biseqt.sequence import Sequence, Alphabet
    >>> A = Alphabet('ACGT')
    >>> S, T = A.parse('TAAGCGT'), A.parse('GGCGTAA')
    >>> seed_index = SeedIndex(S, T, path=':memory:', wordlen=3, alphabet=A)
    >>> list(seed_index.seeds())
    [(4, 2), (3, 1), (0, 4)]
"""
import numpy as np
from scipy.spatial import cKDTree

from .kmers import KmerIndex, KmerDBWrapper


class SeedIndex(KmerDBWrapper):
    """An index for seeds in diagonal coordinates.

    Attributes:
        S (biseqt.sequence.Sequence): The 1st sequence.
        T (biseqt.sequence.Sequence): The 2nd sequence.
        cache (KmerCache): optional :class:`KmerCache` object to use for
            retrieving integer representations of sequences.
    """
    def __init__(self, *seqs, **kw):
        assert(len(seqs)) > 1
        name = '_'.join(S.content_id[:8] for S in seqs)
        super(SeedIndex, self).__init__(name=name, **kw)
        self.kmer_cache = kw.get('kmer_cache', None)
        self.self_comp = len(seqs) == 2 and seqs[0] == seqs[1]
        self.seqs = seqs
        self.d0s = [len(T) - 1 for T in seqs[1:]]

        if self._table_exists():
            self.log('Seeds for %s already indexed, skipping' % name)
        else:
            self.log('Indexing seeds for %s.' % name)
            self._index_seeds()

    @property
    def seeds_table(self):
        """The seeds table name ``seeds_[name]``, cf.
        :attr:`KmerDBWrapper.name`."""
        return 'seeds_' + self.name

    @classmethod
    def to_diagonal_coordinates(cls, *idxs):
        """Convert standard coordinates to diagonal coordinates via:

        .. math::
            \\begin{aligned}
                d & = i - j \\\\
                a & = \min(i, j)
            \\end{aligned}
        """
        ds = [idxs[0] - idxs[k] for k in range(1, len(idxs))]
        a = min(idxs)
        return ds, a

    @classmethod
    def to_ij_coordinates(cls, ds, a):
        """Convert diagonal coordinates to standard coordinates:

        .. math::
            \\begin{aligned}
                i & = a + max(d, 0) \\\\
                j & = a - min(d, 0)
            \\end{aligned}
        """
        idxs_ = [0] + [-d for d in ds]
        diff = a - min(idxs_)
        idxs = [diff + idx for idx in idxs_]
        return tuple(idxs)
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
        with self.connection() as conn:
            conn.cursor().execute("""
                CREATE TABLE %s (
                  'ds' VARCHAR,     -- comma separated diagonal positions
                  'a'  INTEGER      -- antidiagonal position
                );
            """ % self.seeds_table)

        kmer_index_name = '%d_%s' % (self.wordlen, self.name)
        kmer_index = KmerIndex(path=self.path, name=kmer_index_name,
                               wordlen=self.wordlen, alphabet=self.alphabet,
                               log_level=self.log_level,
                               kmer_cache=self.kmer_cache)
        if self.self_comp:
            kmer_index.index_kmers(self.seqs[0])
        else:
            for seq in self.seqs:
                kmer_index.index_kmers(seq)

        kmers = kmer_index.kmers()

        def _records():
            for kmer in kmers:
                hits = kmer_index.hits(kmer)
                from itertools import groupby, product
                hits = {seqid: [c[1] for c in seq_hits]
                        for seqid, seq_hits in groupby(hits,
                                                       key=lambda c: c[0])}
                # only consider kmers present in all sequences
                if len(hits) < len(self.seqs):
                    continue
                # FIXME deal with self_comp: let 1 argument mean self_comp not
                # two identical sequences.
                for idxs in product(*hits.values()):
                    # NOTE we're storing d values and not d_
                    ds, a = self.to_diagonal_coordinates(*idxs)
                    yield ','.join(str(d) for d in ds), a

        with self.connection() as conn:
            cursor = conn.cursor()
            query = 'INSERT INTO %s (ds, a) VALUES (?, ?)' % self.seeds_table
            cursor.executemany(query, _records())
            self.log('Creating SQL index for table %s.' % self.seeds_table)
            cursor.execute('CREATE INDEX %s_diagonal ON %s(ds);' %
                           (self.seeds_table, self.seeds_table))

    def seed_count_by_d_(self):
        """Number of seeds in each diagonal position. Diagonals are
        ordered from smallest :math:`-|T|` to largest :math:`|S|`.
        """
        q = 'SELECT COUNT(a), d_ FROM %s GROUP BY d_' % self.seeds_table
        count_by_d_ = np.zeros(len(self.S) + len(self.T) - 1)
        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(q)
            for count, d_ in cursor:
                count_by_d_[d_] = count
        return count_by_d_

    def seed_count_by_a(self, d_min, d_max):
        """Number of seeds at each antidiagonal position in a givan diagonal
        band.
        """
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

    def count_all_seed_neighbors(self, d_radius, a_radius):
        """For each seed finds the number of seeds in its neighborhood defined
        by:

        .. math::

            V_{(d, a)} = \\{(d', a'): |d - d'| < r_d,  |a - a'| < r_a \\}

        This is done using a Quad-Tree in ``O(m lg m)`` time where m is the
        number of seeds.

        Returns:
            list: tuples ``((d, a), num_neighs)``
        """
        # normalize the two diameters so we can use a standard L∞ neighborhood.
        # typically a_diam is larger, so scale up d values proportionally
        d_coeff = 1. * a_radius / d_radius
        radius = a_radius

        all_seeds = list(self.seeds())
        if not all_seeds:
            return []
        all_seeds_scaled = np.array([[d * d_coeff for d in ds] + [a]
                                     for ds, a in all_seeds])
        quad_tree = cKDTree(all_seeds_scaled)
        all_neighs = quad_tree.query_ball_tree(quad_tree, radius,
                                               p=float('inf'))
        # neighs[i] is the indices of the neighbors of all_seeds[i]; this
        # always contains the seed itself (i.e always: i in neighs[i])
        neigh_counts = [len(neighs) - 1 for neighs in all_neighs]
        return zip(all_seeds, neigh_counts)

    def seeds(self, d_band=None, exclude_trivial=False):
        """Yields all seeds, optionally those within a diagonal band.

        Keyword Args:
            d_band (tuple|None):
                If specified a ``(d_min, d_max)`` tuple restricting the seed
                count to a diagonal band.

        Yields:
            tuple:
                seeds coordinates :math:`(i, j)`.
        """
        query = 'SELECT ds, a FROM %s' % self.seeds_table
        if d_band is not None:
            raise NotImplementedError
            assert len(d_band) == 2, 'need a 2-tuple for diagonal band'
            d_min, d_max = d_band
            query += ' WHERE d_ - %d BETWEEN %d AND %d ' % \
                (self.d0, d_min, d_max)

        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(query)
            for ds, a in cursor:
                ds = [int(i) for i in ds.split(',')]
                # FIXME deal with self_comp
                yield ds, a

    # FIXME change this to d_min, d_max
    def seed_count(self, d_band=None):
        """Counts the number of seeds either in the whole table or in the
        specified diagonal band.

        Args:
            d_band (tuple|None):
                If specified a ``(d_min, d_max)`` tuple restricting the seed
                count to a diagonal band.

        Returns:
            int: Number of seeds found in the entire table or in the specified
            diagonal band.
        """
        query = 'SELECT COUNT(*) FROM %s' % self.seeds_table
        if d_band is not None:
            assert len(d_band) == 2, 'need a 2-tuple for diagonal band'
            d_min, d_max = d_band
            query += ' WHERE d_ - %d BETWEEN %d AND %d ' % \
                (self.d0, d_min, d_max)

        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(query)
            for row in cursor:
                return row[0]
            return 0
