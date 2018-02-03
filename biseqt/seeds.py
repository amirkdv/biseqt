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
    >>> seed_index.index_seeds(max_kmer_score=10)
    >>> seed_index.seeds(1, 2)  # yields all the seeds for sequences 1 and 2
    >>> diag_range = seed_index.highest_scoring_band(1, 2, min_band_score=10)
    >>> seed_index.seeds(1, 2, diag_range)  # only yields seeds in best band
"""

from collections import namedtuple
from itertools import combinations, product
from math import sqrt
import apsw

from .kmers import KmerIndex
from .stochastics import binomial_to_normal, band_radius_calculator


class Seed(namedtuple('Seed', ['id0', 'id1', 'pos0', 'pos1', 'length'])):
    """A ``namedtuple`` representing a matching segment between two sequences.
    A seed is uniquely defined by two sequence identifiers, two integers
    representing the starting positions in each sequence, and the length. For
    instance, the two sequences ``AAACTG`` and ``GCAAACA`` share only one seed
    of length 4, namely ``AAAC``, with starting positions 0 and 2,
    respectively.

    Attributes:
        id0 (int): Identifier of the first sequence id, as it appears in the
            ``sequence`` table.
        id1 (int): Identifier of the second sequence id, as it appears in the
            ``sequence`` table.
        pos0 (int): Starting position of the seed in sequence id0.
        pos1 (int): Starting position of the seed in sequence id1.
        length (int): Length of the matching segment.
    """


class SeedIndex(object):
    """An index for :class:`seeds <Seed>`. Usage involves indexing seeds via
    :func:`index_seeds` and then scoring them via :func:`score_seeds` after
    which for each sequence pair diagonal ranges can be queried for seeds via
    :func:`seeds`.

    Attributes:
        kmer_index (KmerIndex): The kmer index to operate on.
        db (database.DB): The database inherited from :attr:`kmer_index`.
        wordlen (int): The word length inherited from :attr:`kmer_index`.
        seeds_table (str): ``seeds_N`` contains all seeds for each pair of
            scanned sequences (``N`` is :attr:`wordlen`).
        diagonals_table (str): ``diagonals_N`` contains scores and radii for
            each diagonal of the edit graph of each pair of scanned sequences
            (``N`` is :attr:`wordlen`).
    """
    def __init__(self, kmer_index):
        assert isinstance(kmer_index, KmerIndex)
        self.kmer_index = kmer_index
        self.wordlen = self.kmer_index.wordlen
        self.db = self.kmer_index.db
        self.seeds_table = 'seeds_%d' % self.wordlen
        self.kmers_table = 'kmers_%d' % self.wordlen
        self.diagonals_table = 'diagonals_%d' % self.wordlen
        self.connection(reset=True)

    _init_script = """
        CREATE TABLE %s (
          'kmer'  INTEGER,              -- The kmer in integer representation.
          'seq'   INTEGER,              -- REFERENCES sequence(id)
                                        -- but do not declare it to save time
                                        -- on checking referential integrity.
          'pos' INTEGER                 -- the position of kmer in sequence.
        );
        CREATE TABLE %s ( -- seeds table
          'id0'   INTEGER,              -- REFERENCES sequence(id)
                                        -- but do not declare it to save time
                                        -- on checking referential integrity.
          'id1'   INTEGER,              -- REFERENCES kmers_N_hits(id), same.
          'pos0'  INTEGER,              -- starting position in sequence 0.
          'pos1'  INTEGER               -- starting position in sequence 1.
        );
        CREATE TABLE %s ( -- diagonals table
          'id0'   INTEGER,              -- REFERENCES sequence(id)
                                        -- but do not declare it to save time
                                        -- on checking referential integrity.
          'id1'   INTEGER,              -- REFERENCES kmers_N_hits(id), same.
          'diag'  INTEGER,              -- a single diagonal in the edit graph.
          'count' INTEGER DEFAULT NULL, -- number of seeds in this diagonal
          'radius'INTEGER DEFAULT NULL, -- the band radius for the optimal band
                                        -- with this diagonal as center.
          'score' REAL DEFAULT NULL     -- the log(p-value) for the number of
                                        -- seeds on the band centered at this
                                        -- diagonal.
        );
    """

    def initialize(self, conn):
        """Event handler for "db-initialized" (cf. :attr:`DB.events
        <biseqt.database.DB.events>`). Creates two tables:

        * :attr:`seeds_table`,
        * :attr:`diagonals_table`,

        Args:
            conn (sqlite3.Connection): An open connection to operate on.
        """
        conn.cursor().execute(
            self._init_script % (self.kmers_table, self.seeds_table, self.diagonals_table)
        )

    # put the initialization script in the docs
    initialize.__doc__ += '\n\n\t.. code-block:: sql\n\t%s\n' % \
                          '\n\t'.join(_init_script.split('\n'))

    def connection(self, reset=False):
        if reset or self._connection is None:
            self._connection = apsw.Connection(':memory:')
        return self._connection

    # FIXME not needed
    def create_sql_index(self, table=None):
        """Creates a ``(id0, id1)`` SQL index over a given table.

        Args:
            table (str): Either of :attr:`seeds_table` or
                :attr:`diagonals_table`.
        """
        assert table in [self.seeds_table, self.diagonals_table]
        with self.connection() as conn:
            conn.cursor().execute("""
                CREATE INDEX IF NOT EXISTS %s_seqpair ON %s (id0, id1)
            """ % (table, table))

    def log(self, *args, **kwargs):
        """Wraps :func:`log <biseqt.database.DB.log>` of :attr:`db`."""
        self.db.log(*args, **kwargs)

    def index_kmers(self, conn, ids=None):
        def _records():
            for seq, kmers in self.kmer_index.kmers(ids=ids):
                for pos, kmer in enumerate(kmers):
                    yield kmer, seq, pos

        conn.cursor().execute('DROP INDEX IF EXISTS %s_kmer;' % self.kmers_table)
        q = 'INSERT INTO %s (kmer, seq, pos) VALUES (?, ?, ?)' % self.kmers_table
        conn.cursor().executemany(q, _records())
        conn.cursor().execute('CREATE INDEX %s_kmer ON %s (kmer);' % (self.kmers_table, self.kmers_table))

    def index_seeds(self, ids=None, max_kmer_score=None):
        """Indexes all seeds and their diagonal positions. For each kmer with
        :math:`n` hits from *distinct* sequences :math:`n\\choose2` seeds
        are created. If a kmer occurs multiple times in a sequence seeds
        between different positions of the same sequence are not considered.
        If a minimum kmer score is required it is assumed that
        :func:`score_kmers <biseqt.kmers.KmerIndex.score_kmers>` has already
        been called.

        Keyword Args:
            max_kmer_score (float): The maximum score beyond which kmers are
                not considered for seeds. High scoring words are more likely to
                belong to repetitive regions (cf. :func:`score_kmers
                <biseqt.kmers.KmerIndex.score_kmers>`). Default is None in
                which case all kmers are considered.
        """
        def _records(conn):
            q = """
                SELECT GROUP_CONCAT(seq), GROUP_CONCAT(pos)
                FROM %s
                GROUP BY kmer HAVING COUNT(*) > 1
                ORDER BY seq ASC
            """ % self.kmers_table
            cursor = conn.cursor()
            cursor.execute(q)
            for seqids, positions in cursor:
                seqids = eval('[' + seqids + ']')
                positions = eval('[' + positions + ']')
                assert len(seqids) == len(positions)
                hits = [(seqids[i], positions[i]) for i in range(len(seqids))]
                for (id0, pos0), (id1, pos1) in combinations(hits, 2):
                    if id0 == id1:
                        continue
                    yield id0, id1, pos0, pos1

        with self.connection() as conn:
            self.initialize(conn)
            self.index_kmers(conn, ids)
            conn.cursor().executemany(
                """INSERT INTO %s (id0, id1, pos0, pos1) VALUES (?, ?, ?, ?)
                """ % self.seeds_table,
                _records(conn)
            )

    def cum_seed_count(self, id0, id1, len0, len1):
        """Returns the cumulative count of seeds per diagonals. It is assumed
        that :attr:`diagonals table <diagonals_table>` is populated with seed
        counts on each diagonal by :func:`count_seeds_on_diagonals`.

        Args:
            id0 (int): The identifier for the first (vertical) sequence.
            id1 (int): The idenfitier for the second (horizontal) sequence.
            len0 (int): The length of the first sequence.
            len1 (int): The length of the second sequence.

        Returns:
            tuple:
                A tuple of a list of length ``len0 + len1 + 1`` containing the
                cumulative count of seeds upto each diagonal from ``-len1`` to
                ``+len0``, and a list of tuples containing the diagonals and
                their respective radii as already set in the :attr:`digonals
                table <diagonals_table>`.
        """
        query = """
            SELECT diag, radius, count FROM %s
            WHERE id0 = ? AND id1 = ? ORDER BY diag ASC
        """ % self.diagonals_table

        diags = []
        seed_count = [0] * (len1 + len0 + 1)
        prev_diag = 0
        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(query, (id0, id1))
            for diag, radius, count in cursor:
                diags.append((diag, radius))

                diag += len1
                for _diag in range(prev_diag, diag + 1):
                    seed_count[_diag] = seed_count[prev_diag]
                seed_count[diag] += count
                prev_diag = diag

        for _diag in range(prev_diag, len(seed_count)):
            seed_count[_diag] = seed_count[prev_diag]

        return seed_count, diags

    def count_seeds_on_diagonals(self):
        """Populates the :attr:`diagonals table <diagonals_table>` form the
        contents of the :attr:`seeds table <seeds_table>`. Radius and score
        information are left blank."""
        self.create_sql_index(self.seeds_table)
        with self.connection() as conn:
            conn.cursor().execute("""
                INSERT INTO %s (id0, id1, diag, count)
                SELECT id0, id1, pos0 - pos1, COUNT(*) FROM %s
                GROUP BY id0, id1, pos0 - pos1
            """ % (self.diagonals_table, self.seeds_table))

    def calculate_band_radii(self, gap_prob=None, sensitivity=None):
        """Calculates the diagonal band radius for each entry in the
        :attr:`diagonals table <diagonals_table>` using :func:`band_radius()
        <biseqt.stochastics.band_radius>`.

        Keyword Args:
            gap_prob (float): Passed as is, mut be strictly between 0 and 1.
            sensitivity (float): Passed as is, must be strictly between 0 and
                1.
        """
        radius = band_radius_calculator(gap_prob=gap_prob,
                                        sensitivity=sensitivity)
        select = """
            SELECT diag FROM %s WHERE id0 = ? AND id1 = ?
        """ % self.diagonals_table

        update = """
            UPDATE %s SET radius = ? WHERE id0 = ? AND id1 = ? AND diag = ?
        """ % self.diagonals_table

        # each entry is a pair of (id, len) tuples
        # FIXME use the sequence table and the lengths from attrs.
        scanned_seqs = list(self.kmer_index.scanned_sequences())
        seqpairs = combinations(scanned_seqs, 2)

        self.create_sql_index(self.diagonals_table)

        num_pairs = (len(scanned_seqs) * (len(scanned_seqs) - 1)) / 2
        with self.connection() as conn:
            cursor = conn.cursor()
            for (id0, len0), (id1, len1) in seqpairs:
                cursor = conn.cursor()
                diags = [x[0] for x in cursor.execute(select, (id0, id1))]
                recs = ((radius(len0, len1, diag), id0, id1, diag)
                        for diag in diags)
                conn.cursor().executemany(update, recs)

    def score_diagonals(self, ids, gap_prob=None, sensitivity=None):
        """Scores all diagonal bands containing seeds for each pair of
        sequences. It is assumed that all kmers are already indexed and scored
        (cf. :func:`biseqt.kmers.KmerIndex.score_kmers`). First, all seeds are
        indexed via :func:`index_seeds`. Then band radii for each diagonal is
        calculated by :func:`calculate_band_radii`. The score for a band is
        then calculated as the z-score of the observed number of seeds in it
        under the null hypothesis that seeds occur randomly throughout the
        dynamic programming table.

        Keyword Args:
            max_kmer_score (float): Passed to :func:`index_seeds`.
            gap_prob (float): Passed to :func:`calculate_band_radii`.
            sensitivity (float): Passed to :func:`calculate_band_radii`.

        For a fixed pair of sequences with lengths :math:`l_0, l_1`
        a diagonal band centered at :math:`d` and with radius :math:`r`, the
        dimensions of the dynamic programming table is :math:`[0, l_0] \\times
        [0, l_1]`. Under the null hypothesis, namely that the two sequences are
        not related, the number :math:`n` of observed seeds in the diagonal
        band :math:`[d-r, d+r]` is distributed according to a binomial
        distribution:

        .. math::
            n \\sim B(A, p)

        where :math:`A` is the area of the diagonal band:

        .. math::

            A \\simeq 2r\\sqrt{(l_0-|d|)^2 + (l_1-|d|)^2}

        and :math:`p` is the probability that two arbitrary kmers are
        identical, that is:

        .. math::

            p = \\left(\\frac{1}{|\\Sigma|}\\right)^k

        To score a band, the normal approximation of the binomial is considered
        (cf. :func:`binomial_to_normal()
        <biseqt.stochastics.binomial_to_normal>`) and the z-score of the
        observed number of seeds is reported, namely:

        .. math::

                \\mathrm{score} = \\frac{n-Ap}{\\sqrt{Ap(1-p)}}

        A high scoring diagonal band is less likely to be accidental and more
        likely to indicate an overlap between the two sequences.
        """
        #self.kmer_index.score_kmers() FIXME make this run pre-comparison by mapper
        self.index_seeds(ids)
        self.count_seeds_on_diagonals()
        self.calculate_band_radii(gap_prob=gap_prob, sensitivity=sensitivity)
        query = """
            UPDATE %s SET score = ? WHERE id0 = ? AND id1= ? AND diag = ?
        """ % self.diagonals_table

        def score_calculator(len0, len1, diag, radius, num_seeds):
            band_length = sqrt(
                (len0 - abs(diag)) ** 2 + (len1 - abs(diag)) ** 2
            )
            band_area = 2. * radius * band_length
            mu, sd = binomial_to_normal(band_area, .25 ** self.wordlen)
            z_score = (num_seeds - mu) / float(sd)
            return z_score

        # each entry is a pair of (id, len) tuples
        # FIXME use the sequence table and the attrs for lengths
        seqpairs = list(combinations(self.kmer_index.scanned_sequences(ids), 2))
        with self.connection() as conn:
            cursor = conn.cursor()
            for (id0, len0), (id1, len1) in seqpairs:
                seed_count, diags = self.cum_seed_count(id0, id1, len0, len1)

                for diag, radius in diags:
                    center_idx = diag + len1
                    upper_idx = min(center_idx + radius, len(seed_count) - 1)
                    lower_idx = max(center_idx - radius - 1, 0)

                    num = seed_count[upper_idx] - seed_count[lower_idx]
                    score = score_calculator(len0, len1, diag, radius, num)

                    cursor.execute(query, (score, id0, id1, diag))

    def highest_scoring_band(self, id0s, id1s, min_band_score=None):
        """Returns the diagonal range with the highest score (i.e lowest
        p-value, cf. :func:`score_seeds`) for a given pair of sequences.

        Args:
            id0 (int): The identifier for the first (vertical) sequence.
            id1 (int): The idenfitier for the second (horizontal) sequence.

        Keyword Args:
            min_band_score (int|float): The minimum acceptable score; if the
                highest scoring band has a smaller score it will not be
                reported.  Default is None in which case there is a highest
                scoring band as long as there is at least one seed for
                sequences ``id0`` and ``id1``.

        Returns:
            tuple:
                A 2-tuple containing a diagonal range (as a tuple of inclusive
                upper and lower bounds) with the highest score and its score.
                If no acceptable band is found ``(None, None)`` is returned.
        """  # FIXME docs and tests about multiplicity and order of ids; output
             # is 3 tuple
        if not isinstance(id0s, list):
            id0s = [id0s]
        if not isinstance(id1s, list):
            id1s = [id1s]
        assert id0s and id1s

        assert not set(id0s).intersection(set(id1s))
        # DB records of seeds are all such that id0 < id1, flip them if needed:
        args = tuple()
        idpairs = list(product(id0s, id1s))
        for id0, id1 in idpairs:
            args += (id0, id1) if id0 < id1 else (id1, id0)
        id_condition = ' OR '.join('(id0 = ? AND id1 = ?)' for _ in idpairs)
        query = """
            SELECT id0, id1, diag - radius, diag + radius, score FROM %s
            WHERE (%s)
        """ % (self.diagonals_table, id_condition)
        if min_band_score is not None:
            query += ' AND score >= ?'
            args += (float(min_band_score), )
        query += ' ORDER BY score DESC limit 1'

        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(query, args)
            for id0, id1, min_diag, max_diag, score in cursor:
                if id0 in id0s:
                    return (id0, id1), (min_diag, max_diag), score
                else:
                    return (id1, id0), (-max_diag, -min_diag), score

        return (None, None), (None, None), None

    # FIXME badly named
    def min_num_seeds_in_band(self, id0, id1, len0, len1, bin_len=20, diag_range=None):
        seeds = list(self.seeds(id0, id1, diag_range))
        # FIXME make seeds() return in order of anti; then filter out the first
        # and last chunk of table (and NOT the seeds themselves).

        diag = sum(diag_range) / 2
        # [s, 2*s), [2*s, 3*s),  ... [(k-2)*s, (k-1)*s)
        # anti/bin_len = 1, 2, ..., k - 2
        max_anti = min(len0 - diag, len1) + min(diag, 0)
        k = max_anti / bin_len
        if k <= 2:
            return 0, 0
        counts = [0] * (k - 2)
        for seed in seeds:
            anti = min(seed.pos0, seed.pos1)
            #print 'anti =', anti, 'diag =', seed.pos0 - seed.pos1, seed
            k_ = anti / bin_len
            if k_ >= 1 and k_ < k - 1:
                counts[k_ - 1] += 1
        #print id0, id1, len(seeds), counts
        num_0 = sum(1. if c == 0 else 0 for c in counts)
        return num_0, k - 2
        mu, sd = binomial_to_normal(k - 2, .5)
        return (num_0 - mu) / float(sd)
        #return min(counts)

    def num_seeds_in_band(self, id0, id1, diag_range=None):
        q = """
            SELECT SUM(count) FROM %s
            WHERE id0 = ? AND id1 = ? AND diag BETWEEN ? AND ?
        """ % self.diagonals_table
        with self.connection() as conn:
            return conn.cursor().execute(q, (id0, id1) + diag_range).next()[0]

    def seeds(self, id0, id1, diag_range=None):
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
        assert isinstance(id0, int) and isinstance(id1, int) and id0 < id1

        query = """
            SELECT pos0, pos1 FROM %s
            WHERE id0 = ? AND id1 = ?
        """ % (self.seeds_table)
        args = (id0, id1)
        if diag_range is not None:
            assert isinstance(diag_range, tuple)
            assert len(diag_range) == 2 and diag_range[0] <= diag_range[1]
            query += ' AND pos0 - pos1 BETWEEN ? AND ?'
            args += diag_range

        seed_kw = {'id0': id0, 'id1': id1, 'length': self.wordlen}
        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(query, args)
            for pos0, pos1 in cursor:
                yield Seed(pos0=pos0, pos1=pos1, **seed_kw)

    def plot_seeds(self, ax, id0, id1, **kw):
        """Plots all seeds for a given pair of sequences on a matplitlib
        axes object. All keyword arguments are passed as is to
        ``matplotlib.axes.Axes.scatter``.

        Args:
            ax (matplotlib.axes.Axes): The matplotlib axes object to use.
            id0 (int): The identifier for the first (vertical) sequence.
            id1 (int): The idenfitier for the second (horizontal) sequence.
        """
        _kw = {
            's': 1,
            'color': 'k',
        }
        _kw.update(kw)
        seeds = list(self.seeds(id0, id1))
        ax.scatter([s.pos0 for s in seeds], [s.pos1 for s in seeds], **_kw)
        ax.set_aspect('equal')
        ax.set_xlim(0, None)
        ax.set_ylim(0, None)
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_xlabel('Position in read #%d' % id0)
        ax.set_ylabel('Position in read #%d' % id1)
        ax.xaxis.set_label_position('top')
