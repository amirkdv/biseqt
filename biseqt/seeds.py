# -*- coding: utf-8 -*-
"""
.. wikisection:: overview
    :title: Statistical Seed Analysis

    The :mod:`biseqt.seeds` module provides statistical tools for analyzing
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

from scipy.special import erfinv
from collections import namedtuple
from itertools import combinations
from math import sqrt

from .kmers import KmerIndex, binomial_to_normal


# TODO implement band_radius_calculator to avoid calling erfinv too many times.
# TODO move this to random and rename to stochastics
def band_radius(len0, len1, diag, gap_prob=None, sensitivity=None):
    """Calculates the smallest band radius such an overlap alignment, with the
    given gap probability, stays entirely within the diagonal band centered at
    the given diagonal. This is given by:

    .. math::
        r = 2\\sqrt{g(1-g)K}
            \\mathrm{erf}^{-1}\\left(1-\\frac{2\\epsilon}{3}\\right)

    where :math:`g` is the gap probability, :math:`1-\\epsilon` is the desired
    sensitivity, and :math:`K` is the "expected" length of the alignment given
    by:

    .. math::
        K = \\left(\\frac{2}{2 - g}\\right) L

    where :math:`L` is the maximum possible length of the alignment:

    .. math::
        L = \\min(l_0 - d, l_1) + \\min(d, 0)

    with :math:`l_0,l_1` being the length of the sequences (i.e ``len0`` and
    ``len1`` arguments) and :math:`d` the starting diagonal (i.e ``diag``
    argument).

    Diagonals are numbered as follows in the dynamic programming table::

        0 -1 -2 -3 ... -len1
        1
        2
        3
        .
        .
        .
        +len0


    Args:
        len0 (int): Length of the first sequence (the "vertical" sequence in
            the table).
        len1 (int): Length of the second sequence (the "horizontal" sequence in
            the table).
        diag (int): Starting diagonal of alignments to consider.
        gap_prob (float): Probability of indels occuring at any position of an
            alignment.
        sensitivity (float): The probability that an alignment with given gap
            probability remains entirely within the band.
    Returns:
        int: The smallest band radius guaranteeing the required sensitivity.

    """
    assert sensitivity > 0 and sensitivity < 1
    assert gap_prob > 0 and gap_prob < 1

    adjusted_sensitivity = 1 - 2 * (1. - sensitivity) / 3

    max_alignment_length = min(len0 - diag, len1) + min(diag, 0)
    expected_alignment_length = (2 / (2. - gap_prob)) * max_alignment_length
    assert expected_alignment_length >= 0
    radius = 2 * erfinv(adjusted_sensitivity) * sqrt(
        gap_prob * (1 - gap_prob) * expected_alignment_length
    )
    return max(1, int(radius))


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
        self.diagonals_table = 'diagonals_%d' % self.wordlen

        def drop_all_indices(*args, **kwargs):
            self.drop_sql_index(self.seeds_table)
            self.drop_sql_index(self.diagonals_table)

        self.db.add_event_listener('db-initialized', self.initialize)
        self.db.add_event_listener('sequences-loading', drop_all_indices)

    _init_script = """
        CREATE TABLE IF NOT EXISTS %s ( -- seeds table
          'id0'   INTEGER,              -- REFERENCES sequence(id)
                                        -- but do not declare it to save time
                                        -- on checking referential integrity.
          'id1'   INTEGER,              -- REFERENCES kmers_N_hits(id), same.
          'pos0'  INTEGER,              -- starting position in sequence 0.
          'pos1'  INTEGER               -- starting position in sequence 1.
        );
        CREATE TABLE IF NOT EXISTS %s ( -- diagonals table
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
            self._init_script % (self.seeds_table, self.diagonals_table)
        )

    # put the initialization script in the docs
    initialize.__doc__ += '\n\n\t.. code-block:: sql\n\t%s\n' % \
                          '\n\t'.join(_init_script.split('\n'))

    def create_sql_index(self, table=None):
        """Creates a ``(id0, id1)`` SQL index over a given table.

        Args:
            table (str): Either of :attr:`seeds_table` or
                :attr:`diagonals_table`.
        """
        assert table in [self.seeds_table, self.diagonals_table]
        self.log('Creating SQL index for %s.' % table)
        with self.db.connection() as conn:
            conn.cursor().execute("""
                CREATE INDEX IF NOT EXISTS %s_seqpair ON %s (id0, id1)
            """ % (table, table))

    def drop_sql_index(self, table=None):
        """Drops SQL indices created by :func:`create_sql_index`.

        Args:
            table (str): Either of :attr:`seeds_table` or
                :attr:`diagonals_table`.
        """
        assert table in [self.seeds_table, self.diagonals_table]
        self.log('Dropping SQL index for %s.' % table)
        with self.db.connection() as conn:
            conn.cursor().execute("""
                DROP INDEX IF EXISTS %s_seqpair;
            """ % (table))

    def log(self, *args, **kwargs):
        """Wraps :func:`log <biseqt.database.DB.log>` of :attr:`db`."""
        self.db.log(*args, **kwargs)

    def index_seeds(self, max_kmer_score=None):
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
        def _records():
            kmers = self.kmer_index.kmers(max_score=max_kmer_score)
            for kmer, hits, _ in kmers:
                for (id0, pos0), (id1, pos1) in combinations(hits, 2):
                    if id0 == id1:
                        continue
                    yield id0, id1, pos0, pos1

        self.log('Indexing seeds (max_kmer_score=%s).' % max_kmer_score)
        with self.db.connection() as conn:
            conn.cursor().executemany(
                """INSERT INTO %s (id0, id1, pos0, pos1) VALUES (?, ?, ?, ?)
                """ % self.seeds_table,
                _records()
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
        with self.db.connection() as conn:
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

    # FIXME the sequence is:
    # 1. index_seeds()
    # 2. count_seeds_on_diagonals()
    # 3. calculate_band_radii()
    # 4. score_diagonals()
    # simplify it!
    def count_seeds_on_diagonals(self):
        """Populates the :attr:`diagonals table <diagonals_table>` form the
        contents of the :attr:`seeds table <seeds_table>`. Radius and score
        information are left blank."""
        self.create_sql_index(self.seeds_table)
        self.log('Counting seeds on each diagonal of each sequence pair.')
        with self.db.connection() as conn:
            conn.cursor().execute("""
                INSERT INTO %s (id0, id1, diag, count)
                SELECT id0, id1, pos0 - pos1, COUNT(*) FROM %s
                GROUP BY id0, id1, pos0 - pos1
            """ % (self.diagonals_table, self.seeds_table))

    def calculate_band_radii(self, gap_prob=None, sensitivity=None):
        """Calculates the diagonal band radius for each entry in the
        :attr:`diagonals table <diagonals_table>` using :func:`band_radius`.

        Keyword Args:
            gap_prob (float): Passed as is to :func:`band_radius`.
            sensitivity (float): Passed as is to :func:`band_radius`.
        """
        assert sensitivity > 0 and sensitivity < 1
        assert gap_prob > 0 and gap_prob < 1

        band_radius_kw = {'gap_prob': gap_prob, 'sensitivity': sensitivity}

        select = """
            SELECT diag FROM %s WHERE id0 = ? AND id1 = ?
        """ % self.diagonals_table

        update = """
            UPDATE %s SET radius = ? WHERE id0 = ? AND id1 = ? AND diag = ?
        """ % self.diagonals_table

        self.log('Calculating band radii (gap_prob=%s, sensitivity=%s)' %
                 (str(gap_prob), str(sensitivity)))
        # each entry is a pair of (id, len) tuples
        seqpairs = combinations(self.kmer_index.scanned_sequences(), 2)
        self.create_sql_index(self.diagonals_table)
        with self.db.connection() as conn:
            cursor = conn.cursor()
            for (id0, len0), (id1, len1) in seqpairs:
                cursor = conn.cursor()
                diags = [x[0] for x in cursor.execute(select, (id0, id1))]
                recs = ((band_radius(len0, len1, diag, **band_radius_kw),
                         id0, id1, diag) for diag in diags)
                conn.cursor().executemany(update, recs)

    def score_diagonals(self):
        """Scores all diagonal bands containing seeds for each pair of sequences.
        The radius for each diagonal band is assumed to be calculated by
        :func:`calculate_band_radii`. The score for a band is the negative log
        p-value of the observed number of seeds in it under the null hypothesis
        that seeds occur randomly throughout the dynamic programming table.

        Explicitly, for a fixed pair of sequences with lengths :math:`l_0, l_1`
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
        (cf. :func:`biseqt.kmers.binomial_to_normal`) and the z-score of the
        observed number of seeds is reported, namely:

        .. math::

                \\mathrm{score} = \\frac{n-Ap}{\\sqrt{Ap(1-p)}}

        A high scoring diagonal band is less likely to be accidental and more
        likely to indicate an overlap between the two sequences.
        """
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

        self.log('Scoring all diagonals for each pair of sequences.')
        # each entry is a pair of (id, len) tuples
        seqpairs = combinations(self.kmer_index.scanned_sequences(), 2)
        with self.db.connection() as conn:
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

    def highest_scoring_band(self, id0, id1, min_band_score=None):
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
        """
        assert isinstance(id0, int) and isinstance(id1, int) and id0 < id1
        query = """
            SELECT diag - radius, diag + radius, score FROM %s
            WHERE id0 = ? AND id1 = ?
        """ % self.diagonals_table
        args = (id0, id1)
        if min_band_score is not None:
            query += ' AND score >= ?'
            args += (float(min_band_score), )
        query += ' ORDER BY score DESC limit 1'

        with self.db.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(query, args)
            for min_diag, max_diag, score in cursor:
                return (min_diag, max_diag), score

        return (None, None)

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
        with self.db.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(query, args)
            for pos0, pos1 in cursor:
                yield Seed(pos0=pos0, pos1=pos1, **seed_kw)
