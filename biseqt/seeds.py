# -*- coding: utf-8 -*-
"""This module provides statistical tools for analyzing matching segment pairs
(aka seeds) between large numbers of sequences. TODO"""

from collections import namedtuple
from itertools import combinations
from scipy.special import erfinv
from math import sqrt, log

from .kmers import KmerIndex


def band_radius(len0, len1, diag, gap_prob=None, sensitivity=None):
    """Calculates the smallest band radius such an overlap alignment, with the
    given gap probability, stays entirely within the diagonal band centered at
    the given diagonal. This is given by:

    .. math::
        r = 2\\mathrm{erf}^{-1}(\\rho)\\sqrt{g(1-g)K}

    where :math:`g` is the gap probability, :math:`\\rho` is the desired
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

    # FIXME document
    sensitivity = 1./3 + sensitivity * 2./3

    max_alignment_length = min(len0 - diag, len1) + min(diag, 0)
    expected_alignment_length = (2. / (2 - gap_prob)) * max_alignment_length
    assert expected_alignment_length >= 0
    radius = 2 * erfinv(sensitivity) * sqrt(
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
    """
    def __init__(self, kmer_index):
        assert isinstance(kmer_index, KmerIndex)
        self.kmer_index = kmer_index
        self.wordlen = self.kmer_index.wordlen
        self.db = self.kmer_index.db
        self.db.register('initialize', self.initialize)

    _init_script = """
        CREATE TABLE IF NOT EXISTS seeds_%d (
          'id0'   INTEGER,              -- REFERENCES sequence(id)
                                        -- but do not declare it to save time
                                        -- on checking referential integrity.
          'id1'   INTEGER,              -- REFERENCES sequence(id), same.
          'diag'  INTEGER,              -- a single diagonal in the edit graph.
          'radius'INTEGER DEFAULT NULL, -- the diagonal band radius for the two
                                        -- sequences with diag as center.
          'score' REAL DEFAULT NULL,    -- the log(p-value) for the number of
                                        -- seeds on the band centered at this
                                        -- diagonal.
          'seeds' VARCHAR,              -- all seeds on this diagonal separated
                                        -- by ",". Each seed is denoted by
                                        -- "pos0:pos1" containing the starting
                                        -- position on id0 and id1.
          UNIQUE (id0, id1, diag)       -- needed for insert-or-append queries.
        );
    """

    def initialize(self, conn):
        """Event handler for "initialize" (cf. :attr:`DB.events
        <biseqt.database.DB.events>`). Creates a single table:

        * ``seeds_N`` keeping track of all observed :class:`seeds <Seed>`
          between pairs of sequences.

        Args:
            conn (sqlite3.Connection): An open connection to operate on.
        """
        init_script = self._init_script % self.wordlen
        conn.cursor().executescript(init_script)

    # put the initialization script in the docs
    initialize.__doc__ += '\n\n\t.. code-block:: sql\n\t%s\n' % \
                          '\n\t'.join(_init_script.split('\n'))

    def index_seeds(self, max_kmer_score=None):
        """Indexes all seeds and their diagonal positions. For each kmer with
        :math:`n` hits from *distinct* sequences :math:`n\\choose2` seeds
        are created. If a kmer occurs multiple times in a sequence seeds
        between different positions of the same sequence are not considered.

        Keyword Args:
            max_kmer_score: The maximum score beyond which kmers are not
                considered for seeds. High scoring words are more likely to
                belong to repetitive regions (cf. :func:`score_kmers
                <biseqt.kmers.KmerIndex.score_kmers>`). Default is None in
                which case all kmers are considered.
        """
        # this only works if there is a unique constraint on (id0, id1, diag)
        query = """
            INSERT OR REPLACE INTO seeds_%d (id0, id1, diag, seeds)
                SELECT ?, ?, ?,
                    IFNULL(
                        (SELECT seeds FROM seeds_%d
                         WHERE id0 = ? AND id1 = ? AND diag = ?), ""
                    ) || ?
        """ % (self.wordlen, self.wordlen)

        def _records():
            kmers = self.kmer_index.kmers(max_score=max_kmer_score)
            for kmer, hits, _ in kmers:
                for (id0, pos0), (id1, pos1) in combinations(hits, 2):
                    if id0 == id1:
                        continue
                    diag = pos0 - pos1
                    unique = (id0, id1, diag)
                    yield unique + unique + ('%d:%d,' % (pos0, pos1),)

        self.kmer_index.score_kmers(only_missing=True)
        with self.db.connect() as conn:
            conn.cursor().executemany(query, _records())

    def cum_seed_count(self, id0, id1, len0, len1):
        """Returns the cumulative count of seeds per diagonals.

        Args:
            id0 (int): The identifier for the first (vertical) sequence.
            id1 (int): The idenfitier for the second (horizontal) sequence.
            len0 (int): The length of the first sequence.
            len1 (int): The length of the second sequence.

        Returns:
            tuple:
                A list of length ``len0 + len1 + 1`` containing the cumulative
                count of seeds upto each diagonal from ``-len1`` to ``+len0``,
                and a list of all diagonals containing seeds.
        """
        query = """
            SELECT diag, seeds FROM seeds_%d
            WHERE id0 = ? AND id1 = ? ORDER BY diag ASC
        """ % self.wordlen

        diags = []
        seed_count = [0] * (len1 + len0 + 1)
        prev_diag = 0
        with self.db.connect() as conn:
            cursor = conn.cursor()
            cursor.execute(query, (id0, id1))
            # count seeds in diagonal band
            for diag, seeds in cursor:
                diags.append(diag)

                diag += len1
                for _diag in range(prev_diag, diag + 1):
                    seed_count[_diag] = seed_count[prev_diag]
                seed_count[diag] += seeds.count(':')
                prev_diag = diag

        for _diag in range(prev_diag, len(seed_count)):
            seed_count[_diag] = seed_count[prev_diag]
        return seed_count, diags

    def score_seeds(self, only_missing=True, gap_prob=None, sensitivity=None):
        """Scores all diagonal bands containing seeds for each pair of sequences.
        The radius for each diagonal band is calculated using
        :func:`band_radius` and the score for a band is the negative log
        p-value of the observed number of seeds in it under the null hypothesis
        that seeds occur randomly throughout the dynamic programming table.

        Explicitly, for a fixed pair of sequences with lengths :math:`l_0, l_1`
        a diagonal band centered at :math:`d` and with radius :math:`r`, the
        dimensions of the dynamic programming table is :math:`[0, l_0] \\times
        [0, l_1]` and the probability of observing :math:`n` seeds in the
        diagonal band :math:`[d-r, d+r]` is:

        .. math::
            \\Pr(X\\ge n) \\simeq \\left(\\frac{A}{l_0l_1}\\right)^n

        where :math:`A` is the area of the diagonal band:

        .. math::

            A \\simeq 2r\\sqrt{(l_0-|d|)^2 + (l_1-|d|)^2}

        Since the same hypothesis is tested against all diagonal bands (there
        are roughly :math:`l_0+l_1` total hypotheses) a Bonferroni correction
        term is added to give the final formula for the corrected negative log
        p-value:

        .. math::

                - \\log(l_0+l_1) -
                n \\left[\\log(A) - \\log(l_0) - \\log(l_1)\\right]

        A high scoring diagonal band is less likely to be accidental and more
        likely to indicate an overlap between the two sequences.

        Keyword Args:
            gap_prob (float): Passed as is to :func:`band_radius`.
            sensitivity (float): Passed as is to :func:`band_radius`.
        """
        assert sensitivity > 0 and sensitivity < 1
        assert gap_prob > 0 and gap_prob < 1

        query = """
            UPDATE seeds_%d SET score = ?, radius = ?
            WHERE id0 = ? AND id1 = ? AND diag = ?
        """ % self.wordlen

        def score_calculator(len0, len1, diag, radius, num_seeds):
            band_width = 2 * radius
            band_length = sqrt(
                (len0 - abs(diag)) ** 2 + (len1 - abs(diag)) ** 2
            )
            band_area = band_width * band_length
            seed_contribution = log(band_area) - log(len0) - log(len1)
            # TODO consider caching seed contributions for multiples of 100
            # as it barely changes.
            # TODO consider dropping the bonferroni term as it depends only on
            # len0 and len1, and not on what differentiates bands.
            return - log(len0 + len1) - num_seeds * seed_contribution

        band_radius_kw = {'gap_prob': gap_prob, 'sensitivity': sensitivity}
        with self.db.connect() as conn:
            cursor = conn.cursor()
            scanned_sequences = self.kmer_index.scanned_sequences()
            for (id0, len0), (id1, len1) in combinations(scanned_sequences, 2):
                seed_count, diags = self.cum_seed_count(id0, id1, len0, len1)
                # update scores
                for diag in diags:
                    radius = band_radius(len0, len1, diag, **band_radius_kw)
                    center_idx = diag + len1
                    upper_idx = min(center_idx + radius, len(seed_count) - 1)
                    lower_idx = max(center_idx - radius - 1, 0)
                    num_seeds = seed_count[upper_idx] - seed_count[lower_idx]
                    score = score_calculator(len0, len1, diag, radius,
                                             num_seeds)
                    args = (score, radius, id0, id1, diag)
                    cursor.execute(query, args)

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
                Inclusive upper and lower bounds of the highest scoring
                diagonal band, or None if no acceptable band is found.
        """
        assert isinstance(id0, int) and isinstance(id1, int) and id0 < id1
        query = """
            SELECT diag - radius, diag + radius FROM seeds_%d
            WHERE id0 = ? AND id1 = ?
        """ % self.wordlen
        args = (id0, id1)

        if min_band_score is not None:
            query += ' AND score >= ?'
            args += (float(min_band_score), )
        query += ' ORDER BY score DESC limit 1'

        with self.db.connect() as conn:
            cursor = conn.cursor()
            cursor.execute(query, args)
            for min_diag, max_diag in cursor:
                return min_diag, max_diag

            return None

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
            SELECT seeds, score FROM seeds_%d
            WHERE id0 = ? AND id1 = ?
        """ % (self.wordlen)
        args = (id0, id1)
        if diag_range is not None:
            assert isinstance(diag_range, tuple)
            assert len(diag_range) == 2 and diag_range[0] <= diag_range[1]
            query += ' AND diag BETWEEN ? AND ?'
            args += diag_range

        seed_kw = {'id0': id0, 'id1': id1, 'length': self.wordlen}
        with self.db.connect() as conn:
            cursor = conn.cursor()
            cursor.execute(query, args)
            for seeds, score in cursor:
                seeds = tuple(tuple(int(i) for i in seed.split(':'))
                              for seed in seeds.split(',') if seed)
                for pos0, pos1 in seeds:
                    yield Seed(pos0=pos0, pos1=pos1, **seed_kw), score
