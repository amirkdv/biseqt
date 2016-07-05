# -*- coding: utf-8 -*-
"""This module provides tools for k-mer analysis. TODO"""

from math import sqrt, erf, log
import struct

from .database import DB


def binomial_to_normal(n, p):
    """Given the parameters of a binomial distribution, returns the parameters
    of its normal approximation.

    .. math::
        B(n, p) \\simeq \\mathcal{N}(np, \\sqrt{np(1-p)})

    The approximation approaches identity as :math:`n` grows to infinity for
    any fixed :math:`p`.

    Args:
        n (int): First parameter of binomial distribution (i.e number of
            Bernoulli trials).
        p (float): Second parameter of binomial distribution (i.e Bernoulli
            success probability).

    Returns:
        tuple:
            Mean and standard deviation of the approximating normal
            distribution.
    """
    assert p >= 0 and p <= 1 and n > 0
    mu = n * p
    sd = sqrt(n * p * (1 - p))
    return mu, sd


def normal_neg_log_pvalue(mu, sd, x):
    """Gives the negative log p-value of an observation under the null
    hypothesis of normal distribution; that is, given an observation :math:`x`
    from a random variable:

    .. math::
        X \\sim \\mathcal{N}(\\mu, \\sigma)

    this function calculates :math:`-\\log(\\Pr[X\\ge x])` which is:

    .. math::
        -\\log\\left(
            \\frac{1}{2} \\left[1 - \\mathrm{erf} \\left(
                                \\frac{x-\mu}{\\sigma\\sqrt{2}} \\right)
                         \\right]
        \\right)

    Args:
        mu (float): Mean of normal distribution.
        sd (float): Standard deviation of normal distribution.
        x (float): Observation.

    Returns:
        float: A positive real number.
    """
    try:
        return - log(0.5) - log(1 - erf((x - mu) / (sd * sqrt(2))))
    except ValueError:
        # can only happen if the argument to second log is 0
        return float('-inf')


class KmerIndex(object):
    """An index for kmers, their occurences in a body of sequences, and
    :class:`seeds <Seed>`.

    Attributes:
        db (database.DB): The sequence :class:`database <biseqt.database.DB>`.
        wordlen (int): Length of kmers of interest to this index.
    """
    def __init__(self, db, wordlen):
        assert isinstance(db, DB)
        self.db = db
        self.wordlen = wordlen
        db.register('initialize', self.initialize)
        db.register('insert-sequence', self.index_kmers)

        self._digits = '0123456789abcdefghijklmnopqrstuvwxyz'
        assert len(self.db.alphabet) <= len(self._digits), \
            'Maximum alphabet size of %d exceeded' % len(self._digits)
        int_size = 8 * struct.calcsize('P')
        assert self.wordlen < (int_size - 1)/2., \
            'Maximum kmer length %d for %d-bit integers exceeded' % \
            (self.wordlen, int_size)

    _init_script = """
    -- Kmer index initialization script
        CREATE TABLE IF NOT EXISTS kmers_%d (
          'kmer'  INTEGER PRIMARY KEY,  -- The kmer in integer representation.
          'hits'  VARCHAR,              -- The positions where the kmer has
                                        -- been observed; each hit is denoted
                                        -- by "id:pos" and hits are separated
                                        -- by ",".
          'score' REAL DEFAULT NULL     -- The log(p-value) for the number of
                                        -- occurences of this kmer.
        );
        CREATE TABLE IF NOT EXISTS kmer_indexed_%d (
          'id'     INTEGER REFERENCES sequence(id),
                                        -- the id of an indexed sequence
          'length' INTEGER              -- the length of the sequence
        );
    """

    def initialize(self, conn):
        """Event handler for "initialize" (cf. :attr:`DB.events
        <biseqt.database.DB.events>`). Creates two tables:

        * ``kmer_indexed_N`` keeping track of those sequences that have already
          been scanned for kmers.
        * ``kmers_N`` keeping track of observed kmers and their position of
          occurrence in each sequence.

        where ``N`` is the :attr:`word length <wordlen>`.

        Args:
            conn (sqlite3.Connection): An open connection to operate on.
        """
        init_script = self._init_script % (self.wordlen, self.wordlen)
        conn.cursor().executescript(init_script)

    # put the initialization script in the docs
    initialize.__doc__ += '\n\n\t.. code-block:: sql\n\t%s\n' % \
                          '\n\t'.join(_init_script.split('\n'))

    def kmer_as_int(self, contents):
        """Calculates the integer representation of a kmer by treating its
        contents as digits in the base of alphabet length. For instance, in the
        DNA alphabet ``AGA`` becomes :math:`(020)_4` which is 8. Note that
        each kmer gives a unique integer as long as all kmers have the same
        word length (which is the case here). There are two restrictions
        imposed on the word length and alphabet size (enforced in
        :func:`__init__`):

        * The alphabet must be such that all letters can be represented by
          single ASCII characters between ``[0-9a-z]`` (cf. int_). This
          implies a maximum alphabet size of 36.
        * The word length must be such that a single integer can store the
          entire representation of a kmer. This requires that we have:

          .. math::
            k < \\frac{I-1}{2}

          where :math:`k` is the word length and :math:`I` is the number of
          bits allocated for an integer. For instance, on a 64-bit system the
          maximum word length is 31.

        .. _int: https://docs.python.org/2/library/functions.html#int

        .. wikisection:: dev
            :title: Integer Sizes

            All kmers are stored as their integer representation to save on
            space and processing time. Python_ is flexible with the maximum
            size of integers, as integers automatically switch to longs, which
            have "unlimited precision". SQLite_, too, is flexible but has a
            maximum integer cap of 64-bits: integers take 2, 4, 6, or 8 bytes
            per integer dependning on the size it needs.

            The checks resulting from maximum integer size are performed in
            :func:`KmerIndex.__init__ <KmerIndex>` which basically block kmers
            taking more than 64-bit integers to represent.

            .. _python: https://docs.python.org/2/library/stdtypes.html\
                        #numeric-types-int-float-long-complex
            .. _sqlite: https://www.sqlite.org/datatype3.html#section_2
        """
        # document dependence on word length
        as_str = ''.join(self._digits[c] for c in contents)
        return int(as_str, len(self.db.alphabet))

    def scan_kmers(self, seq):
        """A generator for kmer hit tuples of the form ``(kmer, pos)``.

        Args:
            seq (sequence.Sequence): The sequence to be scanned.

        Yields:
            tuple:
                The first element is the kmer in integer form (cf.
                :func:`kmer_as_int`) and the second element is the starting
                position of the kmer in ``seq``.
        """

        for idx in range(len(seq) - self.wordlen + 1):
            kmer = self.kmer_as_int(seq.contents[idx: idx + self.wordlen])
            yield (kmer, idx)

    def index_kmers(self, conn, seq, rec):
        """Event handler for "insert-sequence" (cf. :attr:`database.events
        <biseqt.database.events>`). Indexes all kmers observed in the given
        sequence in ``kmers_N``. For example, for word length 3 and a DNA
        sequence ``ACACAC`` with identifier 1 we get::

            kmer| hits          | score
            4   | 1:0,1:2,1:4   | NULL
            34  | 1:1,1:3,1:4   | NULL

        Args:
            conn (sqlite3.Connection): SQLite connection to operate on.
            seq (sequence.Sequence): The sequence just inserted into the
                database.
            rec (database.Record): The record object corresponding to the
                insertion of ``seq`` with the :attr:`id
                <biseqt.database.Record.id>` field populated.

        .. wikisection:: dev
            :title: Insert-or-Append Queries

            In many places we need to perform an SQL query which either
            inserts a new record or appends the given value to some column
            in case of conflict. This is achived by using one of SQLite's
            conflict resolution mechanisms_ ``ON CONFLICT REPLACE`` or in short
            ``OR REPLACE``. For instance, consider a table::

                id | field
                1  | foo

            where we wish ``INSERT INTO ... (id, field) VALUES (2, 'bar')``
            to give::

                id | field
                1  | foo
                2  | bar

            and ``INSERT INTO ... (id, field) VALUES (1, 'bar')`` to give::

                id | field
                1  | foo,bar

            This can be implemented by using the following query format:

            .. code-block:: sql

                INSERT INTO ... (id, field) VALUES
                SELECT ?, IFNULL(SELECT field FROM ... WHERE id = ?, "") || ?

            invoked like this:

            .. code-block:: python

                id, field = ...
                cursor = sqlite3.connect(...).cursor()
                cursor.execute(query, (id, id, ',' + field))

            Note that this pattern only works if the ``id`` column has a unique
            constraint on it. Otherwise, no conflict will arise to be resolved
            and new values will appear in new records instead of being
            appended to old ones.

            .. _mechanisms: https://www.sqlite.org/lang_conflict.html
        """
        # this only works if there is a unique constraint on the kmer column.
        hit_query = """
            INSERT OR REPLACE INTO kmers_%d (kmer, hits)
                SELECT ?, IFNULL(
                            (SELECT hits FROM kmers_%d WHERE kmer = ?), ""
                          ) || ?
        """ % (self.wordlen, self.wordlen)
        records = ((kmer, kmer, '%d:%d,' % (rec.id, pos))
                   for kmer, pos in self.scan_kmers(seq))
        conn.cursor().executemany(hit_query, records)

        status_query = """
            INSERT INTO kmer_indexed_%d (id, length) VALUES (?, ?)
        """ % self.wordlen
        conn.cursor().execute(status_query, (rec.id, len(seq)))

    def total_length_indexed(self):
        """The total number of letters, among all sequences, indexed so far for
        kmers.

        Returns:
            int
        """
        with self.db.connect() as conn:
            cursor = conn.cursor()
            q = 'SELECT SUM(length) FROM kmer_indexed_%d' % self.wordlen
            cursor.execute(q)
            return int(cursor.next()[0])

    def num_kmers(self):
        """The total number of kmers observed so far.

        Returns:
            int
        """
        with self.db.connect() as conn:
            cursor = conn.cursor()
            cursor.execute('SELECT COUNT(*) FROM kmers_%d' % self.wordlen)
            return int(cursor.next()[0])

    def score_kmers(self, only_missing=True):
        """Calculates the negative log p-value for the number of occurences of
        each kmer under the null hypothesis of a binomial distribution. The
        binomial distribution is approximated by a normal distribution for
        numeric stability (cf. :func:`binomial_to_normal`) and the a Bonferroni
        correction for the total number of kmers (cf. :func:`num_kmers`) is
        applied to the raw negative log p-value given by
        :func:`normal_neg_log_pvalue`. The higher the score of a kmer, the more
        likely it is that it belongs to a repeat structure.

        Keyword Args:
            only_missing (bool): Whether to re-score all kmers or only those
                with a ``NULL`` score. Default is True in which case kmers
                with a score are not re-evaluated.
        """
        N = self.num_kmers()
        L = self.total_length_indexed()
        kmer_probability = 1./(len(self.db.alphabet) ** self.wordlen)
        mu, sd = binomial_to_normal(L, kmer_probability)

        def score_calculator(num_occurrences):
            # bonferroni correction
            return - log(N) + normal_neg_log_pvalue(mu, sd, num_occurrences)

        select = 'SELECT kmer, hits FROM kmers_%d ' % self.wordlen + \
                 ('WHERE score IS NULL' if only_missing else '')

        update = 'UPDATE kmers_%d SET score = ? WHERE kmer = ?' % self.wordlen
        with self.db.connect() as conn:
            select_cursor, insert_cursor = conn.cursor(), conn.cursor()
            select_cursor.execute(select)
            for kmer, hits in select_cursor:
                score = score_calculator(hits.count(':'))
                insert_cursor.execute(update, (score, kmer))

    def scanned_sequences(self):
        """Yields the ids and lengths of all scanned sequences from the
        ``kmer_indexed_N`` table.

        Yields:
            tuple:
                The sequence integer identifier and the length of the sequence.
        """
        query = 'SELECT id, length FROM kmer_indexed_%d' % self.wordlen
        with self.db.connect() as conn:
            cursor = conn.cursor()
            cursor.execute(query)
            for _id, _len in cursor:
                yield _id, _len

    def kmers(self, max_score=None):
        """Lazy-loads the observed kmers, their occurences, and their score.

        Keyword Args:
            max_score (float): The maximum score (which is negative log
                p-value) for kmers to be included (cf.  :func:`score_kmers`).
                Default is None in which case all kmers are included.

        Yields:
            tuple:
                The kmer in integer representation, a list of occurences where
                each occurence is a tuple of sequence id and position, and
                the score for the kmer.
        """
        query = 'SELECT kmer, hits, score from kmers_%s' % self.wordlen
        if max_score is not None:
            query += ' WHERE score < %f' % max_score

        with self.db.connect() as conn:
            cursor = conn.cursor()
            cursor.execute(query)
            for kmer, hits, score in cursor:
                hits = [tuple(int(i) for i in hit.split(':'))
                        for hit in hits.split(',') if hit]
                yield kmer, hits, score
