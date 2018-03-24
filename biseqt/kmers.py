# -*- coding: utf-8 -*-
# FIXME docs
"""
.. wikisection:: overview
    :title: (4) Handling Kmers

    The :mod:`biseqt.kmers` module provides tools for k-mer analysis. Kmers are
    represented as integer (and hence a maximum word length is imposed by
    integer size of the machine).

    >>> from biseqt.sequence import Alphabet, Sequence
    >>> from biseqt.kmers import as_kmer_seq, KmerIndex
    >>> A = Alphabet('ACGT')
    >>> S = A.parse('AAACGCGT')
    >>> wordlen = 3
    >>> as_kmer_seq(S, wordlen)
    [0, 1, 6, 25, 38, 27]
    >>> kmer_index = KmerIndex(path=':memory:', alphabet=A, wordlen=wordlen)
    >>> kmer_index.index_kmers(S)
    1 # the internal 'id' assigned to sequence S
    >>> kmer_index.kmers()
    [0, 1, 6, 25, 27, 38]
    >>> kmer_index.hits(27)
    [(1, 5)] # (seqid, position)

.. wikisection:: dev
    :title: Insert-or-Append Queries

    It may be useful to perform an SQL query which either inserts a new record
    or appends the given value to some column in case of conflict. This is
    achived by using one of SQLite's conflict resolution mechanisms_ ``ON
    CONFLICT REPLACE`` or in short ``OR REPLACE``. For instance, consider a
    table::

        id | field
        1  | foo

    where we wish ``INSERT INTO ... (id, field) VALUES (2, 'bar')`` to give::

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
        conn = apsw.Connection('example.db')
        conn.cursor().execute(query, (id, id, ',' + field))

    Note that this pattern only works if the ``id`` column has a unique
    constraint on it. Otherwise, no conflict will arise to be resolved and new
    values will appear in new records instead of being appended to old ones.

    .. _mechanisms: https://www.sqlite.org/lang_conflict.html


.. wikisection:: dev
    :title: SQLite Performance Tuning

    Tuning strategy naturally depends on the balance between the volume of
    insert vs. select queries and the concurrency requirements. Here we will
    assume:

    * The volume of inserts is much larger than selects,
    * Application logic can be trusted with respecting unique constraints (i.e
      the code creating data does not violate semantic constraints).
    * Usage pattern consists of bulk of inserts followed by bulk of selects
      (e.g. not interleaved).

    Under these circumstances, the following guidelines are suggested:

    * Create indices after bulk inserts not before. For instance, instead of:

      .. code-block:: sql

        CREATE TABLE foo ('f1' int, 'f2' int, UNIQUE(f1, f2))
        -- INSERT INTO foo VALUES ...

      it's more performant to say:

      .. code-block:: sql

        CREATE TABLE foo ('f1' int, 'f2' int)
        -- INSERT INTO foo VALUES ...
        CREATE UNIQUE INDEX foo_index ON foo (f1, f2)

    * If there is no concern about data corruption upon application or
      operating sytem crash journaling can be turned off as well, from `docs
      <journaling_docs>`_:

          | If the application crashes in the middle of a transaction when the
          | OFF journaling mode is set, then the database file will very likely
          | go corrupt.

      To turn off journaling:

      .. code-block:: sql

        PRAGMA journaling_mode = OFF

      Note that turning off journaling breaks rollbacks:

        | The OFF journaling mode disables the rollback journal completely. No
        | rollback journal is ever created and hence there is never a rollback
        | journal to delete. The OFF journaling mode disables the atomic commit
        | and rollback capabilities of SQLite.
    * When a table has a unique integer key it should be declared as ``INTEGER
      PRIMARY KEY`` so that it would take over the default ``rowid`` field.
      This saves space (and thus a small amount of time) on both the field and
      the corresponding index.
    * Foreign key constraints, if enforced, slow down bulk inserts
      significantly. However, by default, `foreign key checks`_ are turned off.
      To turn it on:

      .. code-block:: sql

        PRAGMA foreign_keys = ON;

      Note that this default maybe modified by compile time flags (i.e foreign
      keys may be turned on by default). Furthermore, if foreign keys are
      turned on, consider deferring_ foreign key enforcement to transaction
      commits and and keep in mind that ``pysqlite`` (following python's
      DB-API) fudges with transaction beginning and ends.
    * Larger `page sizes <pagesize_docs>`_ can marginally improve read/write
      performance. To increase the page size:

      .. code-block:: sql

        PRAGMA page_size = 65536

    .. _journaling_docs: https://www.sqlite.org/pragma.html#pragma_journal_mode
    .. _pagesize_docs: https://www.sqlite.org/pragma.html#pragma_page_size
    .. _foreign key checks: https://www.sqlite.org/foreignkeys.html#fk_enable
    .. _deferring: https://www.sqlite.org/foreignkeys.html#fk_deferred
    .. rubric: References

        * http://stackoverflow.com/a/1712873
        * http://codereview.stackexchange.com/q/26822
        * http://stackoverflow.com/q/3134900
"""

import struct
import os
import apsw
import logging

from .util import Logger
from .sequence import Alphabet, Sequence

DIGITS = '0123456789abcdefghijklmnopqrstuvwxyz'


def kmer_as_int(contents, alphabet):
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
    assert isinstance(alphabet, Alphabet)
    # TODO document dependence on word length
    as_str = ''.join(DIGITS[c] for c in contents)
    return int(as_str, len(alphabet))


def as_kmer_seq(seq, wordlen, mask=[]):
    """A generator for kmer hit tuples of the form ``(kmer, pos)``. Kmers
    are represented in integer form (cf. :func:`kmer_as_int`).

    Args:
        seq (sequence.Sequence): The sequence to be scanned.
        wordlen (int): Size of kmers.
        mask (list): A list of sets of integers ``(i_1, ..., i_k)`` which mask
            kmers (represented by ``None``) if the kmer content (set of letters
            appearing in the kmer, represented as integers as in
            :attr:`Sequence.contents`) matches the set.

    Returns:
        list: of integers representing kmers.
    """
    assert isinstance(seq, Sequence)
    assert all(isinstance(lets, set) for lets in mask)
    kmers = []
    for pos in range(len(seq) - wordlen + 1):
        if mask:
            lets = set(seq[pos: pos + wordlen])
            if lets in mask:
                kmers.append(None)
                continue
        kmer = kmer_as_int(seq.contents[pos: pos + wordlen], seq.alphabet)
        kmers.append(kmer)
    return kmers


class KmerDBWrapper(object):
    """Generic wrapper for an SQLite database for Kmers.

    Attributes:
        name (str): String name used as a suffix for table names.
        path (str): Path to the SQLite datbase (or ``:memory:``).
        alphabet (sequence.Alphabet):
            The alphabet for sequences in the database.
        wordlen (int): Length of kmers of interest to this index.
        mask (list): A list of sets of integers which mask kmers (represented
            by ``None``), cf. :func:`as_kmer_seq`.
        init_script (str): SQL script to be executed upon initialization;
            typically creates tables needed by the class.
    """
    def __init__(self, name='', path=':memory:', alphabet=None, wordlen=None,
                 mask=[], log_level=logging.INFO, init_script=None):
        self.name = name
        assert all(isinstance(lets, set) for lets in mask)
        self.mask = mask
        assert isinstance(wordlen, int)
        assert isinstance(alphabet, Alphabet)
        assert len(alphabet) <= len(DIGITS), \
            'Maximum alphabet size of %d exceeded' % len(DIGITS)
        self.alphabet = alphabet
        int_size = 8 * struct.calcsize('P')
        assert wordlen < (int_size - 1)/2., \
            'Maximum kmer length %d for %d-bit integers exceeded' % \
            (wordlen, int_size)
        self.wordlen = wordlen

        if path == ':memory:':
            self.path = path
            relpath = path
        else:
            self.path = os.path.abspath(path)
            assert os.path.exists(self.path) or \
                os.access(os.path.dirname(self.path), os.W_OK), \
                'Database %s is not writable' % self.path
            relpath = os.path.relpath(self.path, os.getcwd())

        log_header = '%d-mer cache (%s)' % (self.wordlen, relpath)
        self._logger = Logger(log_level=log_level, header=log_header)
        self.log_level = log_level
        self._connection = None

        self.init_script = init_script
        if self.init_script:
            with self.connection() as conn:
                conn.cursor().execute(self.init_script)

    # FIXME compare with old-tip and figure out what the deal with resetting is
    def connection(self, reset=False):
        """Provides a SQLite database connection that can be used as a context
        manager. The returned object is always the same connection object
        belonging to the :class:`KmerIndex` instance (otherwise in-memory
        connections would reset the database contents upon every invocation).

        Returns:
            apsw.Connection
        """
        if reset or self._connection is None:
            self._connection = apsw.Connection(self.path)
        return self._connection

    def log(self, *args, **kwargs):
        """Wraps :class:`Logger.log`."""
        self._logger.log(*args, **kwargs)


class KmerCache(KmerDBWrapper):
    """A cache backed by SQLite for representations of sequences as integer
    sequences. Upon initialization the following SQL script is executed

    .. code-block:: sql

        CREATE TABLE seq_kmers (
            'seq' VARCHAR,   -- content identifier of sequence
            'kmres' VARCHAR, -- comma separated representation of sequence
                             -- as integers encoding kmers.
        );

    This implies that any time only one ``KmerCache`` can exist with the same
    path.
    FIXME: using the same database for different word lengths / alphabets is
    quietly accepted (and wrong results returned).
    """
    def __init__(self, name='', **kw):
        self.name = name
        init_script = """
            CREATE TABLE IF NOT EXISTS %s (
                'seq' VARCHAR,   -- content identifier of sequence
                'kmers' VARCHAR  -- comma separated representation of sequence
                                 -- as integers encoding kmers.
            );
        """ % self.kmers_table
        kw['name'] = name
        super(KmerCache, self).__init__(init_script=init_script, **kw)

    @property
    def kmers_table(self):
        """The kmer hits table name ``seq_kmers_[name]``, cf.
        :attr:`KmerDBWrapper.name`."""
        return 'seq_kmers_%s' % self.name

    def cached_seqs(self):
        """Returns content identifiers for all cached sequences."""
        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute('SELECT seq from %s' % self.kmers_table)
            return [x[0] for x in cursor]

    def as_kmer_seq(self, seq):
        """Return the integer representation of a given sequence.

        Args:
            seq (sequence.Sequence): input sequence.

        Returns:
            list: list of integers of length ``n-w+1`` containing kmers in
                input sequence represented as an integer, cf.
                :func:`as_kmer_seq`.
        """
        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(
                'SELECT kmers from %s WHERE seq = ?' % self.kmers_table,
                (seq.content_id,)
            )
            for kmer_seq in cursor:
                return eval(kmer_seq[0])
        # cache miss, translate to kmer sequence
        self.log('producing kmer representation for sequence %s' %
                 seq.content_id[:8])
        kmer_seq = [i for i in as_kmer_seq(seq, self.wordlen, mask=self.mask)
                    if i is not None]
        with self.connection() as conn:
            q = 'INSERT INTO %s (seq, kmers) VALUES (?, ?)' % self.kmers_table
            conn.cursor().execute(q, (seq.content_id, repr(kmer_seq)))
        return kmer_seq


class KmerIndex(KmerDBWrapper):
    """An index backed by SQLite for occurences of kmers in a body of
    sequences. Upon initialization the following script is executated:

    .. code-block:: sql

        CREATE TABLE kmers_[name] (
          'kmer'  INTEGER,      -- The kmer in integer representation.
          'seq'   INTEGER,      -- integer identifier of sequence
          'pos'   INTEGER       -- the position of kmer in sequence.
        );

        CREATE TABLE IF NOT EXISTS kmer_indexed_[name] (
          'seq'  VARCHAR,                           -- content id,
          'seqid' INTEGER PRIMARY KEY AUTOINCREMENT -- integer id.
        );

    Attributes:
        name (str):
        cache (KmerCache): optional :class:`KmerCache` object to use for
            retrieving integer representations of sequences.
    """
    def __init__(self, name='', kmer_cache=None, **kw):
        self.name = name
        init_script = """
            CREATE TABLE IF NOT EXISTS %s (
              'kmer'  INTEGER,      -- the kmer in integer representation.
              'seqid' INTEGER,      -- integer identifier of sequence
                                    -- REFERENCES kmer_indexed(seqid), but not
                                    -- declared to avoid integrity checks
              'pos'   INTEGER       -- the position of kmer in sequence.
            );

            CREATE TABLE IF NOT EXISTS %s (
              'seq'  VARCHAR,                           -- content id,
              'seqid' INTEGER PRIMARY KEY AUTOINCREMENT -- integer id.
            );
        """ % (self.kmers_table, self.log_table)
        kw['name'] = name
        super(KmerIndex, self).__init__(init_script=init_script, **kw)
        if kmer_cache:
            assert isinstance(kmer_cache, KmerCache)
            assert kmer_cache.wordlen == self.wordlen
            assert kmer_cache.alphabet == self.alphabet
            self.kmer_cache = kmer_cache
        else:
            self.kmer_cache = None

    @property
    def kmers_table(self):
        """The kmer hits table name ``kmers_[name]``, cf.
        :attr:`KmerDBWrapper.name`."""
        return 'kmers_' + self.name

    @property
    def log_table(self):
        """The log table name ``kmer_indexed_[name]``, cf.
        :attr:`KmerDBWrapper.name`."""
        return 'kmer_indexed_' + self.name

    def index_kmers(self, seq):
        """ Indexes all kmers observed in the given sequence in
        :attr:`kmers_table`.

        Args:
            seq (sequence.Sequence): The sequence just inserted into the
                database.
            seqid (int): The integer identifier to use for sequence.
        """
        if self.kmer_cache:
            kmer_seq = self.kmer_cache.as_kmer_seq(seq)
        else:
            kmer_seq = as_kmer_seq(seq, self.wordlen, mask=self.mask)

        with self.connection() as conn:
            self.log('indexing %d-mers for sequence %s (%d)' %
                     (self.wordlen, seq.content_id[:8], len(seq)))
            cursor = conn.cursor()
            q = 'SELECT seqid FROM %s WHERE seq = ?' % self.log_table
            for seqid in cursor.execute(q, (seq.content_id,)):
                self.log('sequence %s already indexed, skipping.' %
                         seq.content_id[:8])
                return seqid[0]
            q = """
                INSERT INTO %s (seq) VALUES (?);
                SELECT last_insert_rowid();
            """ % self.log_table
            seqid = cursor.execute(q, (seq.content_id,)).next()[0]
            q = """
                INSERT INTO %s (kmer, seqid, pos)
                VALUES (?,?,?)
            """ % self.kmers_table
            cursor.executemany(q, ((kmer, seqid, pos)
                                   for pos, kmer in enumerate(kmer_seq)))
            return seqid

    def create_sql_index(self):
        """Creates SQL index over the ``kmer`` column of ``kmers`` table."""
        self.log('Creating SQL index for table %s.' % self.kmers_table)
        with self.connection() as conn:
            q = """
                CREATE INDEX IF NOT EXISTS idx_%s ON %s (kmer);
            """ % (self.kmers_table, self.kmers_table)
            conn.cursor().execute(q)

    def hits(self, kmer):
        """Returns all hits of a given kmer in indexed sequences.

        Args:
            kmer (int): kmer of interest.

        Returns:
            list:
                A list of 2-tuples containing sequence ids (int) and positions.
        """
        assert isinstance(kmer, int)
        query = 'SELECT seqid, pos FROM %s WHERE kmer = ?' % self.kmers_table
        with self.connection() as conn:
            return list(conn.cursor().execute(query, (kmer,)))

    def kmers(self):
        """Returns all observed kmers.

        Returns:
            list: list of kmers in integer representation.
        """
        self.create_sql_index()  # FIXME do we need this?
        query = 'SELECT DISTINCT kmer FROM %s' % self.kmers_table
        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(query)
            return [x[0] for x in cursor]

    def drop_data(self):
        """Drop all tables created by this object."""
        with self.connection() as conn:
            conn.cursor().execute('DROP TABLE %s; DROP TABLE %s;' %
                                  (self.kmers_table, self.logs_table))
