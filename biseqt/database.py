# -*- coding: utf-8 -*-
"""
.. wikisection:: overview
    :title: Storing sequence metadata

    The :mod:`biseqt.database` module provides tools for storing sequence
    metadata (and not contents) for analyses spanning many sequences. The main
    entry point is through the :class:`DB` class. For instance if the file
    ``example.fa`` contains::

        > S
        AATCGG

    we can say:

    >>> from biseqt.database import DB
    >>> from biseqt.sequence import Alphabet
    >>> A = Alphabet('ACGT')
    >>> db = DB('example.db', A)
    >>> db.initialize()
    >>> with open('example.fa') as f:
    ...     db.load_fasta(f)
    >>> for record in db.find(condition=lambda r: r.attrs['name'] == 'S'):
    ...     print(record)
    ...     print(db.load_from_record(record))
    Record(id=1, content_id=u'e690f18fc98d4afcba4b8518b88f4d0387a17380', \
    source_file=u'example.fa', source_pos=0, attrs={u'name': u'S'})
    AATCGG

.. wikisection:: dev
    :title: Insert-or-Append Queries

    You may need to perform an SQL query which either inserts a new record or
    appends the given value to some column in case of conflict. This is achived
    by using one of SQLite's conflict resolution mechanisms_ ``ON CONFLICT
    REPLACE`` or in short ``OR REPLACE``. For instance, consider a table::

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
        cursor = sqlite3.connect(...).cursor()
        cursor.execute(query, (id, id, ',' + field))

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

    * When a table has a unique integer key it should be declared as ``INTEGER
      PRIMARY KEY`` so that it would take over the default ``rowid`` field.
      This saves space (and thus a small amount of time) on both the field and
      the corresponding index.
    * Larger `page sizes <pagesize_docs>`_ can marginally improve read/write
      performance. To increase the page size:

      .. code-block:: sql

        PRAGMA page_size = 65536

    .. _journaling_docs: https://www.sqlite.org/pragma.html#pragma_journal_mode
    .. _pagesize_docs: https://www.sqlite.org/pragma.html#pragma_page_size
    .. rubric: References

        * http://stackoverflow.com/a/1712873
        * http://codereview.stackexchange.com/q/26822
        * http://stackoverflow.com/q/3134900
"""

import os
import sqlite3
import json
import logging
from collections import namedtuple

from . import ProgressIndicator
from .io import read_fasta
from .sequence import Alphabet, NamedSequence


class Record(namedtuple('Record', ['id', 'content_id', 'source_file',
                                   'source_pos', 'attrs'])):
    """A ``namedtuple`` which wraps an SQL record from the ``sequence`` table
    of a :class:`DB`. All fields are intended to map directly to columns of the
    table; cf. :func:`DB.initialize`.

    Attributes:
        id (int): The integer identifier assigned by the database.
        content_id (str): The content identifier of the sequence, cf.
            :attr:`NamedSequence.content_id
            <biseqt.sequence.NamedSequence.content_id>`
        source_file (str): The path to the file where this sequence was read.
        source_pos (int): The position in the file where this sequence begins.
        attrs (dict): A dict of arbitrary attributes (stored as JSON and
            unpacked upon load from database).

    .. wikisection:: dev
        :title: Documenting namedtuples

        Namedtuples do not have docstrings by default. There are a variety of
        hacks to make them work with ``help()`` and sphinx.  The following is
        the cleanest way of documenting namedtuple classes in python 2.7 and
        is used throughout ``biseqt``:

        .. code-block:: python

            from collections import namedtuple

            class Example(namedtuple('Example', ['field1', 'field2'])):
                \"\"\"docstring for Example. \"\"\"

        Note that this kills ``__slots__`` but presumably our performance will
        remain unaffected with reasonable usage patterns (i.e as long as we
        don't start writing to ``__dict__`` and only use the namedtuple
        fields). If found necessary, redeclaring ``__slots__`` is necessary
        since it is not inherited by default from ``tuple``:

        .. code-block:: python

            from collections import namedtuple

            class Example(namedtuple('Example', ['field1', 'field2'])):
                \"\"\"docstring for Example. \"\"\"
                __slots__ = ()

        Note that using ``pass`` in the class body voids the docstring.

        .. rubric:: References

        * http://stackoverflow.com/a/1606478
        * http://bugs.python.org/issue16669

    """


class DB(object):
    """Wraps an SQLite database containing sequences and related information.
    This class is responsible for maintaining and initializing database files,
    as well as populating the ``sequence`` table

    Attributes:
        path (str): Path to the SQLite datbase.
        alphabet (Alphabet): The alphabet for sequences in the database.
    """

    events = ['initialize', 'insert-sequence']
    """Events emitted upon special events, cf. :func:`add_event_listener` and
    :func:`emit`. Currently supported events are:

        * ``initialize(db, conn)``: emitted when the initialization script has
          executed; use this event to execute other initialization scripts that
          create tables or configure the database in an idempotent manner.
        * ``insert-sequence(db, conn, seq, rec)``: emitted after a sequence is
          inserted in the sequence table; use this event to perform further
          processing on arriving sequences.

    Example usage:
        >>> from biseqt.database import DB
        >>> from biseqt.sequence import Alphabet
        >>> def callback(db, conn): print('called back')
        >>> db = DB('example.db', Alphabet('ACGT'))
        >>> db.add_event_listener('initialize', callback)
        >>> db.initialize()
        'called back'
    """

    def __init__(self, path, alphabet, log_level=logging.INFO):
        assert isinstance(alphabet, Alphabet)
        self.alphabet = alphabet
        self.processors = {event: [] for event in self.events}

        path = os.path.abspath(path)
        if os.path.exists(path):
            assert os.access(path, os.W_OK), 'Database %s not writable' % path
        else:
            assert os.access(os.path.dirname(path), os.W_OK), \
                'Database %s cannot be created' % path

        self.path = path
        # names of non-id fields of Record
        self._update_fields = [f for f in Record._fields if f != 'id']
        # can only conflict (and thus replace) if content_id already exists
        # which means we will override its metadata.
        self.insert_q = 'INSERT OR REPLACE INTO sequence (%s) VALUES (%s)' % (
            ','.join(self._update_fields),
            ','.join('?' for _ in self._update_fields)
        )

    _init_script = """
    -- Database initialization script

    PRAGMA journal_mode = OFF; -- turn off journaling for performance, this
                               -- means the contents of database cannot be
                               -- trusted after an unexpected crash.

    CREATE TABLE IF NOT EXISTS sequence (
      id            INTEGER PRIMARY KEY ASC, -- internal integer identifier
      content_id    VARCHAR UNIQUE, -- hash of sequence contents
      source_file   VARCHAR, -- path to file containing the sequence
      source_pos    INT,     -- file position where the sequence begins
      attrs         VARCHAR  -- arbitrary attributes in JSON format
    );
    """

    def initialize(self):
        """Initialize the database and emit the ``initialize`` event (cf.
        :attr:`events`). The initialization operation is idempotant, i.e
        initializing an initialized database has no side effect.
        """
        with self.connect() as conn:
            conn.cursor().executescript(self._init_script)
            self.emit('initialize', conn)

    # put the initialization script in the docs
    initialize.__doc__ += '\n\n\t.. code-block:: sql\n\t%s\n' % \
                          '\n\t'.join(_init_script.split('\n'))

    def connect(self):
        """Provides a context manager for an SQLite database connection:

        Returns:
            sqlite3.Connection

        ::

            >>> from biseqt.database import DB
            >>> with DB('example.db').connect() as conn:
            ...     cursor = conn.cursor()
            ...     cursor.execute('SELECT * FROM sequence')
        """
        return sqlite3.connect(self.path)

    def record_to_row(self, record):
        """Converts a :class:`Record` to a tuple that can be inserted in the
        sequence table."""
        row = []
        for field in self._update_fields:
            value = getattr(record, field)
            if field == 'attrs':
                value = json.dumps(value)
            row.append(value)
        return tuple(row)

    def insert(self, seq, source_file=None, source_pos=0, attrs={}):
        """Inserts a sequence in the database and emits ``insert-sequence``
        (cf. :attr:`events`). If the :attr:`content address
        <Record.content_id>` of the sequence matches that of an existing
        sequence, the old record will be overwritten.

        Args:
            seq (sequence.NamedSequence): The sequence to be inserted.
            source_file(str): Used to populate :attr:`Record.source_file`.
            source_pos (int): Used to populate :attr:`Record.source_pos`.
            attrs (dict): populates :attr:`Record.attrs`; default is a
                dictionary containing only the name of the sequence (name is
                always added if not present in the dictionary).

        Returns:
            Record:
                The record inserted in the database with its :attr:`Record.id`
                populated with the newly assigned primary key.
        """
        assert isinstance(seq, NamedSequence) and seq.alphabet == self.alphabet
        kw = {
            'source_file': source_file,
            'source_pos': source_pos,
            'content_id': seq.content_id,
            'attrs': {k: attrs[k] for k in attrs},
        }
        if 'name' not in kw['attrs']:
            kw['attrs']['name'] = seq.name

        rec = Record(id=None, **kw)
        with self.connect() as conn:
            cursor = conn.execute(self.insert_q, self.record_to_row(rec))
            conn.commit()
            # populate the id of the record
            rec = Record(id=cursor.lastrowid,
                         **{k: getattr(rec, k) for k in self._update_fields})
            self.emit('insert-sequence', conn, seq, rec)
        return rec

    def load_fasta(self, f, num=-1, rc=False):
        """Populates the database from an open file containing sequences in
        FASTA format.

        Args:
            f (file): Open file to read from; passed as is to
                :func:`read_fasta <biseqt.io.read_fasta>`.
            num (int): Number of sequences to read from ``f``; passed as is to
                :func:`read_fasta <biseqt.io.read_fasta>`.
            rc (bool): Whether to also include the reverse complement of each
                sequence read from ``f``; default is False and is only allowed
                to be True for DNA sequences.

        Returns:
            list: The :class:`Record` objects corresponding to inserted
                sequences in the database with their :attr:`Record.id`
                populated.
        """
        try:
            path = os.path.abspath(f.name)
        except AttributeError:
            # e.g. StringIO
            path = None

        self.logger.info('Loading sequences from %s' % str(path))
        quiet = self.logger.getEffectiveLevel() > logging.INFO
        indic = ProgressIndicator(num_total=(num if num > 0 else None),
                                  quiet=quiet)
        indic.start()

        if rc:
            assert all(l in self.alphabet for l in 'ACGT')

        inserted = []
        # FIXME what if the source contains reverse complements?
        # similarly, what to do with duplicate records (e.g. reloading a
        # half-loaded fasta file). Note that KmerIndex expects that we only
        # emit insert-sequence when new stuff come in.
        for seq, pos in read_fasta(f, self.alphabet, num=num):
            kw = {'source_file': path, 'source_pos': pos}
            inserted.append(self.insert(seq, **kw))
            if rc:
                seq_rc = seq.reverse().transform(['AT', 'CG'],
                                                 name='(rc) ' + seq.name)
                kw['attrs'] = {'rc_of': seq.content_id}
                inserted.append(self.insert(seq_rc, **kw))
            indic.progress()

        indic.finish()

        return inserted

    def find(self, condition=None, sql_condition=None):
        """Loads :class:`Record` objects satisfying the given conditions.
        Conditions can be specified either as python filtering callables or an
        SQL ``WHERE`` clause.

        Args:
            condition(callable): A python callable that determines whether a
                record should be yielded; default is None which means no
                filtering.
            sql_condition(str): The body of an SQL ``WHERE`` clause; default is
                None which means no filtering.
        """
        q = 'SELECT %s FROM sequence' % ', '.join(Record._fields)
        q += ' WHERE ' + sql_condition if sql_condition else ''

        def _record_factory(cursor, row):
            kw = {Record._fields[i]: row[i] for i in range(len(row))}
            kw['attrs'] = json.loads(kw['attrs'])
            return Record(**kw)

        with self.connect() as conn:
            conn.row_factory = _record_factory
            cursor = conn.execute(q)
            for record in cursor:
                if not condition or condition(record):
                    yield record

    def load_from_record(self, record):
        """Loads a sequence from the original source file given a corresponding
        database :class:`Record`.

        Args:
            record (Record): Record from the sequence table containing an
                accessible :attr:`source_file <Record.source_file>` and
                the right :attr:`source_pos <Record.source_pos>`. If
                :attr:`attrs <Record.attrs>` indicates that the given record
                belongs to a reverse complement, the reverse complement of the
                loaded sequence is returned.
        """
        source_file = record.source_file
        with open(source_file) as f:
            f.seek(record.source_pos)
            seq, pos = next(read_fasta(f, self.alphabet, num=1))
            assert pos == record.source_pos
            if 'rc_of' in record.attrs:
                assert record.attrs['rc_of'] == seq.content_id
                seq = seq.reverse().transform(['AT', 'CG'],
                                              name=record.attrs['name'])
                assert record.content_id == seq.content_id
            assert record.content_id == seq.content_id
            return seq

    def add_event_listener(self, event, func):
        """Registers a callback for the given event.

        Args:
            event (str): A string in :attr:`events`.
            func (callable): The callable to be invoked, for argument list for
                each event see :attr:`events`.
        """
        assert event in self.events
        self.processors[event].append(func)

    def emit(self, event, *args):
        """Emits an event by executing all registered callbacks for it.

        Args:
            event (str): A string in :attr:`events`.
            args (list): The arguments to pass to the callback.
        """
        for func in self.processors[event]:
            func(*args)
