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
    ...     for record in db.find(condition=lambda r: r.attrs['name'] == 'S'):
    ...         print(record)
    ...         print(db.load_from_record(record, f))
    Record(id=1, content_id=u'e690f18fc98d4afcba4b8518b88f4d0387a17380', \
    source_file=u'example.fa', source_pos=0, attrs={u'name': u'S'})
    AATCGG

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

import os
import apsw
import json
import logging
import textwrap
from collections import namedtuple

from .util import Logger, ProgressIndicator
from .sequence import Alphabet, Sequence


class NamedSequence(Sequence):
    """A named version of :class:`biseqt.sequence.Sequence`. The main
    differences are the additional attributes which allow for faster equality
    comparison for large sequences.

    Attributes:
        name (str): The name of the sequence which need not be unique.
    """
    def __init__(self, seq, name=''):
        super(NamedSequence, self).__init__(seq.alphabet, seq.contents)
        self.name = name

    def reverse(self, name=None):
        """Wraps :func:`Sequence.reverse` to make sure a named sequence is
        returned.

        Args:
            name(str): The name to give to the new sequence. Default is None
                in which case ``(reversed)`` is preprended to the original
                sequence name.

        Returns:
            NamedSequence
        """
        rev = super(NamedSequence, self).reverse()
        if name is None:
            name = '(reversed) ' + self.name
        return NamedSequence(rev, name=name)

    def transform(self, mappings={}, name=None):
        """Wraps :func:`Sequence.transform` to make sure a named sequence is
        returned.

        Args:
            mappings (dict|list): As in :func:`Sequence.transform`.
            name (str): The name to give to the new sequence. Default is None
                in which case ``(transformed)`` is prependended to the
                original sequence name.

        Returns:
            NamedSequence
        """
        seq = super(NamedSequence, self).transform(mappings=mappings)
        if name is None:
            name = '(transformed) ' + self.name
        return NamedSequence(seq, name=name)

    def __eq__(self, other):
        return isinstance(other, NamedSequence) and \
               self.content_id == other.content_id and \
               self.name == other.name

    def __repr__(self):
        return 'NamedSequence(Sequence(alphabet=%s, contents=%s), name=%s)' % (
            repr(self.alphabet),
            repr(self.contents),
            repr(self.name),
        )

    def to_fasta(self, width=80):
        """Sequence in FASTA format.

        Returns:
            str
        """
        contents = str(self)
        if width is not None:
            contents = '\n'.join(textwrap.wrap(contents, width))
        return '>%s\n%s\n' % (self.name, contents)


class Record(namedtuple('Record', ['id', 'content_id', 'source_file',
                                   'source_pos', 'attrs'])):
    """A ``namedtuple`` which wraps an SQL record from the ``sequence`` table
    of a :class:`DB`. All fields are intended to map directly to columns of the
    table; cf. :func:`DB.initialize`.

    Attributes:
        id (int): The integer identifier assigned by the database.
        content_id (str): The content identifier of the sequence, cf.
            :attr:`NamedSequence.content_id`.
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


def read_fasta(f, alphabet, num=-1):
    """Lazy generator for content-identifiable :class:`NamedSequence` objects
        loaded from a FASTA source.

    Args:
        f (file): A readable open file or equivalent (e.g.
            :class:`StringIO.StringIO`). The file position is only modified by
            reading; namely, sequences are read starting from the current
            position of the file and after this function returns the file
            position points at the last character of the last yielded sequence.
        alphabet (sequence.Alphabet): The alphabet of the sequences.
        num (int): The number of sequences to read from the file; default is -1
            in which case all sequences are read until end-of-file is reached.

    Yields:
        tuple:
            A :class:`NamedSequence` whose
            :attr:`NamedSequence.name` is the FASTA
            record name, and the starting position of the sequence in ``f`` as
            an integer.
    """
    observed_names = []

    def _parse(raw_seq, name):
        if cur_seq and cur_name:
            assert cur_name not in observed_names, \
                'Duplicate sequence name: %s' % cur_name
            observed_names.append(cur_name)
            return NamedSequence(alphabet.parse(raw_seq), name=cur_name)

    count = cur_pos = 0
    cur_name = cur_seq = ''
    # NOTE `for raw_line in f` uses a read-ahead buffer which makes `f.tell()`
    # useless for remembering where a sequence begins.
    # cf. https://docs.python.org/2/library/stdtypes.html#file.next
    for raw_line in iter(f.readline, ''):
        line = raw_line.strip()
        if not line:
            continue
        if line[0] == '>':
            seq = _parse(cur_seq, cur_name)
            if seq:
                yield seq, cur_pos
                count += 1
                if num > 0 and count >= num:
                    return
            cur_pos = f.tell() - len(raw_line)
            cur_name = line[1:].strip()
            cur_seq = ''
        else:
            cur_seq += line
    seq = _parse(cur_seq, cur_name)
    if seq:
        yield seq, cur_pos


class DB(object):
    """Wraps an SQLite database containing sequences and related information.
    This class is responsible for maintaining and initializing database files,
    as well as populating the ``sequence`` table

    Attributes:
        path (str): Path to the SQLite datbase.
        alphabet (Alphabet): The alphabet for sequences in the database.
    """

    events = ['db-initialized', 'sequence-inserted', 'sequences-loading']
    """Events emitted upon special events, cf. :func:`add_event_listener` and
    :func:`emit`. Currently supported events are:

        * ``db-initialized(conn)``: emitted after the initialization script has
          executed; use this event to execute other initialization scripts that
          create tables or configure the database in an idempotent manner.
        * ``sequences-loading(f)``: emitted before bulk loading sequences.
        * ``sequence-inserted(conn, seq, rec)``: emitted after a sequence is
          inserted in the sequence table; use this event to perform further
          processing on arriving sequences.

    Example usage:
        >>> from biseqt.database import DB
        >>> from biseqt.sequence import Alphabet
        >>> def callback(db, conn): print('called back')
        >>> db = DB('example.db', Alphabet('ACGT'))
        >>> db.add_event_listener('db-initialized', callback)
        >>> db.initialize()
        'called back'
    """

    def __init__(self, path, alphabet, log_level=logging.INFO):
        assert isinstance(alphabet, Alphabet)
        self.alphabet = alphabet
        self.processors = {event: [] for event in self.events}

        if path == ':memory:':
            self.path = path
        else:
            self.path = os.path.abspath(path)
            assert os.path.exists(self.path) or \
                os.access(os.path.dirname(self.path), os.W_OK), \
                'Database %s is not writable' % self.path

        self._logger = Logger(log_level=log_level,
                              header=os.path.relpath(self.path, os.getcwd()))
        self._connection = None

        # names of non-id fields of Record
        self._update_fields = [f for f in Record._fields if f != 'id']
        self.insert_q = """
            INSERT INTO sequence (%s) VALUES (%s);
            SELECT last_insert_rowid() FROM sequence;
        """ % (','.join(self._update_fields),
               ','.join('?' for _ in self._update_fields))

    _init_script = """
    -- Database initialization script

    CREATE TABLE IF NOT EXISTS sequence (
      id            INTEGER PRIMARY KEY ASC, -- internal integer identifier
      content_id    VARCHAR UNIQUE, -- hash of sequence contents
      source_file   VARCHAR, -- path to file containing the sequence
      source_pos    INT,     -- file position where the sequence begins
      attrs         VARCHAR  -- arbitrary attributes in JSON format
    );
    """

    def initialize(self):
        """Initialize the database and emit the ``db-initialized`` event (cf.
        :attr:`events`). The initialization operation is idempotant, i.e
        initializing an initialized database has no side effect.
        """
        with self.connection() as conn:
            conn.cursor().execute(self._init_script)
            self.emit('db-initialized', conn)

    # put the initialization script in the docs
    initialize.__doc__ += '\n\n\t.. code-block:: sql\n\t%s\n' % \
                          '\n\t'.join(_init_script.split('\n'))

    def log(self, *args, **kwargs):
        """Wraps :class:`Logger.log`."""
        self._logger.log(*args, **kwargs)

    def connection(self):
        """Provides a SQLite database connection that can be used as a context
        manager. The returned object is always the same connection object
        belonging to the :class:`DB` instance (otherwise in-memory connections
        would reset the database contents upon every invocation).

        Returns:
            apsw.Connection

        ::

            >>> from biseqt.database import DB
            >>> with DB('example.db').connection() as conn:
            ...     conn.cursor().execute('SELECT * FROM sequence')
        """
        if self._connection is None:
            self._connection = apsw.Connection(self.path)
        return self._connection

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
        """Inserts a sequence in the database and emits ``sequence-inserted``
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
        rec_kw = {
            'source_file': source_file,
            'source_pos': source_pos,
            'content_id': seq.content_id,
            'attrs': {k: attrs[k] for k in attrs},
        }
        if 'name' not in rec_kw['attrs']:
            rec_kw['attrs']['name'] = seq.name

        rec = Record(id=None, **rec_kw)
        with self.connection() as conn:
            cursor = conn.cursor()
            try:
                cursor.execute(self.insert_q, self.record_to_row(rec))
                # populate the id of the record
                rec = Record(id=next(cursor)[0], **rec_kw)
            except apsw.ConstraintError:
                self.log('ignoring duplicate sequence %s' % seq.name,
                         level=logging.WARN)
                return None
            self.emit('sequence-inserted', conn, seq, rec)
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

        self.emit('sequences-loading', f)
        self.log('Loading sequences from %s ...' % str(path))
        indic = ProgressIndicator(num_total=(num if num > 0 else None))
        indic.start()

        if rc:
            assert all(l in self.alphabet for l in 'ACGT')

        inserted = []
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
        self.log('Done loading sequences from %s.' % str(path))

        return inserted

    def find(self, condition=None, sql_condition=None):
        """Loads :class:`Record` objects satisfying the given conditions.
        Conditions can be specified either as python filtering callables or an
        SQL ``WHERE`` clause. Records are always yielded in increasing order of
        their database integer id.

        Args:
            condition(callable): A python callable that determines whether a
                record should be yielded; default is None which means no
                filtering.
            sql_condition(str): The body of an SQL ``WHERE`` clause; default is
                None which means no filtering.
        Yields:
            Record
        """
        q = 'SELECT %s FROM sequence' % ', '.join(Record._fields)
        q += ' WHERE ' + sql_condition if sql_condition else ''
        q += ' ORDER BY id ASC'

        def _record_factory(cursor, row):
            kw = {Record._fields[i]: row[i] for i in range(len(row))}
            kw['attrs'] = json.loads(kw['attrs'])
            return Record(**kw)

        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.setrowtrace(_record_factory)
            for record in cursor.execute(q):
                if not condition or condition(record):
                    yield record

    def load_from_record(self, record, f=None):
        """Loads a sequence from the original source file given a corresponding
        database :class:`Record`.

        Args:
            record (Record): Record from the sequence table containing
                the right :attr:`source_pos <Record.source_pos>` in the source
                file. If :attr:`attrs <Record.attrs>` indicates that the given
                record belongs to a reverse complement, the reverse complement
                of the loaded sequence is returned.
            f (file): An open file handle for the :attr:`source_file
                <Record.source_file>` where the record is to be found. Default
                is None in which case the file will be openned and closed here.
        """
        was_open = False
        if f is None:
            was_open = True
            f = open(record.source_file)

        f.seek(record.source_pos)
        seq, pos = next(read_fasta(f, self.alphabet, num=1))
        assert pos == record.source_pos
        if 'rc_of' in record.attrs:
            assert record.attrs['rc_of'] == seq.content_id
            seq = seq.reverse().transform(['AT', 'CG'],
                                          name=record.attrs['name'])
            assert record.content_id == seq.content_id
        assert record.content_id == seq.content_id

        if was_open:
            f.close()

        return seq

    def add_event_listener(self, event, func):
        """Registers a callback for the given event.

        Args:
            event (str): A string in :attr:`events`.
            func (callable): The callable to be invoked, for argument list for
                each event see :attr:`events`.
        """
        assert event in self.events, 'unrecognized event %s' % event
        self.processors[event].append(func)

    def emit(self, event, *args):
        """Emits an event by executing all registered callbacks for it.

        Args:
            event (str): A string in :attr:`events`.
            args (list): The arguments to pass to the callback.
        """
        for func in self.processors[event]:
            func(*args)
