# -*- coding: utf-8 -*-
import os
import sqlite3
import json
from collections import namedtuple

from .io import read_fasta
from .sequence import Alphabet, NamedSequence


# NOTE The following is the cleanest way of documenting namedtuple classes
# in python 2.7.
#   cf. http://stackoverflow.com/a/1606478
#   cf. http://bugs.python.org/issue16669
class Record(namedtuple('Record', ['id', 'content_id', 'source_file',
                                   'source_pos', 'attrs'])):
    """Wraps an SQL record from the ``sequence`` table of a :class:`DB`. All
    fields are intended to map directly to columns of the table; cf.
    :class:`DB`.

    Attributes:
        id (int): The integer identifier assigned by the database.
        content_id (str): The content identifier of the sequence, cf.
            :attr:`NamedSequence.content_id
            <biseqt.sequence.NamedSequence.content_id>`
        source_file (str): The path to the file where this sequence was read.
        source_pos (int): The position in the file where this sequence begins.
        attrs (dict): A dict of arbitrary attributes (stored as JSON and
            unpacked upon load from database).
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
    """Events emitted upon special events, cf. :func:`register` and
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
        >>> db.register('initialize', callback)
        >>> db.initialize()
        'called back'
    """

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

    def __init__(self, path, alphabet):
        assert isinstance(alphabet, Alphabet)
        self.alphabet = alphabet

        path = os.path.abspath(path)
        if os.path.exists(path):
            assert os.access(path, os.W_OK), 'Database %s not writable' % path
        else:
            assert os.access(os.path.dirname(path), os.W_OK), \
                'Database %s cannot be created' % path

        self.path = path
        self._update_fields = [f for f in Record._fields if f != 'id']
        self.insert_q = 'INSERT INTO sequence (%s) VALUES (%s)' % (
            ','.join(self._update_fields),
            ','.join('?' for _ in self._update_fields)
        )

    def initialize(self):
        """Initialize the database and emit the ``initialize`` event (cf.
        :attr:`events`). The initialization operation is idempotant, i.e
        initializing an initialized database has no side effect.
        """
        with self.connect() as conn:
            conn.cursor().executescript(self._init_script)
            self.emit('initialize', self, conn)

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
        """Populates the database from an iterable of :class:`Record` objects
        created by :func:`create_record`.

        Args:
            seq (sequence.NamedSequence): The sequence to be inserted.
            source_file(str): Used to populate :attr:`Record.source_file`.
            source_pos (int): Used to populate :attr:`Record.source_pos`.
            attrs (dict): populates :attr:`Record.attrs`; default is a
                dictionary containing only the name of the sequence (name is
                always added if not present in the dictionary).
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
            self.emit('insert-sequence', self, conn, seq, rec)
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

        if rc:
            assert all(l in self.alphabet for l in 'ACGT')

        inserted = []
        # TODO what if the source contains reverse complements?
        for seq, pos in read_fasta(f, self.alphabet, num=num):
            kw = {'source_file': path, 'source_pos': pos}
            inserted.append(self.insert(seq, **kw))
            if rc:
                seq_rc = seq.reverse().transform(['AT', 'CG'],
                                                 name='(rc) ' + seq.name)
                kw['attrs'] = {'rc_of': seq.content_id}
                inserted.append(self.insert(seq_rc, **kw))

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
