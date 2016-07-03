# -*- coding: utf-8 -*-
import os
import sqlite3
import json
from collections import namedtuple
from itertools import chain

from .io import read_fasta
from .sequence import Alphabet, NamedSequence


# cf. http://stackoverflow.com/a/1606478
# cf. http://bugs.python.org/issue16669
class Record(namedtuple('Record', ['id', 'content_id', 'source_file',
                                   'source_pos', 'attrs'])):
    """FIXME
    Attributes:
        content_id (str): Blah
    """


def create_record(seq, **kw):
    """FIXME
    """
    assert isinstance(seq, NamedSequence)
    assert 'source_file' in kw
    assert 'source_pos' in kw
    if 'attrs' not in kw:
        kw['attrs'] = {}
    if 'name' not in kw:
        kw['attrs']['name'] = seq.name
    return Record(id=None, content_id=seq.content_id, **kw)


class DB(object):
    """Wraps an SQLite database containing sequences and related information.
    This class is responsible for maintaining and initializing database files,
    as well as populating the ``sequence`` table.

    Attributes:
        path (string): Path to the SQLite datbase.
        alphabet (Alphabet): The alphabet for sequences in the database.
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
        """Initializes a database at the given location. The initialization
        operation is idempotant, i.e initializing and initialized database has
        no side effects.
        """
        assert isinstance(alphabet, Alphabet)
        self.alphabet = alphabet

        if os.path.exists(path):
            assert os.access(path, os.W_OK), 'Database %s not writable' % path
        else:
            assert os.access(os.path.dirname(path), os.W_OK), \
                'Database %s cannot be created' % path

        self.path = path
        with self.connect() as conn:
            conn.cursor().executescript(self._init_script)

    # add the initialization SQL query to docstring so we don't have to
    # duplicate it.
    __init__.__doc__ += '\n\n.. code-block:: sql\n' + _init_script

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

    def populate(self, records):
        """Populates the database from an iterable of :class:`Record` objects.

        Args:
            records (iterable): Records to put in the sequence table; their
                :attr:`Record.id` is ignored.
        """
        fields = Record._fields[1:]  # skip id

        def _to_row(record):
            row = []
            for field in fields:
                value = getattr(record, field)
                if field == 'attrs':
                    value = json.dumps(value)
                row.append(value)
            return tuple(row)

        q = 'INSERT INTO sequence (%s) VALUES (%s)' % \
            (','.join(fields), ','.join('?' for _ in fields))
        with self.connect() as conn:
            conn.cursor().executemany(q, (_to_row(r) for r in records))

    def populate_from_fasta(self, f, num=-1, rc=False):
        """Populates the database from an open file containing sequences in
        FASTA format.

        Args:
            f (file): Open file to read from, passed as is to
                :func:`read_fasta <biseqt.io.read_fasta>`.
            num (int): Number of sequences to read from ``f``; passed as is to
                :func:`read_fasta <biseqt.io.read_fasta>`.
            rc (bool): Whether to also include the reverse complement of each
                sequence read from ``f``; default is False and is only allowed
                to be True for DNA sequences.
        """
        try:
            path = os.path.abspath(f.name)
        except AttributeError:
            # e.g. StringIO
            path = None

        if rc:
            assert all(l in self.alphabet for l in 'ACGT')

        # TODO what if the source contains reverse complements?
        def _recs_from_seq(seq, pos):
            rec = create_record(seq, source_file=path, source_pos=pos)
            yield rec
            if rc:
                compl = seq.reverse().transform(['AT', 'CG'],
                                                name='(rc) ' + seq.name)
                yield create_record(compl, source_file=path, source_pos=pos,
                                    attrs={'rc_of': rec.content_id})

        records = chain(*(_recs_from_seq(seq, pos) for seq, pos
                          in read_fasta(f, self.alphabet, num=num)))
        self.populate(records)

    def find(self, condition=None, sql_condition=None):
        """Loads sequence :class:`Record` objects satisfying the given
        conditions. Conditions can be specified either as python filtering
        callables or an SQL ``WHERE`` clause.

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
