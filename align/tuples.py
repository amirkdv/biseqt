"""Provides tools for *k*-mer analysis."""

import sys
import re
import sqlite3
from Bio import SeqIO
from collections import namedtuple

from . import pw, seq, ProgressIndicator
class Segment(namedtuple('Segment', ['S_id', 'T_id', 'tx'])):
    """Represents an aligned pair of substrings in two sequences. The alignment
    may potentially contain indels. Maximal, exactly-matching segments are
    refered to as "seeds".

    Attributes:
        S_id (int): The id of the "from" sequence as found in the ``seq`` table.
        T_id (int): The id of the "to" sequence as found in the ``seq`` table.
        tx (pw.Transcript):  object representing the alignment.
    """

class TuplesDB(object):
    """Wraps an SQLite database containing tuple indices for sequences. For now,
    only one word length is allowed. The tables are:

        * ``seq`` contains all the provided
          sequences (the "database"), each sequence and its name are mapped to
          an internal integer used for calculations (see :func:`populate`).
        * ``tuples_N`` contains all observed words of length ``N``
          where ``N`` is :attr:`wordlen`.
        * ``tuples_N_hits`` contains all the "hits" (see :func:`index`).

    SQL schemas are as follows (replace `N` with word length):

    .. code-block:: sql

        CREATE TABLE seq (
          'id' integer PRIMARY KEY ASC,
          'name' text,
          'description' text,
          'seq' text
        );
        CREATE TABLE tuples_N (
          'id' integer PRIMARY KEY ASC,
          'tuple' char(N),
          UNIQUE(tuple)
        );
        CREATE TABLE tuples_N_hits (
          'tuple' integer REFERENCES tuples_N(id),
          'seq'   integer REFERENCES seq(id),
          'idx'   integer,
          UNIQUE(seq, idx)
        )

    Attributes:
        db (string): Path to the SQLite datbase.
        wordlen (int): Length of tuples
    """
    def __init__(self, db, alphabet=None):
        assert isinstance(alphabet, seq.Alphabet)
        self.alphabet, self.db = alphabet, db

    def initdb(self):
        """Initializes the database: creates the required tables and
        fails if any of them already exists."""
        with sqlite3.connect(self.db) as conn:
            c = conn.cursor()
            q = """
                CREATE TABLE seq (
                  'id' integer PRIMARY KEY ASC,
                  'name' text,
                  'description' text,
                  'seq' text
                );
            """
            c.execute(q)
        sys.stderr.write('Initialized tuples DB at: %s\n' % self.db)

    def populate(self, fasta_src, lim=-1):
        """Given a FASTA source file, loads all the sequences (up to a limit, if
        specified) into the `seq` table.

        Args:
            fasta_src(str): Path to FASTA source.
            lim (Optional[int]): If positive, will be the number of sequences
                loaded from FASTA source, default is -1.
        """
        def give_seq():
            for idx, seq in enumerate(SeqIO.parse(fasta_src, "fasta")):
                if lim < 0 or idx < lim:
                    yield (str(seq.id), str(seq.description), str(seq.seq))
                else:
                    break
        with sqlite3.connect(self.db) as conn:
            c = conn.cursor()
            q = "INSERT INTO seq (name, description, seq) VALUES (?,?,?)"
            c.executemany(q, give_seq())

    def loadseq(self, seqid):
        """Loads a sequence given its internal numeric seqid.

        Args:
            seqid (int): Sequence ID as found in the ``seq`` table.

        Returns:
            seq.Sequence
        """
        with sqlite3.connect(self.db) as conn:
            c = conn.cursor()
            c.execute('SELECT seq FROM seq WHERE id = ?', (seqid,))
            for row in c:
                return seq.Sequence(str(row[0]), self.alphabet)

    def seqinfo(self):
        """Return a dict of metadata about all sequences keyed by sequence ids
        as found in the ``seq`` table. The output looks like this::

            A = seq.Alphabet('ACGT')
            T = tuples.TuplesDB('genome.db', alphabet=A)
            info = T.seqinfo()
            info[12] #=> {'start': 17, 'length': 479, 'name': u'R12'}

        Returns:
            dict: Metadata dicts in a dict keyed by sequence ID.
        """
        with sqlite3.connect(self.db) as conn:
            c = conn.cursor()
            c.execute("SELECT id, name, LENGTH(seq) FROM seq")
            seqs = {}
            for row in c:
                seqs[row[0]] = {}
                # see seq.make_sequencing_fixture()
                seqs[row[0]]['name'], seqs[row[0]]['start'] = row[1].split('_')
                seqs[row[0]]['start'] = int(seqs[row[0]]['start'][1:])
                seqs[row[0]]['length'] = row[2]
            return seqs

    def seqids(self):
        """Returns a list of all seqids found in the ``seq`` table.

        Returns:
            List[int]
        """
        with sqlite3.connect(self.db) as conn:
            c = conn.cursor()
            c.execute('SELECT id FROM seq')
            return [row[0] for row in c]


class Index(object):
    def __init__(self, tuplesdb, wordlen):
        self.tuplesdb = tuplesdb
        self.wordlen = wordlen
        self.tuples_table = 'tuples_%d' % self.wordlen
        self.hits_table = 'hits_%d' % self.wordlen

    def tup_scan(self, string):
        """A generator for ``(string, idx)`` tuples to scan through any given
        string. For example::

                string = 'ACGTGT'
                tup_scan(string, 5) # => ('ACGTG', 0), ('CGTGT', 1)

        Args:
            string (str): The string to scan.
            wordlen(int): Length of the words.

        Yields:
            tuple: A ``string`` of length ``wordlen`` and a starting position (``int``).
        """
        for idx in range(len(string) - self.wordlen + 1):
            yield (string[idx:idx + self.wordlen], idx)

    def initdb(self):
        """Initializes the database: creates the required tables and
        fails if any of them already exists."""
        with sqlite3.connect(self.tuplesdb.db) as conn:
            c = conn.cursor()
            q = """
                CREATE TABLE %s (
                  'id' integer PRIMARY KEY ASC,
                  'tuple' char(%d),
                  UNIQUE(tuple)
                );
            """ % (self.tuples_table, self.wordlen)
            c.execute(q)
            q = """
                CREATE TABLE %s (
                  'tuple' integer REFERENCES %s(id),
                  'seq'   integer REFERENCES seq(id),
                  'idx'   integer,
                  UNIQUE(seq, idx)
                );
            """ % (self.hits_table, self.tuples_table)
            c.execute(q)
        sys.stderr.write('Initialized index tables %s, %s at: %s\n' %
            (self.tuples_table, self.hits_table, self.tuplesdb.db))

    def index(self):
        """Scans all sequences in the ``seq`` table and records all observed
        tuples in ``tuples_N`` tables and all hits in ``tuples_N_hits`` tables.

        Args:
            hp_condenser (homopolymeic.HpCondenser): If specified, its
                :func:`tup_scan <align.homopolymeric.HpCondenser.tup_scan>` is
                used instead of :func:`this module <tup_scan>`'s generic
                version.

        Note:
            It takes ~3 minutes to index 10-mers of 500 sequences with average
            length 11Kbp into a DB of ~250MB. With 20-mers takes ~4 minutes and
            DB size is ~600MB!
        """
        # We want an Upsert to avoid changing tuple
        # IDs see http://stackoverflow.com/a/4330694 .
        tup_insert_q = """
            INSERT OR REPLACE INTO %s (id, tuple)
            VALUES ((SELECT id FROM %s WHERE tuple = ?), ?)
        """ % (self.tuples_table, self.tuples_table)
        hit_insert_q = """
            INSERT INTO %s (tuple, seq, idx)
                SELECT id, ?, ?  FROM %s WHERE tuple = ?
        """ % (self.hits_table, self.tuples_table)
        sys.stderr.write('indexing sequences:')
        with sqlite3.connect(self.tuplesdb.db) as conn:
            c = conn.cursor()
            # can't nest sqlite queries, pick out IDs first and load each
            # sequence separately:
            c.execute('SELECT id FROM seq')
            ids = [row[0] for row in c]
            for seqid in ids:
                c.execute('SELECT seq FROM seq WHERE id = ?', (seqid,))
                for row in c:
                    string = row[0]
                    break
                c.executemany(
                    tup_insert_q,
                    [(s,s) for s, _ in self.tup_scan(string)]
                )
                c.executemany(
                    hit_insert_q,
                    [(seqid, idx, s) for s, idx in self.tup_scan(string)]
                )
                sys.stderr.write(' ' + str(seqid))
            sys.stderr.write('\n')

    def seeds(self, S_id, T_id):
        """Given two sequence ids, finds all maximal exactly matching segments
        between the two. A segment is in its maximal form if there are no other
        exactly-matching segments whose span is a unit shift in both the query
        and target sequences.

        Note:
            Scores of transcripts for exact matches are left as 0 to avoid
            unnecessary cycles, we populate the score first thing in
            :func:`extend`.
        """
        q = """
            SELECT H1.idx, H2.idx FROM %s H1 -- H1 is S, H2 is T
            INNER JOIN %s H2
            ON H1.tuple = H2.tuple
            WHERE H1.seq = ? AND H2.seq = ?
        """ % (self.hits_table, self.hits_table)
        def row_factory(cursor, row):
            tx = pw.Transcript(row[0], row[1], 0, 'M'*self.wordlen)
            return Segment(S_id=S_id, T_id=T_id, tx=tx)

        with sqlite3.connect(self.tuplesdb.db) as conn:
            conn.row_factory = row_factory
            c = conn.cursor()
            c.execute(q, (S_id, T_id))
            exacts = [r for r in c]

        exacts.sort(key=lambda s: s.tx.idx_S)
        # merge overlapping tuples:
        idx = 0
        while idx < len(exacts):
            cand = idx + 1
            while cand < len(exacts):
                shift_S = exacts[cand].tx.idx_S - exacts[idx].tx.idx_S
                shift_T = exacts[cand].tx.idx_T - exacts[idx].tx.idx_T
                # we know the transcripts are all M's.
                if shift_S == shift_T and shift_S > 0 and \
                   shift_S < len(exacts[idx].tx.opseq):
                    tx = pw.Transcript(
                        exacts[idx].tx.idx_S, exacts[idx].tx.idx_T,
                        0, # score
                        exacts[cand].tx.opseq + 'M'*shift_S # transcript
                    )
                    exacts[idx] = Segment(S_id=S_id, T_id=T_id, tx=tx)
                    exacts.pop(cand)
                else:
                    cand += 1
            idx += 1
        exacts.sort(key=lambda s: len(s.tx.opseq), reverse=True)
        return exacts
