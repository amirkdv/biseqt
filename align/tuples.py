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
        S_id (int): The id of the "from" sequence as found in ``seq``.
        T_id (int): The id of the "to" sequence as found in ``seq``.
        tx (pw.Transcript):  object representing the alignment.
    """


class TuplesDB(object):
    """Wraps an SQLite database containing tuple indices for sequences. For
    all indexing and querying refer to :class:`Index`. For now,
    only one word length is allowed. The sequences are stored in a ``seq``
    table with the following schema:

    .. code-block:: sql

        CREATE TABLE seq (
          'id' integer PRIMARY KEY ASC,
          'name' text,
          'description' text,
          'seq' text
        );

    Attributes:
        db (string): Path to the SQLite datbase.
    """
    def __init__(self, db, alphabet=None):
        assert isinstance(alphabet, seq.Alphabet)
        self.alphabet, self.db = alphabet, db

    def initdb(self):
        """Initializes the database: creates the required tables and
        fails if any of them already exists."""
        sys.stderr.write('Initializing tuples DB at: %s\n' % self.db)
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

    def populate(self, fasta_src, lim=-1):
        """Given a FASTA source file, loads all the sequences (up to a limit,
        if specified) into the ``seq`` table.

        Args:
            fasta_src(str): Path to FASTA source.
            lim (Optional[int]): If positive, will be the number of sequences
                loaded from FASTA source, default is -1.
        """
        sys.stderr.write('Populating tuples DB at: %s\n' % self.db)

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
    # FIXME docs
    """The main responsibility of an Index is to respond to :func:`seeds`
    with a list of :class:`Segment` s given two sequence IDs.
    This is done by recording all "tuples" observed in all sequences
    in the database during :func:`index`. Each index is uniquely defined by
    its :attr:`tuplesdb` and its :attr:`wordlen`. Each such index stores all
    hits in a single database table with the following schema (replace `N`
    with :attr:`wordlen`):

    .. code-block:: sql

        CREATE TABLE tuples_N (
          'tuple' char(N)
          'hits'  varchar,
          UNIQUE(tuple)
        )
        CREATE TABLE seeds_N (
          'S_id' integer REFERENCES seq(id),
          'T_id' integer REFERENCES seq(id),
          'S_idx' integer,
          'T_idx' integer,
          UNIQUE(S_id, T_id, S_idx, T_idx)
        );

    The ``hits`` column of ``tuples_N`` contains a single string with the
    format ``(@<id>:<idx>)*``. For example, a hit at position 12 of sequence
    with ID 57 (according to the ``seq`` table) is represented by ``@57:12``.

    Attributes:
        tuplesdb (tuplesDB): The tuples database.
        wordlen (int): Length of tuples.
    """
    def __init__(self, tuplesdb, wordlen):
        self.tuplesdb = tuplesdb
        self.wordlen = wordlen
        self.tuples_table = 'tuples_%d' % self.wordlen
        self.seeds_table = 'seeds_%d' % self.wordlen

    def initdb(self):
        """Initializes the database: creates the required tables and
        fails if any of them already exists."""
        with sqlite3.connect(self.tuplesdb.db) as conn:
            c = conn.cursor()
            q = """
                CREATE TABLE %s (
                  'tuple' char(%d),
                  'hits'  varchar,
                  UNIQUE(tuple)
                );
            """ % (self.tuples_table, self.wordlen)
            c.execute(q)
            q = """
                CREATE TABLE %s (
                  'S_id' integer REFERENCES seq(id),
                  'T_id' integer REFERENCES seq(id),
                  'S_idx' integer,
                  'T_idx' integer,
                  UNIQUE(S_id, T_id, S_idx, T_idx)
                );
            """ % (self.seeds_table)
            c.execute(q)
            c.execute(
                'CREATE INDEX seeds_ids ON %s (S_id, T_id)' % self.seeds_table
            )
        sys.stderr.write('Initialized index table %s at: %s\n'
                         % (self.tuples_table, self.tuplesdb.db))

    def tup_scan(self, string):
        """A generator for ``(string, idx)`` tuples to scan through any given
        string. For example::

            string = 'ACGTGT'
            tup_scan(string, 5) # => ('ACGTG', 0), ('CGTGT', 1)

        Args:
            string (str): The string to scan.
            wordlen(int): Length of the words.

        Yields:
            tuple: A string of length :attr:`wordlen` and a starting position.
        """
        for idx in range(len(string) - self.wordlen + 1):
            yield (string[idx:idx + self.wordlen], idx)

    # Helper for index(): yields data values to be inserted in the seeds index.
    def _give_seeds(self, cursor):
        # count the total number of tuples so we can report percentage
        # progress:
        cursor.execute('SELECT count(*) FROM %s' % self.tuples_table)
        for row in cursor:
            num_tuples = int(row[0])
            break
        indicator = ProgressIndicator(
            'indexing %d seeds' % num_tuples, num_tuples
        )
        indicator.start()

        cursor.execute('SELECT tuple, hits from %s' % self.tuples_table)
        for tup, hits in cursor:
            # report progess:
            indicator.progress()

            # hits contains multiple entries of the form: "@<seqid>:<idx>".
            assert(hits[0] == '@')
            hits = [tuple(hit.split(':')) for hit in hits[1:].split('@')]
            hits = [(hit[0], int(hit[1])) for hit in hits]
            # not interested in seeds from a sequence to itself:
            if len(set([hit[0] for hit in hits])) < 2:
                continue
            for S_hit_idx in range(len(hits)):
                for T_hit_idx in range(S_hit_idx + 1, len(hits)):
                    S_hit, T_hit = hits[S_hit_idx], hits[T_hit_idx]
                    yield (S_hit[0], T_hit[0], S_hit[1], T_hit[1])

        indicator.finish()

    def index(self):
        """Scans all sequences in the ``seq`` table and records all observed
        tuples in ``tuples_N`` table. Then scans all the observed hits to
        populate the ``seeds_N`` table used by :func:`seeds`. Each "seed"
        record is a 4-tuple ``(S_id, T_id, S_idx, T_idx)`` where the first two
        entries are ID's of sequences (as per the ``seq`` table) and the last
        two entries are integers corresponding to the starting position of the
        seed in the sequence pair.
        """
        hit_ins_q = """
            INSERT OR REPLACE INTO %s (tuple, hits)
            SELECT ?, IFNULL( (SELECT hits FROM %s WHERE tuple = ?), "") || ?
        """ % (self.tuples_table, self.tuples_table)
        with sqlite3.connect(self.tuplesdb.db) as conn:
            # We need multiple cursors: one two read from the seq table
            # (seq_c), one two read/write from/to the tuples table (tuples_c),
            # and one to write to the seeds table (seeds_c).
            tuples_c = conn.cursor()
            seq_c = conn.cursor()
            seeds_c = conn.cursor()

            # count the total number so we can report percentage progress:
            seq_c.execute('SELECT count(*) FROM seq')
            for row in seq_c:
                num_seqs = int(row[0])
                break
            indicator = ProgressIndicator(
                'indexing %d sequences' % num_seqs, num_seqs
            )
            indicator.start()

            # populate the tuples table:

            def _give_tuple():
                for s, idx in self.tup_scan(string):
                    yield (s, s, '@%s:%d' % (seqid, idx))

            seq_c.execute('SELECT id, seq FROM seq')
            for seqid, string in seq_c:
                tuples_c.executemany(hit_ins_q, _give_tuple())
                indicator.progress()

            indicator.finish()

            # populate the seeds table:
            seed_ins_q = """
                INSERT INTO %s (S_id, T_id, S_idx, T_idx)
                VALUES (?, ?, ?, ?)
            """ % self.seeds_table
            seeds_c.executemany(seed_ins_q, self._give_seeds(tuples_c))

    def seeds(self, S_id, T_id):
        """Given two sequence ids, finds all maximal exactly matching segments
        (see :class:`Segment`) between the two. A segment is in its maximal
        form if it cannot be extended in either direction by an exact match.

        Note:
            Scores of transcripts for exact matches are left as 0 to avoid
            unnecessary cycles.
        """
        q = """
            SELECT S_id, T_id, S_idx, T_idx FROM %s
            WHERE S_id = ? AND T_id = ? OR S_id = ? AND T_id = ?
        """ % (self.seeds_table)
        exacts = []
        with sqlite3.connect(self.tuplesdb.db) as conn:
            c = conn.cursor()
            c.execute(q, (S_id, T_id, T_id, S_id))
            for row in c:
                if row[0] == S_id and row[1] == T_id:
                    S_idx, T_idx = row[2], row[3]
                elif row[0] == T_id and row[1] == S_id:
                    T_idx, S_idx = row[2], row[3]
                else:
                    raise RuntimeError("This should not have happend!")
                tx = pw.Transcript(S_idx, T_idx, 0, 'M'*self.wordlen)
                exacts += [Segment(S_id=S_id, T_id=T_id, tx=tx)]

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
                        0,  # score
                        exacts[cand].tx.opseq + 'M'*shift_S  # transcript
                    )
                    exacts[idx] = Segment(S_id=S_id, T_id=T_id, tx=tx)
                    exacts.pop(cand)
                else:
                    cand += 1
            idx += 1
        exacts.sort(key=lambda s: len(s.tx.opseq), reverse=True)
        return exacts
