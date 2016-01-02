"""Provides tools for *k*-mer analysis."""

import sys
import re
import sqlite3
from Bio import SeqIO
from collections import namedtuple
from . import pw, seq, ProgressIndicator, CffiObject, ffi, lib


class Segment(CffiObject):
    """Wraps a C ``segment``: represents an aligned pair of substrings in
    two sequences.

    Attributes:
        S_id (int): The id of the "from" sequence as found in ``seq``.
        T_id (int): The id of the "to" sequence as found in ``seq``.
        tx (align.pw.Transcript): The alignment transctipt.
    """
    def __init__(self, **kw):
        if 'c_obj' in kw:
            self.c_obj = kw['c_obj']
            self.tx = pw.Transcript(c_obj=kw['c_obj'].tx)
        else:
            self.tx = kw['tx']
            self.c_obj = ffi.new('segment*', {
                'S_id': kw['S_id'],
                'T_id': kw['T_id'],
                'tx': kw['tx'].c_obj,
            })

    def __repr__(self):
        return 'Segment(S_id=%d,T_id=%d,tx=%s)' \
            % (self.S_id, self.T_id, self.tx)


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
            sys.stderr.write("Turning off SQLite's journaling for performance: do not trust the contents of the database after a crash.\n")
            c.execute('PRAGMA journal_mode = OFF;')
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
        sys.stderr.write('Loading sequences from: %s\n' % fasta_src)

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
    """The main responsibility of an Index is to respond to :func:`seeds`
    with a list of :class:`Segment` s given two sequence IDs.
    This is done by recording all "tuples" observed in all sequences
    in the database during :func:`index`. Each index is uniquely defined by
    its :attr:`tuplesdb` and its :attr:`wordlen`. Each such index operates
    using 3 database tables with the following schema (replace `N`
    with :attr:`wordlen`). All tables are populated, in the order presented
    below, by :func:`index`.

    .. code-block:: sql

        CREATE TABLE tuples_N (
          'tuple' integer, -- tuple of length N represented in base 3
          'hits'  varchar, -- '@' delimited string of hits with format (@<id>:<idx>)*
          UNIQUE(tuple)
        )
        CREATE TABLE seeds_N (
          -- each record corresponds to a single exactly matching segment
          'S_id' integer REFERENCES seq(id),
          'T_id' integer REFERENCES seq(id),
          'S_idx' integer, -- starting position of seed in S
          'T_idx' integer, -- starting position of seed in T
        );
        CREATE TABLE potential_homologs_N (
          'id' integer REFERENCES seq(id),
          'homologs' varchar, -- comma separated list of sequence IDs, with
                              -- potential homology to sequence with ID 'id'.
        );

    Attributes:
        tuplesdb (tuplesDB): The tuples database.
        wordlen (int): Length of tuples.
        min_seeds_for_homology (int): Minimum number of seeds between two
            sequences that makes them "potential homologs"; default is 1.
    """
    def __init__(self, tuplesdb, wordlen, min_seeds_for_homology=1):
        self.tuplesdb = tuplesdb
        self.wordlen = wordlen
        self.min_seeds_for_homology = min_seeds_for_homology
        self.tuples_table = 'tuples_%d' % self.wordlen
        self.seeds_table = 'seeds_%d' % self.wordlen
        self.potential_homologs_table = 'potential_homologs_%d' % self.wordlen
        self.potential_homologs_q = """
            SELECT S_id, T_id, COUNT(*) AS count
            FROM %s
            GROUP BY S_id, T_id
            HAVING count > %d
        """ % (self.seeds_table, self.min_seeds_for_homology)

    def initdb(self):
        """Initializes the database: creates the required tables and
        fails if any of them already exists."""
        with sqlite3.connect(self.tuplesdb.db) as conn:
            c = conn.cursor()
            q = """
                CREATE TABLE %s (
                  'tuple' integer primary key,
                  'hits'  varchar
                );
            """ % (self.tuples_table)
            c.execute(q)
            q = """
                CREATE TABLE %s (
                  'S_id' integer REFERENCES seq(id),
                  'T_id' integer REFERENCES seq(id),
                  'S_idx' integer,
                  'T_idx' integer
                );
            """ % (self.seeds_table)
            c.execute(q)
            q = """
                CREATE TABLE %s (
                  'id' integer REFERENCES seq(id),
                  'homologs' varchar,
                  UNIQUE(id)
                );
            """ % (self.potential_homologs_table)
            c.execute(q)
        sys.stderr.write('Initialized index tables %s, %s, %s.\n'
            % (self.tuples_table, self.seeds_table, self.potential_homologs_table)
        )

    def tup_scan(self, string):
        """A generator for ``(word, idx)`` tuples to scan through any given
        string. Each k-mer is translated to an integer in the following way: let
        the alphabet (which is accessed through :attr:`tuplesdb`) has length
        :math:`L` and the letters in the alphabet are :math:`l_0,\ldots,l_{L-1}`.
        Letter :math:`l_i` is replaced by the digit :math:`i` and the resulting
        sequence of digits is interpreted in base :math:`L` and reported in its
        equivalent decimal (base 10) representation. Note that:

            * This conversion reduces the required disk space by roughly a third
              and allows for more efficient searching and indexing of the tuples
              table.
            * Having a letter representing zero is OK as long as all represented
              words have the same length which is true here (otherwise ``ACCT``
              would have the same representation as ``CCT``).

        For example::

            string = 'ACGTGT'
            tup_scan(string, 5) # => (110, 0), (443, 1)

        where the integers 110 and 443 are decimal representations for 5-mers
        ``ACGTG`` (which is ``01232`` in base 4) and ``CGTGT`` (which is
        ``12323`` in base 4).

        Args:
            string (str): The string to scan.
            wordlen(int): Length of the words.

        Yields:
            tuple: A string of length :attr:`wordlen` and a starting position.
        """
        digits = {let:idx for idx,let in enumerate(self.tuplesdb.alphabet.letters)}
        for idx in range(len(string) - self.wordlen + 1):
            tup = string[idx:idx + self.wordlen]
            tup = sum(digits[x]*(4**i) for x,i in zip(tup,reversed(range(len(tup)))))
            yield (tup, idx)

    # Helper for index(): yields data values to be inserted in the seeds index.
    def _give_seeds(self, cursor):
        # count the total number of tuples so we can report percentage
        # progress:
        cursor.execute('SELECT count(*) FROM %s' % self.tuples_table)
        for row in cursor:
            num_tuples = int(row[0])
            break
        indicator = ProgressIndicator(
            'Indexing %d observed %d-mers' % (num_tuples, self.wordlen), num_tuples
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
            for S_hit_idx in range(len(hits)):
                for T_hit_idx in range(S_hit_idx + 1, len(hits)):
                    S_hit, T_hit = hits[S_hit_idx], hits[T_hit_idx]
                    # not interested in seeds from a sequence to itself:
                    if S_hit[0] == T_hit[0]:
                        continue

                    yield (S_hit[0], T_hit[0], S_hit[1], T_hit[1])

        indicator.finish()
        sys.stderr.write('Creating SQL index on %s'  % self.seeds_table)
        cursor.execute("""
            CREATE INDEX seeds_ids ON %s (S_id, T_id)
        """ % self.seeds_table)
        sys.stderr.write('.\n')

    def num_potential_homolog_pairs(self, cursor=None):
        """Returns the number of potential homolog pairs in the database.
        Provide a specific cursor object to access the SQLite database if
        accessing this function through a transaction.
        """
        cnt_q = 'SELECT COUNT(*) FROM (%s)' % self.potential_homologs_q
        if cursor is None:
            with sqlite3.connect(self.tuplesdb.db) as conn:
                cursor = conn.cursor()
                cursor.execute(cnt_q)
        else:
            cursor.execute(cnt_q)
        for row in cursor:
            return row[0]

    # Helper for index(): yields data values for the potential homologs index.
    # Note that for each sequence id (which is an integer as per the seq table),
    # only potential homologs with greater sequence IDs are reported to avoid
    # duplicates.
    def _give_potential_homologs(self, cursor):
        num_total = self.num_potential_homolog_pairs(cursor)
        msg = 'Indexing potentially homologous pairs of sequences with at least %d seeds' % self.min_seeds_for_homology
        indicator = ProgressIndicator(msg, num_total)
        indicator.start()
        cursor.execute(self.potential_homologs_q)
        for row in cursor:
            if int(row[2]) < self.min_seeds_for_homology:
                continue;
            if row[0] < row[1]:
                yield (row[0], row[0], str(row[1]) + ',')
            elif row[0] > row[1]:
                yield (row[1], row[1], str(row[0]) + ',')
            else:
                raise RuntimeError("This shouldn't have happened, row=%s", str(row))
            indicator.progress()

        indicator.finish()

    def index(self):
        """Scans all sequences in the ``seq`` table and records all observed
        tuples in ``tuples_N`` table. Then scans all the observed hits to
        populate the ``seeds_N`` table used by :func:`seeds`. Each "seed"
        record is a 4-tuple ``(S_id, T_id, S_idx, T_idx)`` where the first two
        entries are ID's of sequences (as per the ``seq`` table) and the last
        two entries are integers corresponding to the starting position of the
        seed in the sequence pair.

        Finally, the contents of the ``seeds_N`` table is scanned to populate
        ``potential_homologs_N`` for every sequence in the database.
        """
        hit_ins_q = """
            INSERT OR REPLACE INTO %s (tuple, hits)
            SELECT ?, IFNULL( (SELECT hits FROM %s WHERE tuple = ?), "") || ?
        """ % (self.tuples_table, self.tuples_table)
        with sqlite3.connect(self.tuplesdb.db) as conn:
            # We need multiple cursors: one two read from the seq table
            # (seq_c), one two read/write from/to the tuples table (tuples_c),
            # one to read/write from/to the seeds table (seeds_c), and one
            # to write to the seeds cache table (potential_homologs_c)
            tuples_c = conn.cursor()
            seq_c = conn.cursor()
            seeds_c = conn.cursor()
            potential_homologs_c = conn.cursor()

            # count the total number so we can report percentage progress:
            seq_c.execute('SELECT count(*) FROM seq')
            for row in seq_c:
                num_seqs = int(row[0])
                break

            msg = 'Scanning %d sequences for %d-mers' \
                % (num_seqs, self.wordlen)
            indicator = ProgressIndicator(msg, num_seqs)
            indicator.start()

            # populate the tuples table:
            def _give_tuple(string):
                for s, idx in self.tup_scan(string):
                    yield (s, s, '@%s:%d' % (seqid, idx))

            seq_c.execute('SELECT id, seq FROM seq')
            for seqid, seqstr in seq_c:
                tuples_c.executemany(hit_ins_q, _give_tuple(seqstr))
                indicator.progress()

            indicator.finish()

            # populate the seeds table:
            seed_ins_q = """
                INSERT INTO %s (S_id, T_id, S_idx, T_idx)
                VALUES (?, ?, ?, ?)
            """ % self.seeds_table
            seeds_c.executemany(seed_ins_q, self._give_seeds(tuples_c))

            # populate the seeds cache table:
            seed_cache_ins_q = """
                INSERT OR REPLACE INTO %s (id, homologs)
                SELECT ?, IFNULL( (SELECT homologs FROM %s WHERE id = ?), "") || ?
            """ % (self.potential_homologs_table, self.potential_homologs_table)
            potential_homologs_c.executemany(seed_cache_ins_q, self._give_potential_homologs(seeds_c))

    def verify(self):
        """Verifies the database contents to make sure:

            * There seeds table contains no more than one for each 4-tuple
              ``(S_id,T_id,S_idx,T_idx)``.
            * The potential homologs table contains no more than one row for
              each sequence in ``seq``.

            If an inconsistency is found an ``AssertionError`` is raised.

            Note:
                The point of this function is to avoid setting up the above
                constraints on the SQLite tables. The reason is:
                * SQLite does not allow modifying constraints on a table after
                  creation.
                * Having the constraints in place while we are munging seeds
                  cripples performance.
        """
        with sqlite3.connect(self.tuplesdb.db) as conn:
            c = conn.cursor()
            sys.stderr.write('Verifying consistency of database:')
            # potential homologs should have at most a single row per sequence.
            potential_homologs_c.execute("""
                SELECT COUNT(*) FROM %s
            """ % self.potential_homologs_table)
            for row in potential_homologs_c:
                num_records = row[0]
            potential_homologs_c.execute("""
                SELECT COUNT(*) FROM (SELECT DISTINCT id FROM %s)
            """ % self.potential_homologs_table)
            for row in potential_homologs_c:
                assert(num_records == row[0])

            seeds_c.execute("""
                SELECT COUNT(*) FROM %s
            """ % self.seeds_table)
            for row in seeds_c:
                num_records = row[0]
            potential_homologs_c.execute("""
                SELECT COUNT(*) FROM (SELECT DISTINCT S_id,T_id,S_idx,T_idx FROM %s)
            """ % self.seeds_table)
            for row in seeds_c:
                assert(num_records == row[0])

            sys.stderr.write('looks good!\n')

    def potential_homologs(self, seqid):
        """Given a sequence ID, returns a list of sequence IDs, all integers
        greater than ``seqid`` (to avoid duplicates in assembly), that are
        potentially homologous to the given sequence. The result is read
        directly from the ``potential_homologs_N`` table.
        """
        q = """
            SELECT homologs FROM %s WHERE id = ?
        """ % self.potential_homologs_table
        with sqlite3.connect(self.tuplesdb.db) as conn:
            c = conn.cursor()
            c.execute(q, (seqid,))
            for row in c:
                return [int(i) for i in row[0].split(',') if i]

        return []

    def seeds(self, S_id, T_id):
        """Given two sequence ids, finds all exactly matching segments
        (see :class:`Segment`) of length :attr:`wordlen` between the two.
        Segments are not necessarily in maximal form. For purposes of seed
        extension, however, we prefer to not have too many segments that are
        part of a one bigger segments (especially if they do not belong to an
        actual overlap alignment). This can be worked out by using :func:`maximal_seeds`
        which reduces any set of seeds into maximal, necessarily non-overlapping
        segments.

        Args:
            S_id (int): The database ID of the "from" sequence.
            T_id (int): The database ID of the "to" sequence.

        Note:
            Scores of transcripts for exact matches are left as 0 since this
            class does not concern itself with alignment scores.
        """
        q = """
            SELECT S_id, T_id, S_idx, T_idx FROM %s
            WHERE (S_id = ? AND T_id = ? ) OR (S_id = ? AND T_id = ?)
        """ % (self.seeds_table)
        seeds = []
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
                tx = pw.Transcript(
                    S_idx=S_idx, T_idx=T_idx, score=0, opseq='M'*self.wordlen
                )
                seeds += [Segment(S_id=S_id, T_id=T_id, tx=tx)]

        return seeds

    @classmethod
    def maximal_seeds(cls, seeds, S_id, T_id):
        """Given a list of exactly matching segments, reduces them into a list
        of *maximal* exactly matching segments in increasing order of ``S_idx``.
        A segment is in its maximal form if it cannot be extended in either
        direction by an exact match.

        Args:
            list[Segment]: Exactly matching segments, potentially overlapping.
            S_id (int): The database ID of the "from" sequence.
            T_id (int): The database ID of the "to" sequence.

        Returns
            list[Segment]: Maximal segments, guaranteed to not overlap.
        """
        seeds.sort(key=lambda s: s.tx.S_idx)
        # merge overlapping tuples:
        idx = 0
        while idx < len(seeds):
            cand = idx + 1
            while cand < len(seeds):
                shift_S = seeds[cand].tx.S_idx - seeds[idx].tx.S_idx
                shift_T = seeds[cand].tx.T_idx - seeds[idx].tx.T_idx
                # we know the transcripts are all M's.
                if shift_S == shift_T and shift_S > 0 and \
                   shift_S < lib.strlen(seeds[idx].tx.opseq):
                    tx = pw.Transcript(
                        S_idx=seeds[idx].tx.S_idx,
                        T_idx=seeds[idx].tx.T_idx,
                        score=0,
                        opseq=seeds[cand].tx.opseq + 'M'*shift_S
                    )
                    seeds[idx] = Segment(S_id=S_id, T_id=T_id, tx=tx)
                    seeds.pop(cand)
                else:
                    cand += 1
            idx += 1
        return seeds
