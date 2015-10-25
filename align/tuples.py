#!/usr/bin/env python
import sys
import sqlite3
from Bio import SeqIO
from collections import namedtuple

from . import utils
from . import align
from . import scan

Seed = namedtuple('Seed', ['seqid', 'idx', 'idx_q', 'len'])

class MaxConcurrentQueries(RuntimeError):
    pass

class TuplesDB(object):
    """Wraps an Sqlite database containing tuple indices for sequences. For now,
    only one word length is allowed. The tables are:
        * a `seq` table containing all the provided sequences (the "database"),
          each sequence and its name are mapped to an internal integer used
          for calculations (see populate()).
        * a `tuples_N` table for each word length N. This table has two
          purposes:
          1. avoid repeating full tuple strings in records representing hits.
          2. delegate some of the hit-search calculation to SQLite (see Query).
        * a `tuples_N_hits` table for each word length N. This table contains
          tuple hit information about sequences in the "database" (see index()).

    Attributes:
        db(string):             path to the SQLite datbase
        wordlen(int):           length of tuples
        tup_insert_q(string):   an SQL query template that emulates an UPSERT
            operation in the SQL dialect of SQLite; use it to bulk-write tuples
            into `tuples_N` tables.
        hit_insert_q(string):   an SQL query template that records a tuple-hit;
            use it to bulk-write tuple-hits to `tuples_N_hits` tables.

    """
    def __init__(self, db=None, wordlen=10, alphabet=None):
        assert isinstance(alphabet, align.Alphabet)
        self.alphabet = alphabet
        self.db, self.wordlen = db, wordlen
        assert isinstance(self.wordlen, int)
        # we want an Upsert to avoid changing tuple
        # IDs see http://stackoverflow.com/a/4330694 .
        self.tup_insert_q = """
            INSERT OR REPLACE INTO tuples_{} (id, tuple, query)
            VALUES ((SELECT id FROM tuples_{} WHERE tuple = ?), ?, ?)
        """.format(self.wordlen, self.wordlen)
        self.hit_insert_q = """
            INSERT INTO tuples_{}_hits (tuple, seq, idx)
                SELECT id, ?, ?  FROM tuples_{} WHERE tuple = ?
        """.format(self.wordlen, self.wordlen)

    def initdb(self):
        """Initializes the database: creates the required tables and
        fails if any of them already exist"""
        with sqlite3.connect(self.db) as conn:
            c = conn.cursor()
            q = """
                CREATE TABLE seq (
                  'id' integer PRIMARY KEY ASC,
                  'name' text,
                  'seq' text
                );
            """
            c.execute(q)
            q = """
                CREATE TABLE tuples_{} (
                  'id' integer PRIMARY KEY ASC,
                  'tuple' char({}),
                  'query' integer DEFAULT -1,
                  UNIQUE(tuple)
                );
            """.format(self.wordlen, self.wordlen)
            c.execute(q)
            q = """
                CREATE TABLE tuples_{}_hits (
                  'tuple' integer REFERENCES tuples_{}(id),
                  'seq'   integer REFERENCES seq(id),
                  'idx'   integer,
                  UNIQUE(seq, idx)
                );
            """.format(self.wordlen, self.wordlen)
            c.execute(q)
        sys.stderr.write('initialized tuples DB at: %s\n' % self.db)

    def populate(self, fasta_src=None, lim=False):
        """Given a FASTA source file, loads all the sequences (up to a limit,
        if specified) into the `seq` table. No indexing is done; see index().
        """
        def give_seq():
            for idx, seq in enumerate(SeqIO.parse(fasta_src, "fasta")):
                if lim < 0 or idx < lim:
                    yield (str(seq.id), str(seq.seq))
                else:
                    break

        with sqlite3.connect(self.db) as conn:
            c = conn.cursor()
            q = "INSERT INTO seq (name, seq) VALUES (?,?)"
            c.executemany(q, give_seq())

    def index(self):
        """Scans all sequences in the `seq` table and records all observed
        tuples in `tuples_N` tables and all hits in `tuples_N_hits` tables.

        NOTE takes ~3 minutes to index 10-mers of 500 sequences with average
        length 11Kbp into a DB of ~250MB. With 20-mers takes ~4 minutes and DB
        size is ~600MB!!
        """
        sys.stderr.write('indexing sequences:')
        def give_hit(seqid, string):
            for s,idx in scan(string, self.wordlen):
                yield (seqid, idx, s)

        def give_tup(string):
            for s,idx in scan(string, self.wordlen):
                yield (s, s, -1)

        with sqlite3.connect(self.db) as conn:
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
                c.executemany(self.tup_insert_q, give_tup(string))
                c.executemany(self.hit_insert_q, give_hit(seqid, string))
                sys.stderr.write(' ' + str(seqid))
            sys.stderr.write('\n')

    def loadseq(self, seqid):
        """Loads a sequence given its internal numeric seqid.
        """
        with sqlite3.connect(self.db) as conn:
            c = conn.cursor()
            c.execute('SELECT seq FROM seq WHERE id = ?', (seqid,))
            for row in c:
                return str(row[0])

class Query(object):
    """A helper class to deal with queries and seeds against a TuplesDB

    Attributes:
        tuplesdb (TuplesDB): the database to use.
        align_params (align.AlignParams): alignment parameters used to expand
            seeds.

    NOTE Due to the way queries are implemented (see hitsummary(), for example)
    only one query at a time can be using the same tuples database file.
    """
    def __init__(self, qseq, tuplesdb=None, align_params=None):
        self.tuplesdb, self.qseq, self.align_params = tuplesdb, qseq, align_params
        def give_tup():
            for s, idx in scan(self.qseq, self.tuplesdb.wordlen):
                yield (s, s, idx)

        # index the query string:
        with sqlite3.connect(tuplesdb.db) as conn:
            c = conn.cursor()
            q = 'UPDATE tuples_{} SET query = -1'.format(tuplesdb.wordlen)
            c.execute(q)
            c.executemany(tuplesdb.tup_insert_q, give_tup())

        self.S = align.Sequence(qseq, tuplesdb.alphabet)

    def seeds(self, seqid=None):
        """Finds all the seeds in their maximal form given a query string. A
        seed is in its maximal form if there are no other seeds whose span
        is a unit shift in both the query and target sequences.
        """
        seeds = []
        w = self.tuplesdb.wordlen
        q = """
            SELECT T.query, H.idx, H.seq FROM tuples_{}_hits H
            INNER JOIN tuples_{} T ON T.id = H.tuple
            WHERE T.query <> -1
            AND H.seq = ?
        """.format(w, w, w)
        row_factory = lambda c,r: Seed(seqid=r[2], len=w, idx_q=r[0], idx=r[1])
        with sqlite3.connect(self.tuplesdb.db) as conn:
            conn.row_factory = row_factory
            c = conn.cursor()
            c.execute(q, (seqid,))
            seeds = [r for r in c]

        seeds.sort(key=lambda k: k.idx)
        # merge overlapping tuples:
        for idx,s in enumerate(seeds):
            if idx == 0:
                continue
            if seeds[idx-1].idx + seeds[idx-1].len + 1 == s.idx + s.len and \
                seeds[idx-1].idx_q + seeds[idx-1].len + 1 == s.idx_q + s.len:
                seeds[idx-1] = Seed(
                    seqid=seeds[idx-1].seqid,
                    idx=seeds[idx-1].idx,
                    idx_q=seeds[idx-1].idx_q,
                    len=seeds[idx-1].len + 1
                )
                seeds.pop(idx)
        return seeds

    def hitsummary(self):
        """Returns a dict of internal numeric seqids mapped to the number of
        hits they generate from the given query string.
        """
        with sqlite3.connect(self.tuplesdb.db) as conn:
            c = conn.cursor()
            q = """
                SELECT H.seq, count(H.seq) FROM tuples_{}_hits H
                INNER JOIN tuples_{} T
                ON T.id = H.tuple
                WHERE T.query <> -1
                GROUP BY H.seq
                ORDER BY count(H.seq) DESC
            """.format(self.tuplesdb.wordlen, self.tuplesdb.wordlen)
            c.execute(q)
            res = {row[0]:row[1] for row in c}
            return res

    def expand_seed(self, seed, window=20):
        """Expands a seed by repeatedly aligning portions (of length
        `window`) of the query and target sequences until a threshold low score
        is met.
        """
        T = align.Sequence(self.tuplesdb.loadseq(seed.seqid), self.tuplesdb.alphabet)
        S_min_idx, S_max_idx = seed.idx_q + seed.len, seed.idx_q + seed.len + window,
        T_min_idx, T_max_idx = seed.idx   + seed.len, seed.idx   + seed.len + window
        P = align.AlignProblem(
            S=self.S, T=T, params=self.align_params,
            align_type=align.ALIGN_START_ANCHORED,
            S_min_idx=S_min_idx, S_max_idx=S_max_idx,
            T_min_idx=T_min_idx, T_max_idx=T_max_idx
        )
        transcript = P.solve()
        # FIXME inspect transcript and keep extending if we are scoring well.
        if transcript[:3] != 'Err':
            infostr, transcript = transcript.split(':')
            transcript = infostr + ':B' + 'M' * seed.len + transcript[1:]
            utils.print_alignment(self.S, T, transcript, sys.stderr, margin=10)
        return transcript

