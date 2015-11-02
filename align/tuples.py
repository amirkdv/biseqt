#!/usr/bin/env python
import sys
import re
import sqlite3
from Bio import SeqIO
from collections import namedtuple

from . import align, seq

# a pair of equal-length substrings of two sequences
Segment = namedtuple('Segment', ['id_S', 'id_T', 'idx_S', 'idx_T', 'len'])

class MaxConcurrentQueries(RuntimeError):
    pass

def tup_scan(string, wordlen):
    """A generator for (string, idx) tuples to scan through any given string.
    """
    for idx in range(len(string) - wordlen):
        yield (string[idx:idx + wordlen], idx)

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
        assert isinstance(alphabet, seq.Alphabet)
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
                  'description' text,
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

    # TODO allow specifying a homopolymeric.HpCondenser
    def populate(self, fasta_src=None, lim=-1):
        """Given a FASTA source file, loads all the sequences (up to a limit,
        if specified) into the `seq` table. No indexing is done; see index().
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

    def index(self):
        """Scans all sequences in the `seq` table and records all observed
        tuples in `tuples_N` tables and all hits in `tuples_N_hits` tables.

        NOTE takes ~3 minutes to index 10-mers of 500 sequences with average
        length 11Kbp into a DB of ~250MB. With 20-mers takes ~4 minutes and DB
        size is ~600MB!!
        """
        sys.stderr.write('indexing sequences:')
        def give_hit(seqid, string):
            for s,idx in tup_scan(string, self.wordlen):
                yield (seqid, idx, s)

        def give_tup(string):
            for s,idx in tup_scan(string, self.wordlen):
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
                return seq.Sequence(str(row[0]), self.alphabet)

    def seqids(self, info=False):
        """Returns a list of all seqids.
        """
        with sqlite3.connect(self.db) as conn:
            c = conn.cursor()
            if info:
                c.execute("SELECT id, description, LENGTH(seq) FROM seq")
                seqs = {}
                for row in c:
                    seqs[row[0]] = {}
                    seqs[row[0]]['description'] = row[1]
                    seqs[row[0]]['length'] = row[2]
                    # see seq.make_sequencing_fixture()
                    seqs[row[0]]['start'] = int(re.split('\(|\)', row[1])[1].split(':')[1])
                return seqs
            else:
                c.execute('SELECT id FROM seq')
                return [row[0] for row in c]

    def exactly_matching_segments(self, s1, s2):
        """Given two seqids, finds all maximal exactly matching segments between
        the two self.seqid) and another given sequence. A segment is in its
        maximal form if there are no other exactly-matching segments whose span
        is a unit shift in both the query and target sequences.
        """
        segments = []
        q = """
            SELECT H1.idx, H2.idx FROM tuples_{}_hits H1 -- H1 is s1, H2 is s2
            INNER JOIN tuples_{}_hits H2
            ON H1.tuple = H2.tuple
            WHERE H1.seq = ? AND H2.seq = ?
        """.format(self.wordlen, self.wordlen)
        def row_factory(cursor, row):
            return Segment(id_S=s1, id_T=s2, len=self.wordlen,
                idx_S=row[0], idx_T=row[1])

        with sqlite3.connect(self.db) as conn:
            conn.row_factory = row_factory
            c = conn.cursor()
            c.execute(q, (s1,s2))
            segments = [r for r in c]

        segments.sort(key=lambda s: s.idx_S)
        # merge overlapping tuples:
        idx = 0
        while idx < len(segments):
            cand = idx + 1
            while cand < len(segments):
                shift_S = segments[cand].idx_S - segments[idx].idx_S
                shift_T = segments[cand].idx_T - segments[idx].idx_T
                if shift_S == shift_T and shift_S > 0 and shift_S < segments[idx].len:
                    segments[idx] = Segment(
                        id_S=s1, id_T=s2, len=segments[cand].len + shift_S,
                        idx_S=segments[idx].idx_S, idx_T=segments[idx].idx_T,
                    )
                    segments.pop(cand)
                else:
                    cand += 1
            idx += 1
        segments.sort(key=lambda s: s.len, reverse=True)
        return segments

class OverlapFinder(object):
    def __init__(self, S, T, align_params):
        self.S, self.T, self.align_params = S, T, align_params
        self.P = align.AlignProblem(
            S=S, T=T, params=align_params,
            align_type=align.ALIGN_GLOBAL,
            S_min_idx=0, S_max_idx=0, T_min_idx=0, T_max_idx=0
        )

    def extend_one_way_once(self, segment, window, direction='fwd'):
        """Helper method for extend_one_way."""
        # avoid overflows
        if direction == 'fwd':
            window = min(window, min(
                self.S.length - segment.idx_S - segment.len,
                self.T.length - segment.idx_T - segment.len
                )
            )
            self.P.S_min_idx = segment.idx_S + segment.len
            self.P.T_min_idx = segment.idx_T + segment.len
            self.P.S_max_idx = self.P.S_min_idx + window
            self.P.T_max_idx = self.P.T_min_idx + window
        elif direction == 'bwd':
            window = min(window, min(segment.idx_S, segment.idx_T))
            self.P.S_max_idx = segment.idx_S
            self.P.T_max_idx = segment.idx_T
            self.P.S_min_idx = self.P.S_max_idx - window
            self.P.T_min_idx = self.P.T_max_idx - window

        if window == 0:
            return None, None

        score = self.P.solve()
        if direction == 'fwd':
            segment = Segment(
                len=segment.len + window,
                id_S=segment.id_S, id_T=segment.id_T,
                idx_S=segment.idx_S, idx_T=segment.idx_T,
            )
        elif direction == 'bwd':
            segment = Segment(
                len=segment.len + window,
                id_S=segment.id_S, id_T=segment.id_T,
                idx_S=segment.idx_S - window, idx_T=segment.idx_T - window,
            )
        return segment, score

    def extend_one_way(self, segment, drop_threshold, direction='fwd'):
        """Extends a given segment (presumably maximally-exactly-matching) in
        the given direction (fwd or bwd). Extension is done by repeatedly
        performing global alignments on a rolling window along the two sequences
        and stopping once a decrease in score is observed for a certain number
        of times"""
        assert(direction in ['fwd', 'bwd'])
        window = segment.len
        cur_seg, cur_score = segment, 0
        while cur_score >= drop_threshold:
            offset = 0
            seg, score  = self.extend_one_way_once(cur_seg, window, direction=direction)
            if seg is None:
                # hit the end:
                return cur_seg, cur_score
            cur_score += score
            cur_seg = seg

        return None, None

    def extend(self, segments, drop_threshold):
        """Given a number of matching segments for the two sequences finds all
        extended matching gap containing segments by repeatedly aligning a
        rolling frame along the two sequences and dropping a segment once it
        is observed to not be a high-scoring overlap alignment"""
        scores_by_segment = {}
        for segment in segments:
            self.P.S_min_idx, self.P.T_min_idx = segment.idx_S, segment.idx_T
            core_score = self.P.score('B' + 'M'*segment.len)
            window = len(segment)
            fwd_extended, fwd_score = self.extend_one_way(segment, drop_threshold, 'fwd')
            bwd_extended, bwd_score = self.extend_one_way(segment, drop_threshold, 'bwd')
            if fwd_extended and bwd_extended and fwd_score + bwd_score > drop_threshold:
                segment = Segment(
                    id_S=segment.id_S, id_T=segment.id_T,
                    idx_S=bwd_extended.idx_S,
                    idx_T=bwd_extended.idx_T,
                    len=bwd_extended.len + fwd_extended.len - segment.len
                )
                assert(segment.idx_S == 0 or segment.idx_T == 0)
                scores_by_segment[segment] = fwd_score + bwd_score + core_score
        return sorted(scores_by_segment.items(), key=lambda x: x[1], reverse=True)
