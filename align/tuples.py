#!/usr/bin/env python
import sys
import re
import sqlite3
from Bio import SeqIO
from collections import namedtuple

from . import align, seq

# an aligned pair of substrings in two sequences, the alignments are represented
# by a Transcript object (tx) which may contain substitutions or gaps.
Segment = namedtuple('Segment', ['id_S', 'id_T', 'tx'])

class MaxConcurrentQueries(RuntimeError):
    pass

def tup_scan(string, wordlen):
    """A generator for (string, idx) tuples to scan through any given string.
    """
    for idx in range(len(string) - wordlen + 1):
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

class OverlapFinder(object):
    def __init__(self, S, T, align_params, tuplesdb=None):
        self.tuplesdb = tuplesdb
        self.S, self.T, self.align_params = S, T, align_params
        self.P = align.AlignProblem(
            S=S, T=T, params=align_params,
            align_type=align.ALIGN_START_ANCHORED,
            S_min_idx=0, S_max_idx=0, T_min_idx=0, T_max_idx=0
        )

    def exactly_matching_segments(self, s1, s2):
        """Given two seqids, finds all maximal exactly matching segments between
        the two self.seqid) and another given sequence. A segment is in its
        maximal form if there are no other exactly-matching segments whose span
        is a unit shift in both the query and target sequences.

        Note: Scores of transcripts for exact matches are left as 0 to avoid
        unnecessary cycles, we populate the score first thing in extend()
        """
        segments = []
        q = """
            SELECT H1.idx, H2.idx FROM tuples_{}_hits H1 -- H1 is s1, H2 is s2
            INNER JOIN tuples_{}_hits H2
            ON H1.tuple = H2.tuple
            WHERE H1.seq = ? AND H2.seq = ?
        """.format(self.tuplesdb.wordlen, self.tuplesdb.wordlen)
        def row_factory(cursor, row):
            tx = align.Transcript(row[0], row[1], 0, 'M'*self.tuplesdb.wordlen)
            return Segment(id_S=s1, id_T=s2, tx=tx)

        with sqlite3.connect(self.tuplesdb.db) as conn:
            conn.row_factory = row_factory
            c = conn.cursor()
            c.execute(q, (s1,s2))
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
                    tx = align.Transcript(
                        exacts[idx].tx.idx_S, exacts[idx].tx.idx_T,
                        0, # score
                        exacts[cand].tx.opseq + 'M'*shift_S # transcript
                    )
                    exacts[idx] = Segment(id_S=s1, id_T=s2, tx=tx)
                    exacts.pop(cand)
                else:
                    cand += 1
            idx += 1
        exacts.sort(key=lambda s: len(s.tx.opseq), reverse=True)
        return exacts


    def score_gap(self, L):
        if not L:
            return 0
        return self.P.params.gap_open_score + L*self.P.params.gap_extend_score

    def extend_fwd_once(self, segment, window):
        S_len, T_len = self._S_len(segment.tx), self._T_len(segment.tx)
        self.P.type = align.ALIGN_START_ANCHORED
        self.P.S_min_idx = segment.tx.idx_S + S_len
        self.P.T_min_idx = segment.tx.idx_T + T_len
        self.P.S_max_idx = self.P.S_min_idx + window
        self.P.T_max_idx = self.P.T_min_idx + window

        _, score = self.P.solve()
        transcript = self.P.traceback()
        if transcript is None:
            return None

        return Segment(id_S=segment.id_S, id_T=segment.id_T, tx=transcript)

    def extend_bwd_once(self, segment, window):
        self.P.type = align.ALIGN_END_ANCHORED
        self.P.S_max_idx = segment.tx.idx_S
        self.P.T_max_idx = segment.tx.idx_T
        self.P.S_min_idx = self.P.S_max_idx - window
        self.P.T_min_idx = self.P.T_max_idx - window

        _, score = self.P.solve()
        transcript = self.P.traceback()
        if transcript is None:
            return None

        return Segment(id_S=segment.id_S, id_T=segment.id_T, tx=transcript)

    def _S_len(self, transcript):
        return sum([transcript.opseq.count(op) for op in 'DMS'])

    def _T_len(self, transcript):
        return sum([transcript.opseq.count(op) for op in 'IMS'])

    def extend1d(self, segment, drop_threshold, max_succ_drops=3, backwards=False):
        """Extends a given segment (presumably maximally-exactly-matching) in
        the given direction (fwd or bwd). Extension is done by repeatedly
        performing global alignments on a rolling window along the two sequences
        and stopping once a decrease in score is observed for a certain number
        of times"""
        print segment
        window = min(self._S_len(segment.tx), self._T_len(segment.tx)) * 2
        cur_seg, cur_score = segment, 0
        score_history = [segment.tx.score]
        while True:
            if backwards:
                w = min(window, min(cur_seg.tx.idx_S, cur_seg.tx.idx_T))
            else:
                S_wiggle = self.S.length - (cur_seg.tx.idx_S + self._S_len(cur_seg.tx))
                T_wiggle = self.T.length - (cur_seg.tx.idx_T + self._T_len(cur_seg.tx))
                w = min(window, min(S_wiggle, T_wiggle))

            if w == 0:
                # hit the end:
                return cur_seg

            if backwards:
                seg = self.extend_bwd_once(cur_seg, w)
            else:
                seg = self.extend_fwd_once(cur_seg, w)

            if seg is None:
                break

            score_history += [seg.tx.score]
            if len(score_history) > max_succ_drops:
                score_history = score_history[-max_succ_drops:]
            if all([x <= drop_threshold for x in score_history]):
                break

            cur_seg = seg

        return None

    def extend(self, segments, drop_threshold):
        """Given a number of matching segments for the two sequences finds all
        extended matching gap containing segments by repeatedly aligning a
        rolling frame along the two sequences and dropping a segment once it
        is observed to not be a high-scoring overlap alignment"""
        res = []
        for segment in segments:
            self.P.S_min_idx, self.P.T_min_idx = segment.tx.idx_S, segment.tx.idx_T
            core_score = self.P.score(segment.tx.opseq)
            window = len(segment)
            fwd = self.extend1d(segment, drop_threshold)
            bwd = self.extend1d(segment, drop_threshold, backwards=True)
            if fwd and bwd and fwd.tx.score + bwd.tx.score > drop_threshold:
                tx = align.Transcript(
                    idx_S=bwd.tx.idx_S, idx_T=bwd.tx.idx_T,
                    score=core_score + bwd.tx.score + fwd.tx.score,
                    opseq=bwd.tx.opseq[:-len(segment.tx.opseq)] + segment.tx.opseq + fwd.tx.opseq
                )
                assert(bwd.tx.idx_S == 0 or bwd.tx.idx_T == 0)
                res += [Segment(id_S=segment.id_S, id_T=segment.id_T, tx=tx)]
        return sorted(res, key=lambda s: s.tx.score, reverse=True)
