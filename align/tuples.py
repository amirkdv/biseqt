"""Provides tools for *k*-mer analysis."""

import sys
import re
import sqlite3
from Bio import SeqIO
from collections import namedtuple

from . import pw, seq
class Segment(namedtuple('Segment', ['S_id', 'T_id', 'tx'])):
    """Represents an aligned pair of substrings in two sequences. The alignment may
    potentially contain indels. Maximal, exactly-matching segments are refered to as
    "seeds".

    Attributes:
        S_id (int): The id of the "from" sequence as found in the ``seq`` table.
        T_id (int): The id of the "to" sequence as found in the ``seq`` table.
        tx (pw.Transcript):  object representing the alignment.
    """

def tup_scan(string, wordlen):
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
    for idx in range(len(string) - wordlen + 1):
        yield (string[idx:idx + wordlen], idx)

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
    def __init__(self, db, wordlen=10, alphabet=None):
        assert isinstance(alphabet, seq.Alphabet)
        self.alphabet = alphabet
        self.db, self.wordlen = db, wordlen
        assert isinstance(self.wordlen, int)

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
            q = """
                CREATE TABLE tuples_{} (
                  'id' integer PRIMARY KEY ASC,
                  'tuple' char({}),
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

    def populate(self, fasta_src, lim=-1):
        """Given a FASTA source file, loads all the sequences (up to a limit, if
        specified) into the `seq` table. No indexing is done; see ``index()``.

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

    def index(self):
        """Scans all sequences in the ``seq`` table and records all observed
        tuples in ``tuples_N`` tables and all hits in ``tuples_N_hits`` tables.

        Note:
            It takes ~3 minutes to index 10-mers of 500 sequences with average
            length 11Kbp into a DB of ~250MB. With 20-mers takes ~4 minutes and
            DB size is ~600MB!
        """
        # We want an Upsert to avoid changing tuple
        # IDs see http://stackoverflow.com/a/4330694 .
        tup_insert_q = """
            INSERT OR REPLACE INTO tuples_{} (id, tuple)
            VALUES ((SELECT id FROM tuples_{} WHERE tuple = ?), ?)
        """.format(self.wordlen, self.wordlen)
        hit_insert_q = """
            INSERT INTO tuples_{}_hits (tuple, seq, idx)
                SELECT id, ?, ?  FROM tuples_{} WHERE tuple = ?
        """.format(self.wordlen, self.wordlen)
        sys.stderr.write('indexing sequences:')
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
                c.executemany(
                    tup_insert_q,
                    [(s,s) for s, _ in tup_scan(string, self.wordlen)]
                )
                c.executemany(
                    hit_insert_q,
                    [(seqid, idx, s) for s, idx in tup_scan(string, self.wordlen)]
                )
                sys.stderr.write(' ' + str(seqid))
            sys.stderr.write('\n')

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
            T = tuples.TuplesDB('genome.db', alphabet=A, wordlen=5)
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

class OverlapFinder(object):
    """Provided two sequences and a :class:`TuplesDB` provides tools to find
    and extend exactly matching "seeds" into overlap alignments. Both sequences
    must have already been indexed in the tuples database.

    Seeds are *maximal, exactly matching* instances of :class:`Segment`. They
    are expanded by repeated global alignments on small windows moving forward
    and backward from the boundaries of the seed.

    Args:
        S (seq.Sequence): The "from" sequence.
        T (seq.Sequence): The "to" sequence.
        align_params (pw.AlignParams): The alignment parameters.

    Keyword Args:
        tuplesdb (tuples.TuplesDB)
    """
    def __init__(self, S, T, align_params, tuplesdb=None):
        self.tuplesdb = tuplesdb
        self.S, self.T, self.align_params = S, T, align_params
        self.align_problem_kw = {
            'S': self.S,
            'T': self.T,
            'params': self.align_params,
            'align_type': pw.GLOBAL,
        }

    def exactly_matching_segments(self, s1, s2):
        """Given two sequence ids, finds all maximal exactly matching segments
        between the two. A segment is in its maximal form if there are no other
        exactly-matching segments whose span is a unit shift in both the query
        and target sequences.

        Note:
            Scores of transcripts for exact matches are left as 0 to avoid
            unnecessary cycles, we populate the score first thing in
            :func:`extend`.
        """
        segments = []
        q = """
            SELECT H1.idx, H2.idx FROM tuples_{}_hits H1 -- H1 is s1, H2 is s2
            INNER JOIN tuples_{}_hits H2
            ON H1.tuple = H2.tuple
            WHERE H1.seq = ? AND H2.seq = ?
        """.format(self.tuplesdb.wordlen, self.tuplesdb.wordlen)
        def row_factory(cursor, row):
            tx = pw.Transcript(row[0], row[1], 0, 'M'*self.tuplesdb.wordlen)
            return Segment(S_id=s1, T_id=s2, tx=tx)

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
                    tx = pw.Transcript(
                        exacts[idx].tx.idx_S, exacts[idx].tx.idx_T,
                        0, # score
                        exacts[cand].tx.opseq + 'M'*shift_S # transcript
                    )
                    exacts[idx] = Segment(S_id=s1, T_id=s2, tx=tx)
                    exacts.pop(cand)
                else:
                    cand += 1
            idx += 1
        exacts.sort(key=lambda s: len(s.tx.opseq), reverse=True)
        return exacts

    def _extend_fwd_once(self, segment, window):
        """Helper method for ``extend1d``."""
        S_len, T_len = self._S_len(segment.tx.opseq), self._T_len(segment.tx.opseq)
        kw = self.align_problem_kw
        kw.update({
            'S_min_idx': segment.tx.idx_S + S_len,
            'T_min_idx': segment.tx.idx_T + T_len,
        })
        kw.update({
            'S_max_idx': kw['S_min_idx'] + window,
            'T_max_idx': kw['T_min_idx'] + window
        })

        with pw.AlignProblem(**kw) as P:
            score = P.solve()
            assert(score is not None)
            transcript = P.traceback()
            if transcript is None:
                return None

        tx = pw.Transcript(
            idx_S=segment.tx.idx_S,
            idx_T=segment.tx.idx_T,
            score=segment.tx.score + transcript.score,
            opseq=segment.tx.opseq + transcript.opseq
        )
        return Segment(S_id=segment.S_id, T_id=segment.T_id, tx=tx)

    def _extend_bwd_once(self, segment, window):
        """Helper method for ``extend1d``."""
        kw = self.align_problem_kw
        kw.update({
            'S_max_idx': segment.tx.idx_S,
            'T_max_idx': segment.tx.idx_T,
        })
        kw.update({
            'S_min_idx': kw['S_max_idx'] - window,
            'T_min_idx': kw['T_max_idx'] - window
        })

        with pw.AlignProblem(**kw) as P:
            score = P.solve()
            assert(score is not None)
            transcript = P.traceback()
            if transcript is None:
                return None

        tx = pw.Transcript(
            idx_S=transcript.idx_S,
            idx_T=transcript.idx_T,
            score=transcript.score + segment.tx.score,
            opseq=transcript.opseq + segment.tx.opseq
        )
        return Segment(S_id=segment.S_id, T_id=segment.T_id, tx=tx)

    def _S_len(self, opseq):
        return sum([opseq.count(op) for op in 'DMS'])

    def _T_len(self, opseq):
        return sum([opseq.count(op) for op in 'IMS'])

    def _extend1d(self, segment, backwards=False, **kw):
        """Helper method for ``extend()``"""
        drop_threshold = kw['drop_threshold']
        window = kw['window']
        max_succ_drops = kw['max_succ_drops']
        cur_seg = segment
        score_history = [segment.tx.score]
        while True:
            if backwards:
                w = min(window, min(cur_seg.tx.idx_S, cur_seg.tx.idx_T))
            else:
                S_wiggle = self.S.length - (cur_seg.tx.idx_S + self._S_len(cur_seg.tx.opseq))
                T_wiggle = self.T.length - (cur_seg.tx.idx_T + self._T_len(cur_seg.tx.opseq))
                w = min(window, min(S_wiggle, T_wiggle))

            if w == 0:
                # hit the end:
                return cur_seg

            if backwards:
                seg = self._extend_bwd_once(cur_seg, w)
            else:
                seg = self._extend_fwd_once(cur_seg, w)

            if seg is None:
                # no non-empty alignment found.
                break

            score_history += [seg.tx.score - segment.tx.score]
            if len(score_history) > max_succ_drops:
                score_history = score_history[-max_succ_drops:]

            if len(score_history) == max_succ_drops and \
                all([x <= drop_threshold for x in score_history]):
                break

            cur_seg = seg

        return None

    def extend(self, segments, **kw):
        """Given a number of matching segments for the two sequences finds all
        extended matching gap containing segments by repeatedly aligning a
        rolling frame along the two sequences and dropping a segment once it
        is observed to not be a high-scoring overlap alignment.

        Args:
            segment (tuples.Segment)

        Keyword Args:
            drop_threshold (Optional[float]): What constitutes a drop in the
                score from one window to the next, default is 0.
            window (Optional[int]): The size of the rolling window.
            max_succ_drops (Optional[int]): Maximum number of "drops" until the
                segment is dropped (i.e ``None`` is returned).

        Note:
            It seems like we never have two segments for the same pair of
            sequences where one gives an ``S -> T`` edge and the other gives a
            ``T -> S`` edge. Why not just return the first segment that goes all
            the way to the end?
        """
        defaults = {'max_succ_drops': 3, 'drop_threshold': 0, 'window': 10}
        defaults.update(kw)
        kw = defaults
        res = []
        for segment in segments:
            fwd = self._extend1d(segment, **kw)
            bwd = self._extend1d(segment, backwards=True, **kw)
            if fwd and bwd and min(fwd.tx.score, bwd.tx.score) > kw['drop_threshold']:
                assert(bwd.tx.idx_S == 0 or bwd.tx.idx_T == 0)
                opseq = bwd.tx.opseq[:-len(segment.tx.opseq)] + fwd.tx.opseq
                score = self.align_params.score(
                    self.S, self.T, opseq,
                    S_min_idx=bwd.tx.idx_S,
                    T_min_idx=bwd.tx.idx_T
                )
                tx = pw.Transcript(
                    idx_S=bwd.tx.idx_S,
                    idx_T=bwd.tx.idx_T,
                    score=score,
                    opseq=opseq
                )
                res += [Segment(S_id=segment.S_id, T_id=segment.T_id, tx=tx)]
        return sorted(res, key=lambda s: s.tx.score, reverse=True)
