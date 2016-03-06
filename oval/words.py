"""Provides tools for *k*-mer analysis."""

import sys
import os.path
import re
import sqlite3
from math import sqrt, erf, log
from Bio import SeqIO
from collections import namedtuple
from matplotlib import pyplot as plt
from . import pw, seq, ProgressIndicator, CffiObject, ffi, lib


class Index(object):
    """The main responsibility of an Index is to respond to :func:`seeds`
    with a list of :class:`oval.pw.Segment` s given two sequence IDs.
    This is done by recording all "words" observed in all sequences
    in the database during :func:`index`. Each index is uniquely defined by
    its :attr:`seqdb` and its :attr:`wordlen`. Each such index operates
    using 3 database tables with the following schema (replace `N`
    with :attr:`wordlen`). All tables are populated, in the order presented
    below, by :func:`index`.

    .. code-block:: sql

        CREATE TABLE words_N (
          'tuple' integer, -- decimal represetnation of N-tuple taken as a number
                           -- in base L (see tup_scan()).
          'hits'  varchar, -- '@' delimited string of hits with format (@<id>:<idx>)*
        )
        CREATE TABLE seeds_N (
          -- each record corresponds to a single exactly matching segment pair
          'S_id' integer REFERENCES seq(id),
          'T_id' integer REFERENCES seq(id),
          'S_idx' integer, -- starting position of seed in S
          'T_idx' integer, -- starting position of seed in T
        );
        CREATE TABLE scores_N (
          -- each record corresponds to a single shift for a read pair
          'S_id' integer REFERENCES seq(id),
          'T_id' integer REFERENCES seq(id),
          'shift' integer, -- shift of seed (S_idx - T_idx)
          'score' real, -- score of shift (log p-value)
          PRIMARY KEY (S_id, T_id, shift)
        );

    All attributes listed below are keyword arguments to the constructor.

    Attributes:
        seqdb (SeqDB): The sequence database.
        wordlen (int): Length of words, default is 10.
        min_seeds_for_homology (int): Minimum number of seeds between two
            sequences that makes them "potential homologs", default is 1.
        min_word_log_pvalue (float): Minimum allowed pvalue for an observed word
            for it to be considered as a seed; default is -5. This blocks words
            that appear "unusually often" (and therefore potentially belong to
            a repeat region), see :func:`index_seeds`.
    """
    def __init__(self, seqdb, **kwargs):
        self.seqdb = seqdb
        self.wordlen = kwargs.get('wordlen', 10)
        self.min_seeds_for_homology = kwargs.get('min_seeds_for_homology', 1)
        self.min_word_log_pvalue = kwargs.get('min_word_log_pvalue', -5)
        self.t_words = 'words_%d' % self.wordlen
        self.t_seeds = 'seeds_%d' % self.wordlen
        self.t_scores = 'scores_%d' % self.wordlen

    def initdb(self):
        """Initializes the database: creates the required tables and
        fails if any of them already exists."""
        with sqlite3.connect(self.seqdb.db) as conn:
            c = conn.cursor()
            q = """
                CREATE TABLE %s (
                  'tuple' integer primary key,
                  'hits'  varchar
                );
            """ % (self.t_words)
            c.execute(q)
            q = """
                CREATE TABLE %s (
                  'S_id' integer REFERENCES seq(id),
                  'T_id' integer REFERENCES seq(id),
                  'S_idx' integer,
                  'T_idx' integer
                );
            """ % (self.t_seeds)
            c.execute(q)
            q = """
                CREATE TABLE %s (
                  'S_id' integer REFERENCES seq(id),
                  'T_id' integer REFERENCES seq(id),
                  'shift' integer, -- shift of seed (S_idx - T_idx)
                  'score' real, --  = -log(p-value)
                  PRIMARY KEY (S_id, T_id, shift)
                );
            """ % (self.t_scores)
            c.execute(q)
        sys.stderr.write('Initialized index tables %s, %s, %s.\n'
            % (self.words_table, self.seeds_table, self.potential_homologs_table)
        )

    def tup_scan(self, string):
        """A generator for ``(word, idx)`` words to scan through any given
        string. Each k-mer is translated to an integer in the following way: let
        the alphabet (which is accessed through :attr:`seqdb`) has length
        :math:`L` and the letters in the alphabet are :math:`l_0,\ldots,l_{L-1}`.
        Letter :math:`l_i` is replaced by the digit :math:`i` and the resulting
        sequence of digits is interpreted in base :math:`L` and reported in its
        equivalent decimal (base 10) representation. Note that:

            * This conversion reduces the required disk space by roughly a third
              and allows for more efficient searching and indexing of the words
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
        digits = {let:idx for idx,let in enumerate(self.seqdb.alphabet.letters)}
        for idx in range(len(string) - self.wordlen + 1):
            tup = string[idx:idx + self.wordlen]
            tup = sum(digits[x]*(4**i) for x,i in zip(tup,reversed(range(len(tup)))))
            yield (tup, idx)

    def index(self):
        """Populates all index tables, in order: ``words_N``, ``seeds_N``, and
        ``potential_homologs_N`` by calling :func:`index_words`,
        :func:`index_seeds`, and :func:`index_potential_homologs`.
        """
        self.index_words()
        self.index_seeds()

    def index_words(self):
        """Scans all sequences in the ``seq`` table and records all observed
        k-mers (where k is :attr:`wordlen`). Each entry in the ``words_N``
        table records all occurences of a given word (which is represented as
        a decimal integer as per the output of :func:`tup_scan`). For example::

            tuple | hits
            9     | @1:5533@2:1436

        means the word represented by 9 (which, with word length 5, is
        ``AAAGC``) appears in sequence 1 (integer ID as per the ``seq``
        table) at position 5533 (starting at 0) and in sequence 2 at position
        1436.
        """
        hit_ins_q = """
            INSERT OR REPLACE INTO %s (tuple, hits)
            SELECT ?, IFNULL( (SELECT hits FROM %s WHERE tuple = ?), "") || ?
        """ % (self.words_table, self.words_table)
        with sqlite3.connect(self.seqdb.db) as conn:
            seq_c = conn.cursor() # to load sequences from the seq table.
            c = conn.cursor()     # to write words to the words table.

            # count the total number so we can report percentage progress:
            seq_c.execute('SELECT count(*) FROM seq')
            num_seqs = int(seq_c.next()[0])
            msg = 'Scanning %d sequences for %d-mers' \
                % (num_seqs, self.wordlen)
            indicator = ProgressIndicator(msg, num_seqs)
            indicator.start()

            # populate the words table:
            def _give_tuple(string):
                for s, idx in self.tup_scan(string):
                    yield (s, s, '@%s:%d' % (seqid, idx))

            seq_c.execute('SELECT id, seq FROM seq')
            for seqid, seqstr in seq_c:
                c.executemany(hit_ins_q, _give_tuple(seqstr))
                indicator.progress()

            indicator.finish()
            conn.commit()

    def index_seeds(self):
        """Finds all *non-repetitive* (described below) *seeds*
        between all pairs of sequences in ``seq``. A seed is an exactly matching
        segments of length equal to :attr:`wordlen`. This
        function assumes :func:`index_words` has been successfully called. For
        example given the following record in ``words_N``::

            tuple | hits
            9     | @1:5533@2:1436

        this function adds the following record to ``seeds_N``::

            S_id | T_id | S_idx | T_idx
            1    | 2    | 5533  | 1436

        Seeds are considered *repetitive* if they are suspected to the belong
        to a repeat region, calculated as follows.

        Given a database of sequences with total length (i.e sum of
        lengths) :math:`L`, a :math:`k`-mer occuring at :math:`n` places in the
        entire database is repetitive if the natural logarithm of the
        probability of :math:`n` occurences according to a random null
        hypothesis is below :attr:`min_word_log_pvalue` (i.e a high number of
        occurences implies a small p-value and potentially a repetitive
        structure). The random null hypothesis for :math:`X`, the
        random variable corresponding to the number of occurences of any
        :math:`k`-mer in a sequence of length :math:`L`, is known to be
        approximately binomial, that is:

            :math:`X \sim B(L, |\Sigma|^{-k})`

        where :math:`\Sigma` is the alphabet. For large :math:`L` the binomial
        can be approximated by a normal distribution:

            :math:`B(L, p_k) \\simeq \\mathcal{N}(Lp_k, \\sqrt{Lp_k(1-p_k)})`

        This function uses the above approximation to calculate p-values: a
        word occuring at :math:`n` positions in the database has p-value given
        by

            :math:`\\frac{N}{2}\\left[ 1 - \\mathrm{erf}\\left( \\frac{n-\mu}{\\sigma\\sqrt{2}} \\right) \\right]`

        where :math:`\\mu,\\sigma` are the parameters of the approximating
        normal distribution described above and :math:`N`, the total number of
        words, is a Bonferroni correction term.
        """
        with sqlite3.connect(self.seqdb.db) as conn:
            c = conn.cursor()
            q = """
                INSERT INTO %s (S_id, T_id, S_idx, T_idx)
                VALUES (?, ?, ?, ?)
            """ % self.seeds_table
            conn.cursor().executemany(q, self._give_seeds(c))

            sys.stderr.write('Creating SQL index on %s'  % self.seeds_table)
            c.execute("""
                CREATE INDEX seeds_ids ON %s (S_id, T_id)
            """ % self.seeds_table)
            sys.stderr.write('.\n')
            conn.commit()

    def word_pvalue_calculator(self, num_words):
        # Normal approximation (mean and standard deviation) of a binomial distribution
        B2N = lambda n, p: (n*p, sqrt(n * p * (1-p)))
        prob_word = lambda length: (1.0/self.seqdb.alphabet.length) ** length
        cursor.execute('select count(*) from %s' % self.words_table)
        cursor.execute('select sum(length(seq)) from seq')
        L = int(cursor.next()[0]) # total sequence length of all reads
        mu, sd = B2N(L, prob_word(self.wordlen)) # approximating normal distribution of word counts

        pvalue_calculator = lambda x: num_words * 0.5 * (1 - erf((x - mu) / (sd * sqrt(2))))
        return pvalue_calculator

    def band_radius(self, lengths, gap_prob, sensitivity=0.9):
        return 2*sqrt(gap_prob*min(lengths))*erfinv(sensitivity)

    def seed_pvalue_contribution(self, lengths, shift, gap_prob=0.1):
        r = band_radius(lengths, gap_prob)
        A = 2 * r * sqrt(sum(l - abs(shift)**2 for l in lengths))
        return r, sum(log(l) for l in lengths) - log(A), log(sum(lengths))
        FIXME

    # Helper for index(): yields data values to be inserted in the seeds index.
    def _give_seeds(self, cursor):
        word_pvalue_calc = self.word_pvalue_calculator(cursor)
        cursor.execute('select count(*) from %s' % self.words_table)
        N = int(cursor.next()[0]) # total number of words
        indicator = ProgressIndicator(
            'Indexing %d observed %d-mers' % (N, self.wordlen), N
        )
        indicator.start()
        cursor.execute('SELECT tuple, hits from %s' % self.words_table)
        for tup, hits in cursor:
            indicator.progress()
            pvalue = word_pvalue_calc(hits.count('@'))
            if not pvalue or log(pvalue) < self.min_word_log_pvalue:
                continue

            hits = [tuple(hit.split(':')) for hit in hits[1:].split('@')]
            hits = [(hit[0], int(hit[1])) for hit in hits]
            for S_hit_idx in range(len(hits)):
                for T_hit_idx in range(S_hit_idx + 1, len(hits)):
                    S_hit, T_hit = hits[S_hit_idx], hits[T_hit_idx]
                    # not interested in seeds from a sequence to itself:
                    if S_hit[0] == T_hit[0]:
                        continue

                    shift = S_hit[1] - T_hit[1]
                    lengths = [
                        self.seqinfo[S_hit[0]]['length'],
                        self.seqinfo[T_hit[0]]['length']
                    ]
                    r, s, s0 = seed_pvalue_contribution(lengths, shift)
                    FIXME # self.guesses = O(n^2) array of best observed shifts and their scores
                    q = """SELECT shift, score FROM scores
                        WHERE S_id = ? AND T_id = ? AND shift BETWEEN ? AND ?
                    """ FIXME correct index? currently PK(S_id, T_id, shift)
                    othercursor.execute(q, (S_hit[0], T_hit[0], shift + r, shift - r))
                    update_recs = (score + s if score > 0 else s0 for shift, score in othercursor)
                    q = """UPDATE scores SET score = ?
                        WHERE S_id = %d AND T_id = %d AND shift = %d
                    """ % (S_hit[0], T_hit[0], shift)
                    anothercursor.executemany(q, update_recs)

                    yield (S_hit[0], T_hit[0], S_hit[1], T_hit[1])

        indicator.finish()

    def num_potential_homolog_pairs(self, cursor=None):
        """Returns the number of potential homolog pairs in the database.
        Provide a specific cursor object to access the SQLite database if
        accessing this function through a transaction.
        """
        cnt_q = 'SELECT COUNT(*) FROM scores WHERE score > ?'
        if cursor is None:
            with sqlite3.connect(self.seqdb.db) as conn:
                cursor = conn.cursor()
                cursor.execute(cnt_q)
        else:
            cursor.execute(cnt_q)
        for row in cursor:
            return row[0]

    # Helper for index_potential_homolgs(): yields data values for the potential
    # homologs index.
    #
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

    def index_potential_homologs(self):
        with sqlite3.connect(self.seqdb.db) as conn:
            c = conn.cursor()
            # populate the seeds cache table:
            q = """
                INSERT OR REPLACE INTO %s (id, homologs)
                SELECT ?, IFNULL( (SELECT homologs FROM %s WHERE id = ?), "") || ?
            """ % (self.potential_homologs_table, self.potential_homologs_table)
            conn.cursor().executemany(q, self._give_potential_homologs(c))

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
        with sqlite3.connect(self.seqdb.db) as conn:
            c = conn.cursor()
            sys.stderr.write('Verifying database consistency at %s:' % self.seqdb.db)
            # potential homologs should have at most a single row per sequence.
            c.execute('SELECT COUNT(*) FROM %s' % self.potential_homologs_table)
            for row in c:
                num_records = row[0]
            c.execute("""
                SELECT COUNT(*) FROM (SELECT DISTINCT id FROM %s)
            """ % self.potential_homologs_table)
            for row in c:
                assert(num_records == row[0])

            c.execute('SELECT COUNT(*) FROM %s' % self.seeds_table)
            for row in c:
                num_records = row[0]
            c.execute("""
                SELECT COUNT(*) FROM (SELECT DISTINCT S_id,T_id,S_idx,T_idx FROM %s)
            """ % self.seeds_table)
            for row in c:
                assert(num_records == row[0])

            sys.stderr.write(' looks good!\n')

    def potential_homologs(self, seqid):
        """Given a sequence ID, returns a list of sequence IDs, all integers
        greater than ``seqid`` (to avoid duplicates in assembly), that are
        potentially homologous to the given sequence. The result is read
        directly from the ``potential_homologs_N`` table.
        """
        q = """
            SELECT homologs FROM %s WHERE id = ?
        """ % self.potential_homologs_table
        with sqlite3.connect(self.seqdb.db) as conn:
            c = conn.cursor()
            c.execute(q, (seqid,))
            for row in c:
                return [int(i) for i in row[0].split(',') if i]

        return []

    def seeds(self, S_id, T_id):
        """Given two sequence ids, finds all exactly matching segments
        (see :class:`oval.pw.Segment`) of length :attr:`wordlen` between the
        two. Segments are not necessarily in maximal form. For purposes of seed
        extension, however, we prefer to not have too many segments that are
        part of a one bigger segments (especially if they do not belong to an
        actual overlap alignment). This can be worked out by using
        :func:`maximal_seeds` which reduces any set of seeds into maximal,
        necessarily non-overlapping segments.

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
        with sqlite3.connect(self.seqdb.db) as conn:
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
                seeds += [pw.Segment(S_id=S_id, T_id=T_id, tx=tx)]

        return seeds


def plot_word_pvalues(index, path=None, num_bins=500):
    if path is None:
        path = 'word_pvalues.%d.png' % index.wordlen
    with sqlite3.connect(index.seqdb.db) as conn:
        c = conn.cursor()
        calc = index.word_pvalue_calculator(c)
        c.execute('select hits from %s' % index.words_table)
        log_pvalues = []
        for row in c:
            pvalue = calc(row[0].count('@'))
            if pvalue:
                log_pvalues += [log(pvalue)]

        plt.grid(True)
        plt.hist(log_pvalues, num_bins, normed=True, histtype='step', cumulative=True, color='k')
        plt.xlabel('Corrected log(p-value) of %d-words' % index.wordlen)
        plt.ylabel('Cumulative distribution')
        plt.savefig(path, dpi=300)


def maximal_seeds(seeds, S_id, T_id):
    """Given a list of exactly matching segments, reduces them into a list
    of *maximal* exactly matching segments in increasing order of ``S_idx``.
    A segment is in its maximal form if it cannot be extended in either
    direction by an exact match.

    Args:
        list[pw.Segment]: Exactly matching segments, potentially overlapping.
        S_id (int): The database ID of the "from" sequence.
        T_id (int): The database ID of the "to" sequence.

    Returns
        list[pw.Segment]: Maximal segments, guaranteed to not overlap.
    """
    seeds.sort(key=lambda s: s.tx.S_idx)
    # merge overlapping seeds:
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
                seeds[idx] = pw.Segment(S_id=S_id, T_id=T_id, tx=tx)
                seeds.pop(cand)
            else:
                cand += 1
        idx += 1
    return seeds
