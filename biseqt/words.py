"""Provides tools for *k*-mer analysis."""

import sys
import os.path
import re
import sqlite3
from math import sqrt, erf, log
from scipy.special import erfinv
from Bio import SeqIO
from collections import namedtuple
from matplotlib import pyplot as plt
from . import pw, seq, ProgressIndicator, CffiObject, ffi, lib


class Index(object):
    """The main responsibility of an Index is to respond to :func:`seeds`
    with a list of :class:`biseqt.pw.Segment` s given two sequence IDs.
    This is done by recording all "words" observed in all sequences
    in the database during :func:`index`. Each index is uniquely defined by
    its :attr:`seqdb` and its :attr:`wordlen`. Each such index operates
    using 3 database tables with the following schema (replace `N`
    with :attr:`wordlen`). All tables are populated, in the order presented
    below, by :func:`index`.

    .. code-block:: sql

        CREATE TABLE words_N (
          'tuple' integer primary key, -- decimal represetnation of N-tuple taken as a number
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
        self.sensitivity_erfinv = erfinv(kwargs.get('band_sensitivity', 0.9))
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
            c.execute('CREATE index tuple ON %s(tuple)' % self.t_words)
            q = """
                CREATE TABLE %s (
                  'S_id' integer REFERENCES seq(id),
                  'T_id' integer REFERENCES seq(id),
                  'S_idx' integer,
                  'T_idx' integer
                );
            """ % (self.t_seeds)
            c.execute(q)
            # NOTE cf. most_significant_shift
            #q = """
                #CREATE TABLE %s (
                  #'S_id' integer REFERENCES seq(id),
                  #'T_id' integer REFERENCES seq(id),
                  #'shift' integer, -- shift of seed (S_idx - T_idx)
                  #'score' real, --  = -log(p-value)
                  #PRIMARY KEY (S_id, T_id, shift)
                #);
            #""" % (self.t_scores)
            #c.execute(q)
        sys.stderr.write('Initialized index tables %s, %s, %s.\n'
            % (self.t_words, self.t_seeds, self.t_scores)
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
            tup = sum(digits[x]*(len(digits)**i) for x,i in zip(tup,reversed(range(len(tup)))))
            yield (tup, idx)

    def index(self):
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
        # NOTE "OR REPLACE" means "ON CONFLICT REPLACE". We need a conflict to
        # occur for this one-shot query to work so tuple must have a unique
        # constraint on it.
        hit_ins_q = """
            INSERT OR REPLACE INTO %s (tuple, hits)
            SELECT ?, IFNULL( (SELECT hits FROM %s WHERE tuple = ?), "") || ?
        """ % (self.t_words, self.t_words)
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
            c.execute('SELECT tuple, hits from %s' % self.t_words)
            q = """
                INSERT INTO %s (S_id, T_id, S_idx, T_idx)
                VALUES (?, ?, ?, ?)
            """ % self.t_seeds
            conn.cursor().executemany(q, self.scan_seeds(c, conn.cursor()))

            sys.stderr.write('Creating SQL index on %s'  % self.t_seeds)
            c.execute("""
                CREATE INDEX seeds_ids ON %s (S_id, T_id)
            """ % self.t_seeds)
            sys.stderr.write('.\n')
            # FIXME doesn't it autocommit?
            conn.commit()

    def word_pvalue_calculator(self):
        N = self.num_words

        with sqlite3.connect(self.seqdb.db) as conn:
            c = conn.cursor()
            c.execute('select sum(length(seq)) from seq')
            L = int(c.next()[0]) # total sequence length of all reads

        # Normal approximation (mean and standard deviation) of a binomial distribution
        B2N = lambda n, p: (n*p, sqrt(n * p * (1-p)))
        prob_word = lambda length: (1.0/len(self.seqdb.alphabet)) ** length
        mu, sd = B2N(L, prob_word(self.wordlen)) # approximating normal distribution of word counts

        pvalue_calculator = lambda x: N * 0.5 * (1 - erf((x - mu) / (sd * sqrt(2))))
        return pvalue_calculator

    @property
    def num_words(self):
        with sqlite3.connect(self.seqdb.db) as conn:
            c = conn.cursor()
            c.execute('select count(*) from %s' % self.t_words)
            return int(c.next()[0])

    def band_radius(self, readlens, d, g):
        # max alignment length
        L = min(readlens[0] - d, readlens[1]) + min(d, 0)
        # expected alignment length
        K = (2.0/(2-g)) * L
        r = 2 * self.sensitivity_erfinv * sqrt(g * (1-g) * K)
        return int(r)

    def seed_pvalue_contribution(self, readlens, shift, gap_prob=None):
        assert(gap_prob > 0 and gap_prob < 1)
        r = self.band_radius(readlens, shift, gap_prob)
        A = 2 * r * sqrt(sum((l - abs(shift))**2 for l in readlens))
        # (radius, contribution of one seed, starting value)
        return r, sum(log(l) for l in readlens) - log(A), log(sum(readlens))

    # Helper for index(): yields data values to be inserted in the seeds index.
    def scan_seeds(self, hit_recs, cursor):
        word_pvalue_calc = self.word_pvalue_calculator()
        seqinfo = self.seqdb.seqinfo()
        N = self.num_words
        indicator = ProgressIndicator(
            'Indexing seeds from %d observed %d-mers' % (N, self.wordlen), N
        )
        indicator.start()
        for tup, hits in hit_recs:
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
                    S_id, T_id = int(S_hit[0]), int(T_hit[0])
                    S_idx, T_idx = int(S_hit[1]), int(T_hit[1])
                    if seqinfo[S_id]['name'] == seqinfo[T_id]['name']:
                        # don't care about seeds from S to S or S'
                        continue

                    # TODO see most_significant_shift below.
                    #shift = S_hit[1] - T_hit[1]
                    #lengths = [
                        #seqinfo[S_id]['length'],
                        #seqinfo[T_id]['length']
                    #]
                    # (radius, contribution of one seed, base)
                    #r, s, s0 = self.seed_pvalue_contribution(lengths, shift)
                    #q = """INSERT OR REPLACE INTO %s (S_id, T_id, shift, score)
                        #SELECT ?, ?, ?, IFNULL(
                            #(SELECT ? + score FROM %s
                             #WHERE S_id = ? AND T_id = ? AND shift = ?
                            #), ?
                        #)
                    #""" % (self.t_scores, self.t_scores)
                    #_give_rec = lambda d: (S_id, T_id, d, s, S_id, T_id, d, s0)
                    #cursor.executemany(q, (_give_rec(d) for d in range(shift-r, shift+r+1)))
                    yield (S_id, T_id, S_idx, T_idx)

        indicator.finish()

    # FIXME takes too long (also cf. scan_seeds) ; see discovery.most_significant_shift
    # Would HDF make this faster? http://docs.h5py.org
    #def most_significant_shift(self, S_id, T_id):
        #with sqlite3.connect(self.seqdb.db) as conn:
            #q = """
                #SELECT shift, score FROM %s
                #WHERE S_id = ? AND T_id = ? ORDER BY score DESC LIMIT 1
            #""" % self.t_scores
            #c = conn.cursor()
            #c.execute(q, (S_id, T_id))
            #row = c.next()
            #return (int(row[0]), float(row[1]))

    def seeds(self, S_id, T_id):
        """Given two sequence ids, finds all exactly matching segments
        (see :class:`biseqt.pw.Segment`) of length :attr:`wordlen` between the
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
        """ % (self.t_seeds)
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
        calc = index.word_pvalue_calculator()
        log_pvalues = []
        N = index.num_words
        indicator = ProgressIndicator('Calculating p-values for all %d %d-mers' % (N, index.wordlen), N)
        indicator.start()
        c.execute('select hits from %s' % index.t_words)
        for row in c:
            pvalue = calc(row[0].count('@'))
            indicator.progress()
            if pvalue:
                log_pvalues += [log(pvalue)]

        indicator.finish()

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
            # TODO rename this shift_X to something else, this is not "shift" as
            # everywhere else
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
