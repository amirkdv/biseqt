import random
import sys
from hashlib import sha1
from bisect import bisect_left
from Bio import SeqIO, Seq, SeqRecord
import os.path
import sqlite3
from . import ffi, lib, CffiObject, ProgressIndicator

from itertools import chain


class Alphabet(object):
    """A sequence alphabet.

    Args:
        letters (str|List[str]): The letters of the alphabet.
    """
    def __init__(self, letters):
        if isinstance(letters, str):
            letters = [c for c in letters]
        self.letters = letters
        assert(len(set(len(l) for l in self.letters)) == 1)
        self.letlen = len(letters[0])

    def __len__(self):
        return len(self.letters)

    def randstr(self, length, **kw):
        """Generates a random string of the specified length from this
        alphabet. Optionally a probability distribution of letters may be
        specified.

        Note:
            This function and all its dependents (i.e :func:`randseq` and
            :func:`Sequence.mutate`) should not be used in condensed alphabets
            since they may generate invalid sequences like ``A1T3T5``.

        Args:
            length(int): Length of generated string in number of letters (and
                not necessarily number of characters)

        Keyword Args:
            letters_dist(Optional[dict]): The probability distribution of
                letters as a list of probabilities in order of letters in
                :attr:`letters`. Default is uniform.
        Returns:
            str: A random string.
        """
        cummulative_dist = kw.get('letters_dist', [
            1.0/len(self.letters) for _ in range(len(self.letters))
        ])
        for idx, prob in enumerate(cummulative_dist):
            cummulative_dist[idx] = cummulative_dist[idx-1] + prob if idx else prob
        return ''.join(self.letters[bisect_left(cummulative_dist, random.randint(0, 1000)/1000.0)] for _ in range(length))

    def randseq(self, length, **kw):
        """Generates a random :class:`Sequence` of the specified length from
        this alphabet. All arguments are identical to :func:`randstr`.

        Returns:
            seq.Sequence: A random sequence.
        """
        return Sequence(self.randstr(length, **kw), self)


class Sequence(object):
    """Wraps a C ``char[]`` and its corresponding ``int*`` of letter indices.
    Note that this is *not* a :class:`CffiObject` subclass since there's no
    underlying C struct. String indexing and slicing are supported::

        x = seq.Sequence('A2C3A2', seq.Alphabet(['A2', 'C3']))
        x[1]   # => 'C3'
        x[:2]  # => 'C3A2'
        x[:-2] # => 'A2'

    Attributes:
        alphabet (seq.Alphabet): The alphabet this sequence belongs to.
        c_idxseq  (cffi.cdata): The C ``int*`` which is actually used for all
            operations involving alignments; this is an array containing the
            indices of each letter in the sequence in the alphabet.
    """
    def __init__(self, letlist, alphabet):
        """Initializes the sequence object: translates all letters to integers
        corresponding to the position of each letter in the alphabet.

        Args:
            letlist(str|list): The contents of the sequence as a sequence of
                letters. In the special case where the alphabet letter length is
                1 a string will also be accepted.
            alphabet(seq.Alphabet): The alphabet to which this sequence belongs.
        """
        if not isinstance(letlist, list) and not isinstance(letlist, str):
            raise ValueError('`letlist` must be a list of letters.')
        letlen = len(alphabet.letters[0])
        if isinstance(letlist, str):
            if letlen != 1:
                msg = '`letlist` can only be a string if the alphabet' + \
                    'letter length is 1 (it is %d).' % alphabet.letlen
                raise ValueError(msg)
        self.letlist = list(letlist)
        self.length = len(self.letlist)
        self.alphabet = alphabet
        # build the sequence as an int array of positions in alphabet
        let_pos = {let: pos for pos, let in enumerate(alphabet.letters)}
        self.c_idxseq = ffi.new('int[]', [let_pos[let] for let in letlist])

    def __repr__(self):
        N, L = self.length, self.alphabet.letlen
        return ''.join([self.__getitem__(i) for i in range(self.length)])

    def __len__(self):
        return self.length

    def __getitem__(self, key):
        if isinstance(key, int):
            if key < self.length:
                start, finish = key, key + 1
            else:
                raise IndexError('Sequence index out of range')
        elif isinstance(key, slice):
            start, finish, _ = key.indices(self.length)
        else:
            raise TypeError(
                'Sequence indices must be integers not {}' % type(key).__name__
            )

        return ''.join(self.letlist[start:finish])

    def mutate(self, **kw):
        """Returns a mutant of this sequence with specified probabilities. The
        sequence is scanned and copied to the mutated sequence where at each
        position:

        - the current letter will be replaced by an arbitrary letter with a
          distribution that depends on the original letter.
        - with a certain fixed probability a gap may be openned with random
          length with geometric distribution.

        Accordingly, a transcript is generated which corresponds to the
        performed edit sequence.

        Note that there is a bound on the precision of the effective
        proabilities, see ``precision``.

        Keyword Args:
            subst_probs(List[List[float]]): the probability distribution for
                each pair of possible substitutions such that
                :attr:`c_obj.subst_probs[i][j]` is the probability of lette
                *i* (integer index) of the alphabet being substituted by letter
                *j* (integer index).
            go_prob(Optional[float]): probability of a gap starting at any
                position, default 0 (use 1 for linear gap penalty).
            ge_prob(Optional[float]): Bernoulli success probability of the gap
                extension distribution, default 0.
            insert_dist(Optional[List[float]]): the distribution passed to
                :func:`Alphabet.randstr` when inserting arbitrary strings;
                default is uniform.
            precision(Optional[float]): As in :func:`Alphabet.randstr` and the
                same applies to gap probabilities. Default is 0.001.

        Returns:
            tuple: The mutant (a :class:`Sequence`)and corresponding transcript
                (an :class:`biseqt.pw.Transcript`).
        """
        subst_probs = kw['subst_probs']
        go_prob, ge_prob = kw.get('go_prob', 0), kw.get('ge_prob', 0)
        assert(go_prob <= ge_prob)
        insert_dist = kw.get('insert_dist', None)
        precision = kw.get('precision', 0.001)
        N = 1/precision
        assert(precision > 0)
        assert(go_prob < 1)
        assert(ge_prob < 1)
        T, opseq, op, k = '', '', None, 0
        while k < self.length:
            if op:
                opseq += op
            if op == 'D' and random.randint(0, N) < ge_prob * N:
                # with probability ge_prob extend the deletion stretch:
                op, k = 'D', k + 1
                continue
            elif op == 'I' and random.randint(0, N) < ge_prob * N:
                # with probability ge_prob extend the insertion stretch:
                op, k = 'I', k
                T += self.alphabet.randstr(
                    1, letter_dist=insert_dist, precision=precision
                )
                continue
            else:
                go_cutoff = go_prob * N / 2.0
                if str(op) not in 'ID' and random.randint(0, N) < go_cutoff:
                    op, k = 'D', k + 1
                elif str(op) not in 'ID' and random.randint(0, N) < go_cutoff:
                    op, k = 'I', k
                    T += self.alphabet.randstr(
                        1, letter_dist=insert_dist, precision=precision
                    )
                else:
                    # math/sub with probability 1 - g_o
                    T += self.alphabet.randstr(
                        1,
                        letters_dist=subst_probs[self.c_idxseq[k]],
                        precision=precision
                    )
                    op, k = ('M' if T[-1] == self[k] else 'S'), k + 1

        opseq += op
        return (Sequence(T, self.alphabet), opseq)

    def randread(self, **kw):
        """Generates a random collection of lossy reads from the current
        sequence. Each read is generated by reading a substring with Gaussian
        length and mutating it with given probabilites. Additionally, two
        substrings, both of lenght ``len_mean``, are included one covering the
        very beginning and the other the very end of the sequence. For
        example::

            subst_probs = { ... }
            N = 10000
            genome = seq.Alphabet('ACGT').randseq(N)
            reads = genome.randread(subst_probs=subst_probs, coverage=40)

        Keyword Args:
            subst_probs(List[List[float]]): As in :func:`~mutate`.
            go_prob(Optional[float]): As in :func:`~mutate`.
            ge_prob(Optional[float]): As in :func:`~mutate`.
            coverage (Optional[float]): The expected number of times each
                letter in the sequence appears in the entire read collection.
                Default is 5.
            len_mean (Optional[float]): The mean of the normal distribution of
                read lengths. Default is 100.
            len_sd (Optional[float]): The standard deviation of the normal
                distribution or read lengths. Default is 1.

        Yields:
            tuple: of :class:`Sequence` and ``int`` (starting position).
        """
        subst_probs = kw['subst_probs']
        go_prob, ge_prob = kw.get('go_prob', 0), kw.get('ge_prob', 0)
        coverage = kw.get('coverage', 5)
        len_mean, len_sd = kw.get('len_mean', 100), kw.get('len_sd', 1)
        N = self.length
        num = int(1.0 * N * coverage/len_mean)
        msg = 'generating %d random reads with length ~ ' +  \
            '(mean=%.2f, s.d.=%.2f) for a %d nucl. long genome'
        indicator = ProgressIndicator(msg % (num + 2, float(len_mean), float(len_sd), N), num + 2)
        indicator.start()

        # include a read that reaches the begenning:
        read, _ = Sequence(self[:len_mean], self.alphabet).mutate(**kw)
        yield read, 0
        indicator.progress()

        for i in range(num):
            # minimum read leangth is 10, and max is N-1
            length = max(10, min(N-1, int(random.gauss(len_mean, len_sd))))
            start = random.randint(0, N-length)
            read = ''.join([self[k] for k in range(start, start + length)])
            read = Sequence(read, self.alphabet)
            read, _ = read.mutate(**kw)
            yield read, start
            indicator.progress()

        # include a read that reaches the end:
        read, _ = Sequence(self[-len_mean:], self.alphabet).mutate(**kw)
        yield read, N - len_mean
        indicator.progress()
        indicator.finish()


READ = 0
REFERENCE = 1

class SeqDB(object):
    """Wraps an SQLite database containing sequences and potentiall word indices
    (see :class:`biseqt.words.Index`) or assembly data structures. The sequences
    are stored in a ``seq`` table.

    Attributes:
        db (string): Path to the SQLite datbase.
        alphabet (Alphabet): The alphabet for sequences in the database.
    """
    def __init__(self, db, alphabet=None):
        assert isinstance(alphabet, Alphabet)
        self.alphabet, self.db = alphabet, db
        # populate the seqinfo cache
        self.seqinfo(use_cache=False)

    def initdb(self):
        """Initializes the database: creates the required tables and
        fails if any of them already exists."""
        sys.stderr.write('Initializing sequence DB at: %s\n' % self.db)
        with sqlite3.connect(self.db) as conn:
            c = conn.cursor()
            sys.stderr.write("Turning off SQLite's journaling for performance: do not trust the contents of the database after a crash.\n")
            c.execute('PRAGMA journal_mode = OFF;')
            q = """
                CREATE TABLE seq (
                  id   INTEGER PRIMARY KEY ASC,
                  name VARCHAR(40), -- content address (SHA1)
                  type INTEGER,     -- READ or REFERENCE (module constants)
                  orig TEXT,        -- original sequence identifier
                  rc   INTEGER DEFAULT 0, -- whether sha(seq) == name (0) or sha(seq.rc) == name (1)
                  seq  TEXT         -- the contents of the sequence
                );
            """
            c.execute(q)

    def populate(self, fasta_src, seq_type):
        """Given a FASTA source file, loads all the sequences into the database.

        Args:
            fasta_src(str): Path to FASTA source.
            seq_type(int): One of the READ or REFERENCE module constants
                indicating the type of the provided sequences.
        """
        sys.stderr.write('Loading sequences from: %s\n' % fasta_src)

        # leave the rc reference empty for now
        q = """
            INSERT INTO seq (name, type, orig, rc, seq) VALUES (?, ?, ?, ?, ?)
        """
        name = lambda x: sha1(str(x)).hexdigest()
        seqgen = lambda: enumerate(SeqIO.parse(fasta_src, 'fasta'))
        # FIXME debug
        recs = chain(
            ((name(r.seq), seq_type, r.id, 0, str(r.seq)) for idx,r in seqgen()), # if idx < 100),
            ((name(r.seq), seq_type, r.id, 1, str(r.reverse_complement().seq)) for idx,r in seqgen()) # if idx < 100)
        )
        with sqlite3.connect(self.db) as conn:
            c = conn.cursor()
            c.executemany(q, recs)

        # update the seqinfo cache
        self.seqinfo(use_cache=False)

    def loadseq(self, seqid):
        """Loads a sequence given its internal numeric seqid.

        Args:
            seqid (int): Sequence ID as found in the ``seq`` table.

        Returns:
            Sequence
        """
        with sqlite3.connect(self.db) as conn:
            c = conn.cursor()
            c.execute('SELECT seq FROM seq WHERE id = ?', (seqid,))
            for row in c:
                return Sequence(str(row[0]), self.alphabet)

    # FIXME docs
    def seqinfo(self, use_cache=True):
        """Return a dict of metadata about all sequences keyed by sequence ids
        as found in the ``seq`` table. The output looks like this::

            A = Alphabet('ACGT')
            T = SeqDB('genome.db', alphabet=A)
            info = T.seqinfo()
            info[12] #=> {'start': 17, 'length': 479, 'name': u'R12'}

        Returns:
            dict: Metadata dicts in a dict keyed by sequence ID.
        """
        if not os.path.exists(self.db):
            # initdb() has not been called for this path.
            return None

        if use_cache and self._seqinfo is not None:
            return self._seqinfo

        self._seqinfo = {}
        with sqlite3.connect(self.db) as conn:
            c = conn.cursor()
            c.execute("SELECT id, name, type, orig, rc, LENGTH(seq) FROM seq")
            for row in c:
                self._seqinfo[row[0]] = {
                    'name':   row[1],
                    'type':   row[2],
                    'orig':   row[3],
                    'rc':     row[4],
                    'length': row[5],
                }

        return self._seqinfo

    # FIXME docs
    def seqids(self, seq_type='both'):
        """Returns a list of all seqids found in the ``seq`` table.

        Returns:
            List[int]
        """
        condition = 'WHERE type = %d' % int(seq_type) if seq_type != 'both' else ''
        with sqlite3.connect(self.db) as conn:
            c = conn.cursor()
            c.execute('SELECT id FROM seq ' + condition)
            return [row[0] for row in c]

# FIXME this should become a method in SeqDB
def make_sequencing_fixture(genome_file, reads_file, genome_length=1000, **kw):
    """Helper method for tests. Generates a random genome and a random sequence
    of reads from it.

    Args:
        genome_file(str): Path to FASTA file to write genome to.
        reads_file(str): Path to FASTA file to write sequencing reads to.

    Keyword Args:
        genome_length(int): Length of random genome.

    Additionally, all keyword arguments of :func:`Sequence.randread` are
    allowed and passed as is.
    """
    A = Alphabet('ACGT')
    if 'go_prob' not in kw:
        kw['go_prob'] = 0.1
    if 'go_prob' not in kw:
        kw['ge_prob'] = 0.5
    G = A.randseq(genome_length)
    seqrec = SeqRecord.SeqRecord(
        Seq.Seq(str(G)), id='genome', description="(full correct genome)"
    )
    SeqIO.write([seqrec], genome_file, 'fasta')
    readrecs = []
    for idx, (read, start) in enumerate(G.randread(**kw)):
        seqid = 'R%s_P%d' % (str(uuid4())[:8], start)
        readrecs += [
            SeqRecord.SeqRecord(Seq.Seq(str(read)), id=seqid)
        ]
    sys.stderr.write('saving reads to %s.\n' % reads_file)
    SeqIO.write(readrecs, reads_file, 'fasta')
