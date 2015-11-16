import random
from Bio import SeqIO, Seq, SeqRecord

from . import ffi, lib, CffiObject

class Alphabet(CffiObject):
    """Wraps a C ``sequence_alphabet*``.

    Attributes:

        c_obj (cffi.cdata): points to a ``sequence_alphabet`` struct.

    Args:
        letters (str|List[str]): The letters of the alphabet.

    Raises:
        AssertionError: If all the letters of the alphabet are not of the same
            length.
    """
    def __init__(self, letters):
        if isinstance(letters, str):
            letters = [c for c in letters]
        assert(len(set([len(s) for s in letters])) == 1)
        # each letter string in the alphabet must be "owned" by an object
        # that's kept alive.
        self._c_letters_ka = [ffi.new('char[]', letters[i]) for i in range(len(letters))]
        self._c_alph_ka = ffi.new('char *[]', self._c_letters_ka)
        self.c_obj = ffi.new('sequence_alphabet*', {
            'length': len(letters),
            'letter_length': len(letters[0]),
            'letters': self._c_alph_ka
        })

    def __len__(self):
        return len(self.letters)

    def __getattr__(self, name):
        if name == 'letters':
            N, L = self.length, self.letter_length
            return [''.join([self.c_obj.letters[i][j] for j in range(L)]) for i in range(N)]
        else:
            return super(Alphabet, self).__getattr__(name)

    def randstr(self, length, **kw):
        """Generates a random string of the specified length from this
        alphabet. Optionally a probability distribution of letters may be
        specified.

        Note:
            Arbitrary precision on probability distributions is not supported.
            The default precision is 0.001, modify this with :arg:`precision`.

        Args:
            length(int): Length of generated string in number of letters (and
                not necessarily number of characters)

        Keyword Args:
            letters_dist(Optional[dict]): The probability distribution of
                letters as a list of probabilities in order of letters in
                :attr:`c_obj.letters`. Default is uniform.

            precision(Optional[float]): The maximum precision of the probability
                distribution. Default is 0.001.
        Returns:
            str: A random string.
        """
        letters_dist = kw.get('letters_dist', [1.0/self.length for _ in range(self.length)])
        precision = kw.get('precision', 0.001)
        assert(precision > 0)
        assert(abs(1-sum(letters_dist)) < precision)
        space = []
        for i in range(self.length):
            N = 1.0/precision
            space += [ffi.string(self._c_letters_ka[i])] * int(N * length * letters_dist[i])
        return ''.join(random.sample(space, length))


    def randseq(self, length, **kw):
        """Generates a random :class:`Sequence` of the specified length from
        this alphabet. All arguments are identical to :func:`randstr`.

        Returns:
            seq.Sequence: A random sequence.
        """
        return Sequence(self.randstr(length, **kw), self)

class Sequence():
    """Wraps a C ``char[]`` and its corresponding ``int*`` of letter indices.
    Note that this is *not* a :class:`CffiObject` subclass since there's no
    underlying C struct. String indexing and slicing are supported::

        x = seq.Sequence('A2C3A2', seq.Alphabet(['A2', 'C3']))
        x[1]   # => 'C3'
        x[:2]  # => 'C3A2'
        x[:-2] # => 'A2'

    Attributes:
        alphabet (seq.Alphabet): The alphabet this sequence belongs to.
        c_charseq (cffi.cdata): Points to the underlying C ``char[]``.
        c_idxseq  (cffi.cdata): The C ``int*`` which is actually used for all
            operations involving alignments; this is an array containing the
            indices of each letter in the sequence in the alphabet.
    """
    def __init__(self, string, alphabet):
        assert(len(string) % alphabet.letter_length == 0)
        self.alphabet = alphabet
        self.length = len(string)/alphabet.letter_length
        self.c_charseq = ffi.new('char[]', string)
        self.c_idxseq = lib.idxseq_from_charseq(alphabet.c_obj, self.c_charseq, self.length)

    def __repr__(self):
        N, L = self.length, self.alphabet.letter_length
        return ''.join([self.__getitem__(i) for i in range(self.length)])

    def __len__(self):
        return self.length

    def __getitem__(self, key):
        if isinstance(key, int):
            if key < self.length:
                s, f = key, key + 1
            else:
                raise IndexError('Sequence index out of range')
        elif isinstance(key, slice):
            s, f, _ = key.indices(self.length)
        else:
            raise TypeError('Sequence indices must be integers not {}'.format(type(key).__name__))
        return ffi.string(self.c_charseq)[s*self.alphabet.letter_length: f*self.alphabet.letter_length]

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
        proabiliites, see :arg:`precision`.

        Keyword Args:
            subst_probs(List[List[float]]): the probability distribution for
                each pair of possible substitutions such that
                :attr:`c_obj.subst_probs[i][j]` is the probability of letter *i*
                (integer index) of the alphabet being substituted by letter *j*
                (integer index).
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
                (an :class:`align.pw.Transcript`).
        """
        subst_probs = kw['subst_probs']
        go_prob, ge_prob =  kw.get('go_prob', 0), kw.get('ge_prob', 0)
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
                T += self.alphabet.randstr(1,
                    letter_dist=insert_dist,
                    precision=precision
                )
                continue
            else:
                # with probability go_prob start a gap unless previous op was
                # a gap itself.
                if str(op) not in 'ID' and random.randint(0, N) < ge_prob * N / 2.0:
                    op, k = 'D', k + 1
                elif str(op) not in 'ID' and random.randint(0,N) < ge_prob * N / 2.0:
                    op, k = 'I', k
                    T += self.alphabet.randstr(1,
                        letter_dist=insert_dist,
                        precision=precision
                    )
                else:
                    # math/sub with probability 1 - g_o
                    T += self.alphabet.randstr(1,
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
        very beginning and the other the very end of the sequence. For example::

            subst_probs = {
                ... # snip
            }
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
            len_var (Optional[float]): The variance of the normal
                distribution or read lengths. Default is 1.

        Yields:
            tuple: the read (a :class:`Sequence`) and its starting position (``int``).
        """
        subst_probs = kw['subst_probs']
        go_prob, ge_prob =  kw.get('go_prob', 0), kw.get('ge_prob', 0)
        coverage = kw.get('coverage', 5)
        len_mean, len_var = kw.get('len_mean', 100), kw.get('len_var', 1)
        N = self.length

        # include a read that reaches the begenning:
        read, _ = Sequence(self[:len_mean], self.alphabet).mutate(**kw)
        yield read, 0

        num = int(1.0 * N * coverage/len_mean)
        for i in range(num):
            # minimum read leangth is 10, and max is N-1
            length = max(10, min(N-1, int(random.gauss(len_mean, len_var))))
            start = random.randint(0, N-length)
            x = Sequence(''.join([self[k] for k in range(start, start + length)]), self.alphabet)
            read, _ = x.mutate(**kw)
            yield read, start

        # include a read that reaches the end:
        read, _ = Sequence(self[-len_mean:], self.alphabet).mutate(**kw)
        yield read, N - len_mean

def make_sequencing_fixture(genome_file, reads_file, genome_length=1000, **kw):
    """Helper method for tests. Generates a random genome and a random sequence
    of reads from it.

    Args:
        genome_file(str): Path to FASTA file to write genome to.
        reads_file(str): Path to FASTA file to write sequencing reads to.

    Keyword Args:
        genome_length(int): Length of random genome.

    Additionally, all keyword arguments of :func:`Sequence.randread` are allowed
    and passed as is.
    """
    A = Alphabet('ACGT')
    if 'go_prob' not in kw:
        kw['go_prob'] = 0.1
    if 'go_prob' not in kw:
        kw['ge_prob'] = 0.5
    if 'subst_probs' not in kw:
        kw['subst_probs'] = [[0.94 if k==i else 0.02 for k in range(4)] for i in range(4)]
    G = A.randseq(genome_length)
    seqrec = SeqRecord.SeqRecord(Seq.Seq(str(G)), id='genome', description="(full correct genome)")
    SeqIO.write([seqrec], genome_file, 'fasta')
    readrecs = []
    for idx, (read, start) in enumerate(G.randread(**kw)):
        readrecs += [
            SeqRecord.SeqRecord(Seq.Seq(str(read)), id="R%d_P%d" % (idx+1, start))
        ]
    SeqIO.write(readrecs, reads_file, 'fasta')
