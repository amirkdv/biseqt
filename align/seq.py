import numpy as np
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

    def randstr(self, length, letters_dist=None):
        """Generates a random string of the specified length from this
        alphabet. Optionally a probability distribution of letters may be
        specified.

        Args:
            length(int): Length of generated string in number of letters (and
                not necessarily number of characters)

            letters_dist(Optional[dict]): The probability distribution of
                letters as a list of probabilities in order of letters in
                :attr:`c_obj.letters`. Default is uniform.

        Returns:
            str: A random string.
        """
        if letters_dist is None:
            letters_dist = [1.0/self.length for _ in range(self.length)]

        assert(abs(1-sum(letters_dist)) < 0.001)
        space = []
        for i in range(self.length):
            # NOTE this effectively sets the max precision of subst_probs to .01
            space += [ffi.string(self._c_letters_ka[i])] * int(100 * length * letters_dist[i])
        return ''.join(random.sample(space, length))


    def randseq(self, length, letters_dist=None):
        """Generates a random :class:`Sequence` of the specified length from
        this alphabet. Optionally a probability distribution of letters may be
        specified.

        Args:
            length(int): Length of generated string in number of letters (and
                not necessarily number of characters)

            letters_dist(Optional[dict]): The probability distribution of
                letters as a list of probabilities in order of letters in
                :attr:`c_obj.letters`. Default is uniform.

        Returns:
            seq.Sequence: A random sequence.
        """
        return Sequence(self.randstr(length, letters_dist), self)

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
        c_idxseq  (cffi.cdata): points to the actually used ``int*``.
    """
    def __init__(self, string, alphabet):
        assert(len(string) % alphabet.letter_length == 0)
        self.alphabet = alphabet
        self.length = len(string)/alphabet.letter_length
        global lib
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

    def mutate(self, go_prob=0, ge_prob=0, subst_probs=None, insert_dist=None):
        """Returns a mutant of this sequence with specified probabilities. The
        sequence is scanned and copied to the mutated sequence where at each
        position:

        - the current letter will be replaced by an arbitrary letter with a
          distribution that depends on the original letter.
        - with a certain fixed probability a gap may be openned with random
          length with geometric distribution.

        Accordingly, a transcript is generated which corresponds to the
        performed edit sequence.

        Args:
             S(seq.Sequence): original sequence.
             go_prob(Optional[float]): probability of a gap starting at any
                 position, default 0 (use 1 for linear gap penalty).
             ge_prob(Optional[float]): Bernoulli success probability of the gap
                extension distribution, default 0.
             subst_probs(List[List[float]]): the probability distribution for
                each pair of possible substitutions such that
                :attr:`c_obj.subst_probs[i][j]` is the probability of letter *i*
                (integer index) of the alphabet being substituted by letter *j*
                (integer index).
             insert_dist(List[float]): the distribution passed to
                 :func:`Alphabet.randstr` when inserting arbitrary strings;
                 default is uniform.

        Returns:
            tuple: The mutant (a :class:`Sequence`)and corresponding transcript
                (an :class:`align.pw.Transcript`).
        """
        assert(go_prob < 1)
        assert(ge_prob < 1)
        T = ''
        k = 0
        opseq = ''
        while k < self.length:
            if go_prob:
                # NOTE max precision for gap_open is .01
                if random.randint(0, 100) < go_prob * 100:
                    # deletion of length with geometric distribution, but not
                    # more than we can actually delete:
                    length = min(np.random.geometric(1 - ge_prob), self.length - k)
                    opseq += 'D' * length
                    k += length
                    continue
                if random.randint(0, 100) < go_prob * 100:
                    length = np.random.geometric(1 - ge_prob)
                    # insert of length with geometric distribution:
                    opseq += 'I' * length
                    T += self.alphabet.randstr(length, insert_dist)
                    continue
            # no gap, substitute:
            T += self.alphabet.randstr(1, subst_probs[self.c_idxseq[k]])[0]
            opseq += 'M' if T[-1] == self[k] else 'S'
            k += 1
        return (Sequence(T, self.alphabet), opseq)

    def randread(self, subst_probs=None, go_prob=0, ge_prob=0, coverage=5, len_mean=100, len_var=25):
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

        Args:
            subst_probs(List[List[float]]): As in :func:`~mutate`.
            go_prob(Optional[float]): As in :func:`~mutate`.
            ge_prob(Optional[float]): As in :func:`~mutate`.
            coverage (Optional[float]): The expected number of times each
                letter in the sequence appears in the entire read collection.
                Default is 5.
            len_mean (Optional[float]): The mean of the normal distribution of
                read lengths.
            len_var (Optional[float]): The variance of the normal
                distribution or read lengths.

        Yields:
            tuple: the read (a :class:`Sequence`) and its starting position (``int``).
        """
        N = self.length

        # include a read that reaches the begenning:
        yield Sequence(self[:len_mean], self.alphabet), 0

        num = int(1.0 * N * coverage/len_mean)
        for i in range(num):
            # minimum read leangth is 10, and max is N-1
            length = max(10, min(N-1, int(np.random.normal(len_mean, len_var))))
            start = np.random.randint(0, N-length)
            x = Sequence(''.join([self[k] for k in range(start, start + length)]), self.alphabet)
            read, _ = x.mutate(go_prob=go_prob, ge_prob=ge_prob, subst_probs=subst_probs)
            yield read, start

        # include a read that reaches the end:
        yield Sequence(self[-len_mean:], self.alphabet), N - len_mean

def make_sequencing_fixture(genome_file, reads_file, genome_length=1000, **kw):
    """Helper method for tests. Generates a random genome and a random sequence
    of reads from it. All keyword arguments are passed as is to
    :func:`Sequence.randread`.

    Args:
         genome_file(str): Path to FASTA file to write genome to.
         reads_file(str): Path to FASTA file to write sequencing reads to.
         genome_length(int): Length of random genome.
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
