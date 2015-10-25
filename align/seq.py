import numpy as np
import random

from . import ffi, lib, utils, CffiObject

class Alphabet(CffiObject):
    """Wraps a C `sequence_alphabet*`.

    Attributes:
        c_obj (cffi.cdata): points to a sequence_alphabet struct.
        _c_letters_ka (list[cffi.cdata]): has ownership of (keeps alive) C
            pointers to each letter (which is a `char[]`) of the alphabet.
        _c_alph_ka (cffi.cdata): has ownership of (keeps alive) the C
            pointer to the full substitution matrix.
    """
    def __init__(self, alphabet):
        if isinstance(alphabet, str):
            alphabet = [c for c in alphabet]
        assert(len(set([len(s) for s in alphabet])) == 1)
        # each letter string in the alphabet must be "owned" by an object
        # that's kept alive.
        self._c_letters_ka = [ffi.new('char[]', alphabet[i]) for i in range(len(alphabet))]
        self._c_alph_ka = ffi.new('char *[]', self._c_letters_ka)
        self.c_obj = ffi.new('sequence_alphabet*', {
            'length': len(alphabet),
            'letter_length': len(alphabet[0]),
            'letters': self._c_alph_ka
        })

    def __getattr__(self, name):
        if name == 'letters':
            N, L = self.length, self.letter_length
            return [''.join([self.c_obj.letters[i][j] for j in range(L)]) for i in range(N)]
        else:
            return super(Alphabet, self).__getattr__(name)


class Sequence(CffiObject):
    """Wraps a C `char[]` and keeps its length. Placeholder for potential
    additions.

    Attributes:
        alphabet (Alphabet)
        c_charseq (cffi.cdata): points to the underlying C char[].
        c_idxseq  (cffi.cdata): points to the actually used int*.
    """
    def __init__(self, string, alphabet):
        assert(len(string) % alphabet.letter_length == 0)
        global lib
        self.c_charseq = ffi.new('char[]', string)
        self.length = len(string)/alphabet.letter_length
        self.c_idxseq = lib.idxseq_from_charseq(alphabet.c_obj, self.c_charseq, self.length)
        self.alphabet = alphabet

    def __repr__(self):
        N, L = self.length, self.alphabet.letter_length
        return ''.join([self.__getitem__(i) for i in range(self.length)])

    def __getitem__(self, key):
        if isinstance(key, int) and key < self.length:
            return ''.join([self.alphabet._c_letters_ka[self.c_idxseq[key]][i] for i in range(self.alphabet.letter_length)])
        else:
            raise KeyError

def rawrand(length, dist):
    space = []
    for k in dist.keys():
        # NOTE this effectively sets the maximum precision of error rates to .01
        space += [k] * int(100 * length * dist[k])
    return ''.join(random.sample(space, length))

def randseq(length, alphabet, dist=None):
    """Generates a random sequence of the specified length within the alphabet
    specified by the keys in the distribution matrix. The distribution matrix
    should look like this for nucleotide sequences:
        {'A': 0.25,
         'C': 0.25,
         'G': 0.25,
         'T': 0.25}

    :param length(int): length of generated sequence.
    :param dist(dict): Optional; keys are the letters in the alphabet and values
        are the probability of an arbitrary nucleatoride being each letter.
        Default is uniform.
    """
    if dist is None:
        letters, L = alphabet.letters, alphabet.length
        dist = {k:1.0/L for k in letters}

    assert abs(1-sum(dist.values())) < 0.001
    return Sequence(rawrand(length, dist), alphabet)

# TODO support hompolymeric-specific gap parameters
def mutate(S, gap_open=0.1, gap_continue=0.5, error_rates=None, insert_dist=None):
    """Mutates a given sequence with specified probabilities. The sequence is
    scanned and copied to the mutated sequence where at each position:
    * the current letter will be replaced by an arbitrary letter with a
        distribution that depends on the original letter.
    * with a certain fixed probability a gap may be openned with random length
        with geometric distribution.
    Accordingly, an opseq (see `align.solve()`) is generated which corresponds
    to the performed edit sequence.

    :param S(seq.Sequence): original sequence.
    :param gap_open(float): probability of a gap starting at any position.
    :param gap_continue(float): Bernoulli success probability of the gap
        extension distribution
    :param error_rates(dict): a letter-by-letter error probability matrix. For
        nucleotides, for example, it should have the following structure:

            {'A':{'A': 0.7,
                  'C': 0.1,
                  'G': 0.1,
                  'T': 0.1},
             'C':{...
             ...
            }
    :param insert_dist(dict): the distribution passed to randseq() when
        inserting arbitrary strings; default is uniform.
    """
    assert(all([k in error_rates for k in S.alphabet.letters]))
    if insert_dist is None:
        letters, L = S.alphabet.letters, S.alphabet.length
        insert_dist = {k:1.0/L for k in letters}
    T = ''
    k = 0
    transcript = ''
    while k < S.length:
        if gap_open is not None:
            assert(gap_open < 1) # if not none, gap_open is assumed to be the gap probability
            assert(gap_continue < 1)
            if random.randint(0, 1000) < gap_open * 1000:
                length = np.random.geometric(1 - gap_continue)
                if random.choice([0,1]):
                    # deletion
                    transcript += 'D' * length
                    k += length
                    continue
                else:
                    # insertion
                    transcript += 'I' * length
                    T += rawrand(length, insert_dist)
                    k += 1
                    continue
        T += rawrand(1, error_rates[S[k]])[0]
        transcript += 'M' if T[-1] == S[k] else 'S'
        k += 1
    return (Sequence(T, S.alphabet), transcript)

# TODO use homopolymeric-specific gap parameters
def randread(genome, error_rates=None, coverage=40, len_mean=6000, len_var=1000):
    """Generates a random collection of lossy reads from a given genome.

    :param genome(str): the "true" original genome.
    :param coverage (float): the expected number of times each letter in the
        sequence appears in the entire read collection.
    :param len_mean (float):  the mean of the normal distribution of read lengths.
    :param len_var (float):   the variance of the normal distribution or read lengths.
    """
    N = genome.length
    num = int(1.0*N*coverage/len_mean)
    for i in range(num):
        length = max(10, min(N-1, int(np.random.normal(len_mean, len_var))))
        start = np.random.randint(0, N-length)
        x = Sequence(''.join([genome[k] for k in range(start,start+length)]), genome.alphabet)
        read,_ = mutate(x, gap_open=0.1, error_rates=error_rates)
        yield read
