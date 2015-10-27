import numpy as np
import random
from termcolor import colored
from collections import namedtuple


from . import ffi, lib, CffiObject

"""A tuple wrapper for alignment strings. Solutions to the alignment problem are
represented by transcript strings with the following format:

    (<Si,Tj>),<score>:...

Si and Tj are integers specifying the positiong along each string where
the alignment begins. Score is the score of the transcript to 2 decimal
places. What follows the ':' is a sequence of "ops" defined as follows:
    B begin
    M match
    S substitution
    I insert
    D delete
All op sequences begin with a B and insertion/deletions are meant to
mean "from S to T".
"""
# FIXME when libalign::solve can return only scores make this in to a class
Transcript = namedtuple('Transcript', ['S_idx', 'T_idx', 'score', 'opseq'])

def parse_transcript(transcript):
    """Parses a raw transcript as given by `libalign.so` and returns an
    Transcript tuple. Raw transcripts are expected to have the following
    format:

        (<Si,Tj>),<score>:<opseq>
    """

    infostr, opseq = transcript.split(':', 1)
    indices, score = infostr.rsplit(',', 1)
    S_idx, T_idx = indices[1:-1].split(',') # skip the open/close parens
    return Transcript(S_idx=int(S_idx), T_idx=int(T_idx), opseq=opseq, score=float(score))

def print_transcript(S, T, transcript, f, width=120, margin=20, colors=True):
    """Pretty prints a given transcript to f.

    :param S: "from" sequence (anything that casts to the correct str object).
    :param T: "to" sequence (like S).
    :param transcript(Transcript): the transcript output of libalign.
    :param f: file handle to write the output to; for standard output use
        `sys.stdout`.
    :param width(optional): terminal width (int) used for wrapping; default 120.
    :param margin(optional): length (int) of leading and trailing sequences in S
        and T before and after the alignment; default 20.
    :param colors(optional): whether or not (truthy) to use colors in output;
        default True.
    """
    assert(S.alphabet.letter_length == T.alphabet.letter_length)
    assert(S.alphabet.letters == T.alphabet.letters)
    letlen = S.alphabet.letter_length
    S_idx, T_idx = transcript.S_idx, transcript.T_idx

    slines = tlines = []
    sline = tline = ''

    def print_lines(sline, tline, f):
        maxlen = max(len(sline), len(tline))
        sline, tline = sline.rjust(maxlen), tline.rjust(maxlen)
        f.write('%s\n%s\n' % (sline,tline))

    def new_line(sline, tline, _S_idx, _T_idx, f):
        print_lines(sline, tline, f)
        sline, tline = 'S[%d]: ' % _S_idx, 'T[%d]: ' % _T_idx
        return (max(len(sline), len(tline)), sline, tline)

    # The pre margin:
    pre_margin = min(margin, max(S_idx, T_idx))
    sline = 'S[%d]: ' % max(0, S_idx - pre_margin + 1)
    tline = 'T[%d]: ' % max(0, T_idx - pre_margin + 1)
    counter = max(len(sline), len(tline))
    for i in reversed(range(1, pre_margin)):
        if counter >= width:
            counter, sline, tline = new_line(sline, tline, S_idx+i, T_idx+i, f)
        sline += S[S_idx-i] if i <= S_idx else ' '
        tline += T[T_idx-i] if i <= T_idx else ' '
        counter += letlen

    # The alignment itself:
    for i,op in enumerate(transcript.opseq):
        if op == 'B':
            continue
        if counter >= width:
            counter, sline, tline = new_line(sline, tline, S_idx, T_idx, f)
        if op in 'MS':
            s, t = S[S_idx], T[T_idx]
            S_idx += 1
            T_idx += 1
        elif op == 'I':
            s, t = '-', T[T_idx]
            T_idx += 1
        elif op == 'D':
            s, t = S[S_idx], '-'
            S_idx += 1
        on_color = color = None
        if colors:
            if op in 'MS':
                color = 'green' if op == 'M' else 'red'
            elif op == 'I':
                on_color = 'on_red'
            elif op == 'D':
                on_color = 'on_red'
        sline += colored(s, color=color, on_color=on_color)
        tline += colored(t, color=color, on_color=on_color)
        counter += letlen

    # The post margin:
    post_margin = min(margin, max(S.length-S_idx, T.length-T_idx))
    for i in range(post_margin):
        if counter >= width:
            counter, sline, tline = new_line(
                sline, tline, S_idx+i, T_idx+i, f)
        sline += S[S_idx+i] if S_idx + i < S.length else ' '
        tline += T[T_idx+i] if T_idx + i < T.length else ' '
        counter += letlen

    print_lines(sline, tline, f)

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

    def randstr(self, length, letters_dist=None):
        """Generates a random sequence of the specified length from this
        alphabet. Optionally a discrete distribution may be specified.

        :param length(int): length of generated sequence.
        :param letters_dist(dict): Optional; The probability distribution of
            letters as a list of probabilities in order of letters in
            self.c_obj.letters. Default is uniform.
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
        return Sequence(self.randstr(length, letters_dist), self)

class Sequence():
    """Wraps a C char[] and its corresponding `idx_seq int[]`. Note that this
    is *not* a CffiObject subclass since there's no underlying C struct. String
    indexing and slicing are supported:

        x = seq.Sequence("A2C3A2", seq.Alphabet(["A2", "C3"]))
        x[1]   # => "C3"
        x[:2]  # => "C3A2"
        x[:-2] # => "A2"

    Attributes:
        alphabet (Alphabet)
        c_charseq (cffi.cdata): points to the underlying C char[].
        c_idxseq  (cffi.cdata): points to the actually used int*.
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
        * the current letter will be replaced by an arbitrary letter with a
            distribution that depends on the original letter.
        * with a certain fixed probability a gap may be openned with random length
            with geometric distribution.
        Accordingly, an opseq (see `align.solve()`) is generated which corresponds
        to the performed edit sequence.

        :param S(seq.Sequence): original sequence.
        :param go_prob(float): probability of a gap starting at any position.
        :param ge_prob(float): Bernoulli success probability of the gap
            extension distribution
        :param subst_probs(list[list]): the probability distribution for each
            pair of possible substitutions such that subst_probs[i][j] is the
            probability of letter i (integer index) of the alphabet being
            substituted by letter j (integer index).
        :param insert_dist(list): the distribution passed to Alphabet.randstr()
            when inserting arbitrary strings; default is uniform.
        """
        assert(go_prob < 1)
        assert(ge_prob < 1)
        T = ''
        k = 0
        opseq = 'B'
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

    def randread(self, subst_probs=None, go_prob=0, ge_prob=0, coverage=40, len_mean=6000, len_var=1000):
        """Generates a random collection of lossy reads from the current sequence.

        :param subst_probs(list[list]): as in mutate(), letter-by-letter
            probabilities of substitutions by the seqeuencing process.
        :param go_prob(float): as in mutate(), probability that at any point in any
            read a gap will be opened (use 1 for linear gap penalty)
        :param ge_prob(float): as in mutate(), probability that at any point in any
            gap in a read, the gap will be extended.
        :param coverage (float): the expected number of times each letter in the
            sequence appears in the entire read collection.
        :param len_mean (float):  the mean of the normal distribution of read lengths.
        :param len_var (float):   the variance of the normal distribution or read lengths.
        """
        N = self.length
        num = int(1.0 * N * coverage/len_mean)
        for i in range(num):
            length = max(10, min(N-1, int(np.random.normal(len_mean, len_var))))
            start = np.random.randint(0, N-length)
            x = Sequence(''.join([self[k] for k in range(start, start + length)]), self.alphabet)
            read, _ = x.mutate(go_prob=go_prob, ge_prob=ge_prob, subst_probs=subst_probs)
            yield read
