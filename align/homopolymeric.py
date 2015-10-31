from math import log10, floor
from .seq import Alphabet

def hp_tokenize(string):
    """Generator for homopolymeric substrings in a given sequences. Each value
    is a (char, num) tuple.
    """
    counter = 0
    cur = string[0]
    while counter < len(string):
        if string[counter] == cur:
            counter += 1
        else:
            yield string[0], counter
            string = string[counter:]
            counter = 0
            cur = string[0]
    # left overs:
    if counter and len(string):
        yield string[0], counter

class HpCondensor(object):
    """Transforms a sequence back and forth to an alternative alphabet by
    collapsing all homopolymeric substrings into single "letters". For example:

        T = HpCondensor()
        T.condense("AACCCCGGT") #=> A2C4G2T1
        T.expand("A2C4G2T1")   #=> AACCCCGGT

    Attributes:
        maxlen (int): if max len is truthy, all homopolymeric sequences
            longer than maxlen are treated as if their length was maxlen.
        letlen (int):
            maxlen is used to decide the length of letters in the new
            alphabet (they have to be constant for all letters).
    """
    def __init__(self, maxlen=9):
        assert maxlen > 0
        self.letlen = int(floor(log10(maxlen))) + 2
        self.maxlen = int(maxlen)

    def translate_alphabet(self, source='ACGT'):
        for char in source:
            for num in range(1, self.maxlen + 1):
                letter = char + str(min(num, self.maxlen)).rjust(self.letlen - 1, '0')
                yield char, num, letter

    def condense(self, string):
        condensed = ''
        for char, num in hp_tokenize(string):
            num = str(min(num, self.maxlen)).rjust(self.letlen - 1, '0')
            condensed += char + num
        return condensed

    def expand(self, string):
        """The inverse of condense(). For exmaple:

            condense("A2C4G2") #=> AACCCCGGT

        Note: if the original sequence contains homopolymeric substrings longer
        than the maxlen provided to condense() the expad(condense()) is not identity.

        :param string(str): condensed sequence
        :return str: original sequence
        """
        assert len(string) % self.letlen == 0
        orig = ''
        for i in range(len(string)/self.letlen):
            letter = string[self.letlen*i: self.letlen*i+self.letlen]
            char, num = letter[0], int(letter[1:])
            orig +=  char * num
        return orig

    # TODO allow translating Transcript objects (requires translatin scores)
    def expand_opseq(self, S, T, opseq):
        S, T = str(S), str(T)
        """Expands a given sequence of edit ops (string of B/M/S/I/D) generated
        for the condensed versions of S and T to the equivalent op sequence for
        S and T.
        """
        assert all([s in 'BMISD' for s in opseq])
        orig = ''
        tokens_S = hp_tokenize(S)
        tokens_T = hp_tokenize(T)

        char_S, num_S = tokens_S.next()
        char_T, num_T = tokens_T.next()
        for op in opseq:
            if (None,None) in [char_S, char_T, num_S, num_T]:
                raise ValueError('The transcript does not match the sequences')
            if op == 'B':
                orig += 'B'
            if op == 'M':
                orig += 'M' * min(num_S, self.maxlen)
                if num_S > self.maxlen and num_T > self.maxlen:
                    orig += 'M' * (min(num_S, num_T) - self.maxlen)
                    if num_S > num_T:
                        orig += 'D' * (num_S - num_T)
                    elif num_T > num_S:
                        orig += 'I' * (num_S - num_T)
                elif num_S > self.maxlen and num_T == self.maxlen:
                    orig += 'D' * (num_S - self.maxlen)
                elif num_T > self.maxlen and num_S == self.maxlen:
                    orig += 'I' * (num_T - self.maxlen)
                char_S, num_S = next(tokens_S, (None,None))
                char_T, num_T = next(tokens_T, (None,None))
            if op == 'S':
                if char_S == char_T:
                    orig += 'M' * min(num_T, num_S)
                else:
                    orig += 'S' * min(num_T, num_S)

                if num_T > num_S:
                    orig += 'I' * (num_T - num_S)
                elif num_S > num_T:
                    orig += 'D' * (num_S - num_T)
                char_S, num_S = next(tokens_S, (None,None))
                char_T, num_T = next(tokens_T, (None,None))
            if op == 'I':
                orig += 'I' * num_T
                char_T, num_T = next(tokens_T, (None,None))
            if op == 'D':
                orig += 'D' * num_S
                char_S, num_S = next(tokens_S, (None,None))

        return orig

    def translate_subst_scores(self, subst_scores, alphabet='ACGT',
        hp_go_score=0, hp_ge_score=-0.5, go_score=-3, ge_score=-2):
        alphabet_d = [x for x in self.translate_alphabet(alphabet)]
        L = len(alphabet_d)
        scores = [[None for _ in range(L)] for _ in range(L)]
        for i in range(L):
            ci, ni, _ = alphabet_d[i]
            for j in range(L):
                cj, nj, _ = alphabet_d[j]
                ki, kj = alphabet.index(ci), alphabet.index(cj)
                if ci == cj:
                    scores[i][j] = min(ni, nj) * subst_scores[ki][kj] + hp_go_score + hp_ge_score * abs(ni-nj)
                else:
                    scores[i][j] = min(ni, nj) * subst_scores[ki][kj] + go_score + ge_score* abs(ni - nj)
        return scores
