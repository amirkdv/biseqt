from math import log10, floor
from . import seq, align, hp_tokenize

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
    def __init__(self, alphabet, maxlen=9):
        assert maxlen > 0
        self.letlen = int(floor(log10(maxlen))) + 2
        self.maxlen = int(maxlen)
        self.src_alphabet = alphabet
        self.dst_alphabet = self.translate_alphabet(alphabet)

    def translate_alphabet(self, source):
        letters = []
        for char in source.letters:
            for num in range(1, self.maxlen + 1):
                letter = char + str(min(num, self.maxlen)).rjust(self.letlen - 1, '0')
                letters += [letter]
        return seq.Alphabet(letters)

    def condense(self, sequence):
        condensed = ''
        for char, num in hp_tokenize(str(sequence)):
            num = str(min(num, self.maxlen)).rjust(self.letlen - 1, '0')
            condensed += char + num
        return seq.Sequence(condensed, self.dst_alphabet)

    def expand_sequence(self, sequence):
        """The inverse of condense(). For exmaple:

            condense("A2C4G2") #=> AACCCCGGT

        Note: if the original sequence contains homopolymeric substrings longer
        than the maxlen provided to condense() the expad(condense()) is not identity.

        :param string(str): condensed sequence
        :return str: original sequence
        """
        string = str(sequence)
        assert len(string) % self.letlen == 0
        orig = ''
        for i in range(len(string)/self.letlen):
            letter = string[self.letlen*i: self.letlen*i+self.letlen]
            char, num = letter[0], int(letter[1:])
            orig +=  char * num
        return seq.Sequence(orig, self.src_alphabet)

    def expand_transcript(self, S, T, transcript):
        S, T = str(S), str(T)
        """Expands a given transcript for condensed versions of S and T to the
        equivalent op sequence for S and T.
        """
        opseq = ''
        idx_S = transcript.idx_S * self.dst_alphabet.letter_length
        idx_T = transcript.idx_T * self.dst_alphabet.letter_length

        tokens_S = hp_tokenize(S[idx_S:])
        tokens_T = hp_tokenize(T[idx_T:])

        char_S, num_S = tokens_S.next()
        char_T, num_T = tokens_T.next()
        for op in transcript.opseq:
            if (None,None) in [char_S, char_T, num_S, num_T]:
                raise ValueError('The transcript does not match the sequences')
            if op == 'B':
                opseq += 'B'
            if op == 'M':
                opseq += 'M' * min(num_S, self.maxlen)
                if num_S > self.maxlen and num_T > self.maxlen:
                    opseq += 'M' * (min(num_S, num_T) - self.maxlen)
                    if num_S > num_T:
                        opseq += 'D' * (num_S - num_T)
                    elif num_T > num_S:
                        opseq += 'I' * (num_S - num_T)
                elif num_S > self.maxlen and num_T == self.maxlen:
                    opseq += 'D' * (num_S - self.maxlen)
                elif num_T > self.maxlen and num_S == self.maxlen:
                    opseq += 'I' * (num_T - self.maxlen)
                char_S, num_S = next(tokens_S, (None,None))
                char_T, num_T = next(tokens_T, (None,None))
            if op == 'S':
                if char_S == char_T:
                    opseq += 'M' * min(num_T, num_S)
                else:
                    opseq += 'S' * min(num_T, num_S)

                if num_T > num_S:
                    opseq += 'I' * (num_T - num_S)
                elif num_S > num_T:
                    opseq += 'D' * (num_S - num_T)
                char_S, num_S = next(tokens_S, (None,None))
                char_T, num_T = next(tokens_T, (None,None))
            if op == 'I':
                opseq += 'I' * num_T
                char_T, num_T = next(tokens_T, (None,None))
            if op == 'D':
                opseq += 'D' * num_S
                char_S, num_S = next(tokens_S, (None,None))

        return align.Transcript(idx_S=idx_S, idx_T=idx_T,
            score=transcript.score, opseq=opseq)

    def translate_subst_scores(self, subst_scores, alphabet=None,
        hp_go_score=None, hp_ge_score=None, go_score=None, ge_score=None):
        L = len(self.dst_alphabet)
        scores = [[None for _ in range(L)] for _ in range(L)]
        for i in range(L):
            let = self.dst_alphabet.letters[i]
            ci, ni = let[0], int(let[1:])
            for j in range(L):
                let = self.dst_alphabet.letters[j]
                cj, nj = let[0], int(let[1:])
                ki, kj = alphabet.letters.index(ci), alphabet.letters.index(cj)
                if ci == cj:
                    scores[i][j] = min(ni, nj) * subst_scores[ki][kj] + hp_go_score + hp_ge_score * abs(ni-nj)
                else:
                    scores[i][j] = min(ni, nj) * subst_scores[ki][kj] + go_score + ge_score* abs(ni - nj)
        return scores

    def translate_align_params(self, align_params, hp_go_score=None, hp_ge_score=None):
        assert(align_params.alphabet.letters == self.src_alphabet.letters)
        subst_scores_d = self.translate_subst_scores(
            align_params.subst_scores,
            alphabet=align_params.alphabet,
            hp_go_score=hp_go_score,
            hp_ge_score=hp_ge_score,
            go_score=align_params.gap_open_score,
            ge_score=align_params.gap_extend_score
        )
        return align.AlignParams(
            alphabet=self.dst_alphabet,
            subst_scores=subst_scores_d,
            go_score=align_params.gap_open_score,
            ge_score=align_params.gap_extend_score,
            max_diversion=align_params.max_diversion
        )
