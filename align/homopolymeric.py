from math import log10, floor
from . import seq, pw, hp_tokenize

class HpCondenser(object):
    """Transforms a sequence back and forth to an alternative alphabet by
    collapsing all homopolymeric substrings into single "letters".

    Args:
        alphabet (seq.Alphabet): The source alphabet.
        maxlen   (int): Maximum length of homopolymeric substrings. Longer
            hompolymeric substrings are considered to have this length.

    Attributes:
        src_alphabet (seq.Alphabet): The source alphabet.
        dst_alphabet (seq.Alphabet): The destination (condensed) alphabet.
        maxlen (int): All homopolymeric sequences longer than this are
            treated as if their length was maxlen.
        letlen (int):
            maxlen is used to decide the length of letters in the new
            alphabet (they have to be constant for all letters).

    Note:
        All operations are prefixed by ``condense_`` or ``expand_`` where
        the former means translating *to* the condensed alphabet world and the
        latter means translating *from* the condensed alphabet world.
    """
    def __init__(self, alphabet, maxlen=9):
        assert maxlen > 0
        self.letlen = int(floor(log10(maxlen))) + 2
        self.maxlen = int(maxlen)
        self.src_alphabet = alphabet
        # build the destination alphabet
        letters = []
        for char in alphabet.letters:
            for num in range(1, self.maxlen + 1):
                letter = char + str(min(num, self.maxlen)).rjust(self.letlen - 1, '0')
                letters += [letter]
        self.dst_alphabet = seq.Alphabet(letters)

    def condense_sequence(self, sequence):
        """Translates a given sequence into the condensed alphabet. For
        example::

            Tr = HpCondensor(seq.Alphabet('ACGT'))
            Tr.condense_sequence('AACCCCGGT') #=> 'A2C4G2'

        Args:
            sequence (seq.Sequence): The sequence to translate.

        Returns:
            seq.Sequence: The translated sequence in condensed alphabet.
        """
        assert(sequence.alphabet.letters == self.src_alphabet.letters)
        condensed = ''
        for char, num in hp_tokenize(str(sequence)):
            num = str(min(num, self.maxlen)).rjust(self.letlen - 1, '0')
            condensed += char + num
        return seq.Sequence(condensed, self.dst_alphabet)

    def expand_sequence(self, sequence):
        """The inverse of condense_sequence(). For exmaple::

            Tr = HpCondensor(seq.Alphabet('ACGT'))
            Tr.expand_sequence('A2C4G2') #=> 'AACCCCGGT'

        Note:
            If the original sequence contains homopolymeric substrings longer
            than self.maxlen then ``expand(condense(.))`` is not identity.

        Args:
            string (str): The condensed sequence.

        Returns:
            str: The original sequence in source alphabet.
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
        """Expands a given transcript for condensed versions of S and T to the
        equivalent transcript for S and T.

        Args:
            S (seq.Sequence): "From" sequence of the transcript.
            T (seq.Sequence): "To" sequence of the transcript.
            transcript (pw.Transcript): The transcript for condensed
                sequences.

        Returns:
            pw.Transcript: The equivalent transcript for original sequences.

        Note:
            Although ``expand(condense())`` can be lossy for homopolymeric
            substrings longer than :attr:`maxlen`, ``expand_transcript()`` does not
            have an issue with them since the original sequences (i.e ``S`` and
            ``T``) are available.
        """
        S, T = str(S), str(T)
        opseq = ''
        idx_S = transcript.idx_S * self.dst_alphabet.letter_length
        idx_T = transcript.idx_T * self.dst_alphabet.letter_length

        tokens_S = hp_tokenize(S[idx_S:])
        tokens_T = hp_tokenize(T[idx_T:])

        char_S, num_S = tokens_S.next()
        char_T, num_T = tokens_T.next()
        for op in transcript.opseq:
            if (None, None) in [(char_S, char_T), (num_S, num_T)]:
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

        return pw.Transcript(idx_S=idx_S, idx_T=idx_T,
            score=transcript.score, opseq=opseq)

    def _condense_subst_scores(self, subst_scores, **kw):
        """Helper method for condense_align_params."""
        alphabet = kw['alphabet']
        go_score, ge_score = kw['go_score'], kw['ge_score']
        hp_go_score, hp_ge_score = kw['hp_go_score'], kw['hp_ge_score']

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

    def condense_align_params(self, align_params, hp_go_score=0, hp_ge_score=0):
        """Translates alignment parameters to one that applies to the condensed
        alphabet. Translation is done based on homopolymeric indel scores which
        are treated separately from ordinary indels.

        Args:
            align_params (pw.AlignParams): Alignment parameters for the
                source alphabet.
            hp_go_score (float): Alignment score (in source alphabet) for
                homopolymeric gap open. Use 0 for linear gap penalty for
                homopolymeric indels.
            hp_ge_score (float): Alignment score (in source alphabet) for
                hompolymeric gap extension.

        Returns:
            pw.AlignParams: Alignment parameters for the destination
                alphabet.
        """
        assert(align_params.alphabet.letters == self.src_alphabet.letters)
        subst_scores_d = self._condense_subst_scores(
            align_params.subst_scores,
            alphabet=align_params.alphabet,
            hp_go_score=hp_go_score,
            hp_ge_score=hp_ge_score,
            go_score=align_params.gap_open_score,
            ge_score=align_params.gap_extend_score
        )
        return pw.AlignParams(
            alphabet=self.dst_alphabet,
            subst_scores=subst_scores_d,
            go_score=align_params.gap_open_score,
            ge_score=align_params.gap_extend_score,
            max_diversion=align_params.max_diversion
        )
