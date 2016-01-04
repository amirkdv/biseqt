from math import log10, floor
from . import seq, pw, hp_tokenize
from bisect import bisect_left, bisect_right

class HpCondensedSequence(seq.Sequence):
    def __init__(self, string, alphabet, hp_positions):
        self.hp_positions = hp_positions
        super(HpCondensedSequence, self).__init__(string, alphabet)

class HpCondenser(object):
    """Transforms a sequence back and forth to an alternative alphabet by
    collapsing all homopolymeric substrings into single "letters".

    Args:
        src_alphabet (seq.Alphabet): The source alphabet.
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
    def __init__(self, src_alphabet, maxlen=9):
        assert maxlen > 0
        self.letlen = int(floor(log10(maxlen))) + 2
        self.maxlen = int(maxlen)
        self.src_alphabet = src_alphabet
        # build the destination alphabet
        letters = []
        for char in src_alphabet.letters:
            for num in range(1, self.maxlen + 1):
                num = str(min(num, self.maxlen)).rjust(self.letlen - 1, '0')
                letter = char + num
                letters += [letter]
        self.dst_alphabet = seq.Alphabet(letters)
        # Calculate the number of parititons of integers up to self.maxlen,
        # see condense_subst_probs
        partitions = [set() for _ in range(self.maxlen + 1)]
        partitions[1] = set([(1,)])
        for n in range(1, self.maxlen + 1):
            partitions[n].add((n,))
            for x in range(1, n):
                for part in partitions[n - x]:
                    partitions[n].add(tuple(sorted((x, ) + part)))

        self._num_partitions = [
            len(partitions[n]) for n in range(self.maxlen + 1)
        ]

    def condense_sequence(self, sequence):
        """Translates a given sequence into the condensed alphabet.
        For example::

            Tr = HpCondensor(seq.Alphabet('ACGT'))
            Tr.condense_sequence('AACCCCGGT') #=> 'A2C4G2'

        Args:
            sequence (seq.Sequence): The sequence to translate.

        Returns:
            seq.Sequence: The translated sequence in condensed alphabet.
        """
        assert(sequence.alphabet.letters == self.src_alphabet.letters)
        condensed = []
        hp_positions = []
        for char, num, pos in hp_tokenize(str(sequence)):
            condensed += [self._condense_hp_stretch(char, num)]
            hp_positions += [pos]
        hp_positions += [len(sequence)]
        return HpCondensedSequence(condensed, self.dst_alphabet, hp_positions)

    def _condense_hp_stretch(self, char, num):
        return char + str(min(num, self.maxlen)).rjust(self.letlen - 1, '0')

    def condense_align_params(self, align_params, hp_gap_score=0):
        """Translates alignment parameters to one that applies to the condensed
        alphabet. Translation is done based on homopolymeric indel scores which
        are treated separately from ordinary indels. Letting :math:`x,y` denote
        original alphabet letters and :math:`x_i,y_j` denote condensed alphabet
        letters, the substitution cost is the following when :math:`x=y`:

            :math:`S(x_i \\rightarrow x_j) = \\min(i,j) S(x \\rightarrow x) + G_h(|i-j|)`

        where :math:`G_h(n) = ng_h` is the linear homopolymeric gap penalty.
        The substitution score is the following when :math:`x \\ne y`:

            :math:`S(x_i \\rightarrow y_j) = \\min(i,j) S(x \\rightarrow y) + G(|i-j|)`

        where :math:`G(\cdot)` is the usual (non-homopolymeric) gap penalty.

        The returned alignment parameters also use *content-dependent gap
        scores* provided by ``libalign`` to express the fact that, for example,
        an insertion of ``A5`` is more costly than an insertion of ``A1``.

        Args:
            align_params (pw.AlignParams): Alignment parameters for the
                source alphabet.
            hp_gap_score (float): Alignment score (in source alphabet) for
                homopolymeric gaps. This is incorporated as a linear gap
                extension score in the substitution scores of homopolymeric
                stretches (see above).

        Returns:
            pw.AlignParams: Alignment parameters for the destination alphabet.
        """
        assert(align_params.alphabet.letters == self.src_alphabet.letters)

        # subst_scores gets rebuilt everytime it's accessed, fetch it once:
        subst_scores = align_params.subst_scores
        L = len(self.dst_alphabet)
        subst_scores_d = [[None for _ in range(L)] for _ in range(L)]
        gap_scores_d = [None for _ in range(L)]
        for i in range(L):
            let = self.dst_alphabet.letters[i]
            ci, ni = let[0], int(let[1:])
            gap_scores_d[i] = ni * align_params.gap_extend_score
            for j in range(L):
                let = self.dst_alphabet.letters[j]
                cj, nj = let[0], int(let[1:])
                ki = align_params.alphabet.letters.index(ci)
                kj = align_params.alphabet.letters.index(cj)
                if ci == cj:
                    subst_scores_d[i][j] = min(ni, nj) * subst_scores[ki][kj] \
                        + hp_gap_score * abs(ni-nj)
                else:
                    subst_scores_d[i][j] = min(ni, nj) * subst_scores[ki][kj]
                    if ni != nj:
                        subst_scores_d[i][j] += align_params.gap_open_score \
                            + align_params.gap_extend_score * abs(ni - nj)

        return pw.AlignParams(
            alphabet=self.dst_alphabet,
            subst_scores=subst_scores_d,
            go_score=align_params.gap_open_score,
            content_dependent_gap_scores=gap_scores_d,
            max_diversion=align_params.max_diversion
        )

    # TODO does this apply to arbitrary segments or does it require them
    # to be exactly matching (i.e seeds).
    def condense_seed(self, S, T, seed):
        """Transform the given seed to those that apply to the condensed
        sequences. This requires updating the coordinates (start position) and
        the transcript. The transcripts are only modified in their length (the
        score is left as is and note that the entire opseq is necessarily M's).

        Before translation, the seed may be shortened, if necessary, on both
        ends to ensure that its beginning and end land on the boundary of
        homopolymeric stretches (otherwise they are meaningless in the condensed
        alphabet).

        Args:
            S (HpCondensedSequence): The "from" sequence in the condensed
                alphabet.
            T (HpCondensedSequence): The "to" sequence in the condensed
                alphabet.
            seed (pw.Segment): The seed to be translated.
        """
        S_idx = bisect_right(S.hp_positions, seed.tx.S_idx)
        T_idx = bisect_right(T.hp_positions, seed.tx.T_idx)
        assert(S_idx < len(S.hp_positions) and T_idx < len(T.hp_positions))
        # do they start at an h.p. boundary?
        S_on_hp = S.hp_positions[S_idx-1] == seed.tx.S_idx
        T_on_hp = T.hp_positions[T_idx-1] == seed.tx.T_idx
        if S_on_hp and T_on_hp and S_idx*T_idx > 0:
            # S_idx and T_idx are at the start of the next h.p. stretch,
            # if both start on h.p. we can step back:
            S_idx -= 1
            T_idx -= 1

        S_end = bisect_left(S.hp_positions, seed.tx.S_idx + len(seed.tx.opseq))
        T_end = bisect_left(T.hp_positions, seed.tx.T_idx + len(seed.tx.opseq))
        assert(S_end < len(S.hp_positions) and T_end < len(T.hp_positions))

        # do they end at an h.p. boundary?
        S_on_hp = S.hp_positions[S_end] == seed.tx.S_idx + len(seed.tx.opseq)
        T_on_hp = T.hp_positions[T_end] == seed.tx.T_idx + len(seed.tx.opseq)
        if not S_on_hp or not T_on_hp:
            S_end -= 1
            T_end -= 1

        assert(S_end - S_idx == T_end - T_idx)
        length = S_end - S_idx
        if not length:
            return None
        return pw.Segment(S_id=seed.S_id, T_id=seed.T_id,
            tx=pw.Transcript(
                S_idx=S_idx, T_idx=T_idx, score=seed.tx.score, opseq=length*'M'
            )
        )
