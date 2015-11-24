from math import log10, floor
from . import seq, pw, hp_tokenize, tuples


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
        condensed = ''.join(
            [self._condense(*x) for x in hp_tokenize(str(sequence))]
        )
        return seq.Sequence(condensed, self.dst_alphabet)

    def _condense(self, char, num):
        return char + str(min(num, self.maxlen)).rjust(self.letlen - 1, '0')

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
            orig += char * num
        return seq.Sequence(orig, self.src_alphabet)

    def expand_transcript(self, S, T, transcript):
        """Expands a given transcript for condensed versions of S and T to the
        equivalent transcript for S and T. The score is left untouched.

        Args:
            S (seq.Sequence): "From" sequence of the transcript in the source
                alphabet.
            T (seq.Sequence): "To" sequence of the transcript in the source
                alphabet.
            transcript (pw.Transcript): The transcript for condensed
                sequences.

        Returns:
            pw.Transcript: The equivalent transcript for original sequences.

        Note:
            Although ``expand(condense())`` can be lossy for homopolymeric
            substrings longer than :attr:`maxlen`, ``expand_transcript()`` does
            not have an issue with them since the original sequences (i.e ``S``
            and ``T``) are available.
        """
        S, T = str(S), str(T)
        opseq = ''
        tokens_S = hp_tokenize(S)
        tokens_T = hp_tokenize(T)

        char_S, num_S = tokens_S.next()
        char_T, num_T = tokens_T.next()
        # calculate the original S_idx and T_idx
        S_idx, T_idx = 0, 0
        cnt = 0
        while cnt < transcript.S_idx:
            S_idx += num_S
            cnt += 1
            char_S, num_S = tokens_S.next()
        cnt = 0
        while cnt < transcript.T_idx:
            T_idx += num_T
            cnt += 1
            char_T, num_T = tokens_T.next()

        # translate the opseq
        for op in transcript.opseq:
            if (None, None) in [(char_S, num_S), (char_T, num_T)]:
                raise ValueError('The transcript does not match the sequences')

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
                char_S, num_S = next(tokens_S, (None, None))
                char_T, num_T = next(tokens_T, (None, None))
            if op == 'S':
                if char_S == char_T:
                    opseq += 'M' * min(num_T, num_S)
                else:
                    opseq += 'S' * min(num_T, num_S)

                if num_T > num_S:
                    opseq += 'I' * (num_T - num_S)
                elif num_S > num_T:
                    opseq += 'D' * (num_S - num_T)
                char_S, num_S = next(tokens_S, (None, None))
                char_T, num_T = next(tokens_T, (None, None))
            if op == 'I':
                opseq += 'I' * num_T
                char_T, num_T = next(tokens_T, (None, None))
            if op == 'D':
                opseq += 'D' * num_S
                char_S, num_S = next(tokens_S, (None, None))

        return pw.Transcript(
            S_idx=S_idx, T_idx=T_idx, score=transcript.score, opseq=opseq
        )

    def condense_subst_probs(self, **kw):
        """Translates the substitution probabilities in the source alphabet
        to substitution probabilities in the destination (condensed) alphabet.
        Letting :math:`x,y` denote original alphabet letters and
        :math:`x_i,y_j` denote condensed alphabet letters, the translation
        formula is the following when the length of letters are identical:

            :math:`\Pr(x_i \\rightarrow y_i) = \Pr(x \\rightarrow y)^i(1-g_h)^{i-1}`

        where :math:`g_h` is the homopolymeric gap probability (only a linear
        model is supported). When the length of letters differ:

            :math:`\Pr(x_i \\rightarrow y_j) = \\pi(i)\Pr(x \\rightarrow y)^{\\min(i,j)} (1-g_h)^{\\min(i,j)-1}g_h^{|i-j|}`

        where :math:`\\pi(\\cdot)` is the integer partition function.

        Note:
            The calculations here may have serious errors. In fact, the
            calculated probabilities as described above don't necessarily add
            up to 1! Returned probability matrix is normalized in each row
            to make sure the output is not terribly wrong.
        """
        subst_probs = kw['subst_probs']
        hp_gap_prob = kw['hp_gap_prob']
        assert(hp_gap_prob > 0)
        N, L = len(self.src_alphabet), len(self.dst_alphabet)
        letters_dist = kw.get('letters_dist', [1.0/N for _ in range(N)])
        subst_probs_d = [[None for _ in range(L)] for _ in range(L)]
        for i in range(L):
            let = self.dst_alphabet.letters[i]
            ci, ni = let[0], int(let[1:])  # e.g A31 gives ci = 'A' and ni = 31
            for j in range(L):
                let = self.dst_alphabet.letters[j]
                cj, nj = let[0], int(let[1:])
                ki = self.src_alphabet.letters.index(ci)
                kj = self.src_alphabet.letters.index(cj)
                subst_probs_d[i][j] = subst_probs[ki][kj] ** min(ni, nj) * \
                    (1-hp_gap_prob) ** (min(ni, nj) - 1)
                if ni != nj:
                    subst_probs_d[i][j] *= self._num_partitions[ni] * \
                        hp_gap_prob ** abs(ni - nj)
        # FIXME probabilities don't add up to 1, normalize:
        for idx, row in enumerate(subst_probs_d):
            s = sum(row)
            subst_probs_d[idx] = [x/s for x in row]
        return subst_probs_d

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

        Args:
            align_params (pw.AlignParams): Alignment parameters for the
                source alphabet.
            hp_gap_score (float): Alignment score (in source alphabet) for
                homopolymeric gaps. This is incorporated in the substitution
                costs of homopolymeric stretches. Only a linear gap model is
                well-defined.

        Returns:
            pw.AlignParams: Alignment parameters for the destination alphabet.
        """
        assert(align_params.alphabet.letters == self.src_alphabet.letters)

        # subst_scores gets rebuilt everytime it's accessed, fetch it once:
        subst_scores = align_params.subst_scores
        L = len(self.dst_alphabet)
        subst_scores_d = [[None for _ in range(L)] for _ in range(L)]
        for i in range(L):
            let = self.dst_alphabet.letters[i]
            ci, ni = let[0], int(let[1:])
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
            ge_score=align_params.gap_extend_score,
            max_diversion=align_params.max_diversion
        )

    def condense_seed(self, S, T, seed):
        """Condenses a seed into a :class:`tuples.Segment` for the
        corresponding condensed sequences. Note that this process is not
        well-defined unless certain assumptions is made about segments. Here we
        require that the provided seed is an exactly matching seed, at least in
        some condensed alphabet.

        Args:
            S (seq.Sequence): The "from" sequence in original alphabet.
            T (seq.Sequence): The "to" sequence in original alphabet.
            seed (tuples.Segment): An exactly matching segment, at least in
                some condensed alphabet.

        Returns:
            tuples.Segment: Corresponding segment translated such that it
                applies to condensed sequences as generated by us.
        """
        assert(set(seed.tx.opseq).issubset(set('IMD')))
        tokens_S = hp_tokenize(str(S))
        tokens_T = hp_tokenize(str(T))
        # calculate the condensed S_idx and T_idx
        S_min_idx, cnt = 0, 0
        while cnt < seed.tx.S_idx:
            char_S, num_S = tokens_S.next()
            S_min_idx += 1
            cnt += num_S
        assert(cnt == seed.tx.S_idx)

        T_min_idx, cnt = 0, 0
        while cnt < seed.tx.T_idx:
            char_T, num_T = tokens_T.next()
            T_min_idx += 1
            cnt += num_T
        assert(cnt == seed.tx.T_idx)
        i = 0
        opseq = ''
        # Heavy usage of the assumption that the opseq belongs to an exactly
        # matching seed follows (hopefully it works too if the index is built
        # with a different HpCondenser)
        while i < len(seed.tx.opseq):
            char_T, num_T = tokens_T.next()
            char_S, num_S = tokens_S.next()
            assert(char_S == char_T)
            if self._condense(char_S, num_S) == self._condense(char_T, num_T):
                opseq += 'M'
            else:
                opseq += 'S'
            i += max(num_S, num_T)

        tx = pw.Transcript(
            S_idx=S_min_idx, T_idx=T_min_idx, score=seed.tx.score, opseq=opseq
        )
        return tuples.Segment(S_id=seed.S_id, T_id=seed.T_id, tx=tx)


class HpCondensedIndex(tuples.Index):
    """Adds homopolymeric-condensed indexing support to tuples indices.
    Example usage::

        A = seq.Alphabet('ACGT')
        Tr = homopolymeric.HpCondenser(A, maxlen=9)
        B = tuples.TuplesDB('path/to/file', alphabet=Tr.dst_alphabet)
        HpI = homopolymeric.HpCondensedIndex(B, 5, hp_condenser=Tr)
    """
    def __init__(self, *args, **kwargs):
        self.hp_condenser = kwargs['hp_condenser']
        super(HpCondensedIndex, self).__init__(*args)

    def tup_scan(self, string):
        """Similar to :func:`align.tuples.Index.tup_scan` except a tuple of
        length N is taken to mean N letters in the condensed alphabet.
        The yielded indices also refer to positions in the condensed
        sequence. For example::

            Tr = HpCondenser(seq.Alphabet('ACGT'), maxlen=3)
            Idx = HpCondensedIndex(Tr)
            string = 'AAACCCCGGTGGT'
            Idx.tup_scan(string, 5) # => ('A3C3G2T1G2', 0), ('C3G2T1G2T1', 1)

        This modification alone is enough to ensure
        :func:`index() <align.tuples.Index.index>` works properly in the
        condensed alphabet.
        """
        tup = []
        idx = [0]
        for char, num in hp_tokenize(string):
            if len(tup) == self.wordlen:
                yield ''.join(tup), idx[0]
                tup.pop(0)
                idx.pop(0)
            tup += [self.hp_condenser._condense(char, num)]
            idx += [(idx[-1] if idx else 0) + 1]
        if tup:
            yield ''.join(tup), idx[0]

    def seeds(self, S_id, T_id):
        """Wraps parent's :func:`seeds() <align.tuples.Index.seeds>` to
        translate all seed transcripts back to original alphabet."""
        condensed_seeds = super(HpCondensedIndex, self).seeds(S_id, T_id)
        res = []
        S, T = self.tuplesdb.loadseq(S_id), self.tuplesdb.loadseq(T_id)
        for seed in condensed_seeds:
            tx = self.hp_condenser.expand_transcript(S, T, seed.tx)
            res += [tuples.Segment(S_id=S_id, T_id=T_id, tx=tx)]
        return res
