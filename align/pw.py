"""Provides various pairwise sequence alignment algorithms. The following class
constants are inherited from ``libalign`` and are used to indicate the type of
alignment problem:

- ``GLOBAL``: Global alignment problem, i.e Needeman-Wunsch.
- ``LOCAL``: Local alignment problem, i.e Smith-Waterman.
- ``OVERLAP``: Find a prefix-suffix alignment in any direction, typically used
  for assembly.
- ``START_ANCHORED``: Find a local alignment demanding that it begins at the
  starting position of both sequences.
- ``END_ANCHORED``: Find a local alignment demanding that it ends at the ending
  position of both sequences.
"""
from math import log
import re
from termcolor import colored
from contextlib import contextmanager

from . import ffi, lib, seq, CffiObject
from . import hp_tokenize

GLOBAL = lib.GLOBAL
LOCAL = lib.LOCAL
START_ANCHORED = lib.START_ANCHORED
END_ANCHORED = lib.END_ANCHORED
OVERLAP = lib.OVERLAP

class AlignParams(CffiObject):
    """Wraps the C struct ``align_params``, see ``libalign.h``.

    Attributes:
        alphabet (Alphabet): alphabet used to represent the sequences.
        c_obj (cffi.cdata): points to the underlying `align_params` struct.

    Additionally all struct members of ``align_params`` can be read (and not
    written) as usual attributes.
    """
    def __init__(self, alphabet=None, subst_scores=[],
        go_score=0, ge_score=0, max_diversion=-1):
        self.alphabet = alphabet
        # each row in the subst matrix must be "owned" by an object that's kept
        # alive.
        self._c_subst_rows_ka = [ffi.new('double[]', subst_scores[i]) for i in range(self.alphabet.length)]
        self._c_subst_full_ka = ffi.new('double *[]', self._c_subst_rows_ka)
        self.c_obj = ffi.new('align_params*', {
            'alphabet': self.alphabet.c_obj,
            'subst_scores': self._c_subst_full_ka,
            'gap_open_score': go_score,
            'gap_extend_score': ge_score,
            'max_diversion': max_diversion,
        })

    @classmethod
    def subst_scores_from_probs(cls, alphabet, **kw):
        """Converts a substitution probability matrix to a substitution score
        matrix using a null-hypothesis letters distribution. The scores are
        natural logs of odds ratios:

        :math:`S(i,j) = \log[(1-g)\Pr(a_j|a_i)] - \log[\Pr(a_j)]`

        where :math:`S(i,j)` is the substitution score of letter :math:`a_i` to
        letter :math:`a_j` and :math:`g` is the gap probability.

        Note:
            Only a linear gap model is supported since otherwise the
            substitution scores must depend on context. For example, if the gap
            extension probability is higher than gap open probability, then the
            probability of :math:`A \\rightarrow C` is lower when observed
            immediately after a gap.

        Args:
            alphabet(seq.Alphabet): the underlying alphabet, needed since
                probabilities are in order of letter index in alphabet.

        Keyword Args:
            subst_probs(List[List[float]]): As in :func:`align.seq.Sequence.mutate`.
            letter_dist(Optional[List[float]]): Probability distributions of
                each letter of the alphabet in the null (random) hypothesis,
                default is uniform.
            gap_prob(Optional[float]): The gap probability, default is 0.

        Returns:
            List[List[float]]: Substitution score matrix for given alphabet,
                as expected by :func:`AlignParams`.
        """
        L = alphabet.length
        subst_probs = kw['subst_probs']
        letter_dist = kw.get('letter_dist', [1.0/L for k in range(L)])
        gap_prob = kw.get('gap_prob', 0)
        subst_scores = [[0 for _ in range(L)] for _ in range(L)]
        for i in range(L):
            assert(abs(1-sum([subst_probs[i][j] for j in range(L)])) < 0.001)
            for j in range(L):
                assert(subst_probs[i][j] > 0)
                assert(letter_dist[i] * letter_dist[j] != 0)
                subst_scores[i][j] = log(1-gap_prob) + \
                    log(subst_probs[i][j]) - log(letter_dist[j])
        return subst_scores

    @classmethod
    def gap_scores_from_probs(cls, go_prob, ge_prob):
        """Converts gap open/extend probabilities to gap open/extend scores
        in an affine gap penalty scheme. The probabilites are taken to mean the
        following:

            * The gap open probability :math:`g_o` is the probability of a
              single indel following a substitution/match.
            * The gap extend probability :math:`g_e` is the probability of a
              single indel following an indel of the same kind.

        Note:
            In the above sense the score (log likelihood) of a gap of
            length :math:`n \\ge 1` is :math:`\log g_o + (n-1)\log g_e`.
            This differs by the 1 offset from textbook definitions of the
            affine gap penalty (and from what ``libalign`` expects). The two are
            equivalent since the above gap penalty function can be rewritten as
            :math:`\log {g_o \over g_e} + n \log g_e`. These are precisely the
            scores this function returns.

            Consequently, to ensure that the gap open score is not positive,
            we require that the gap open probability be less than the gap
            extend probability.

        """
        assert(go_prob <= ge_prob)
        return log(go_prob/ge_prob), log(ge_prob)

    def __getattr__(self, name):
        """Allow attributes to access members of the underlying ``align_params``
        struct. Additionally provides a ``subst_scores`` attribute which builds
        a python list of lists from the corresponding C data structure.
        """
        if name == 'subst_scores':
            idx = range(self.alphabet.length)
            return [[self.c_obj.subst_scores[i][j] for j in idx] for i in idx]
        else:
            return super(AlignParams, self).__getattr__(name)

    def score(self, S, T, opseq, S_min_idx=0, T_min_idx=0):
        """Calculates the score for an arbitray opseq over given sequences.
        Opseqs are allowed to be partial alignments (i.e finishing before
        reaching the end of frame)::

            C = AlignParams(...)
            C.score('ACCTT', 'AGCTTA', 'MSMMMD')

        Args:
            S (seq.Sequence): The "from" sequence of alignment.
            T (seq.Sequence): The "to" sequence of alignment.
            opseq  (str): The edit transcript of the form ``(M|S|I|D)+``.

        Keyword Args:
            S_min_idx(Optional[int]): The starting position of the opseq in
                ``S``, default is 0.
            T_min_idx(Optiona[int]): The starting position of the opseq in
                ``T``, default is 0.
        """
        score = 0.0
        i, j = S_min_idx, T_min_idx
        for op,num in hp_tokenize(opseq):
            if op in 'MS':
                for k in range(num):
                    score += self.subst_scores[S.c_idxseq[i+k]][T.c_idxseq[j+k]]
                i, j = i + num, j + num
            elif op in 'ID':
                score += self.gap_open_score + self.gap_extend_score * num
                if op == 'I':
                    j = j + num
                else:
                    i = i + num
            else:
                raise ValueError('Invalid edit operation: %c' % op)
        return score

class AlignProblem(CffiObject):
    """Wraps the C struct ``align_problem`` and provides a context manager to
    solve and potentially traceback an alignment problem. Example::

        A = seq.Alphabet('ACGT')
        S, T = A.randseq(100), A.randseq(100)
        C = align.AlignParams(
            ... # snip
        )
        with align.AlignProblem(S, T, C, align_type=align.GLOBAL) as P:
            score = P.solve()
            transcript = P.traceback()

        transcript.pretty_print(S, T, sys.stdout)

    All arguments (keyword and not) become attributes with identical names.

    Args:
        S (seq.Sequence): The "from" sequence.
        T (seq.Sequence): The "to" sequence.
        params (pw.AlignParams): Alignment parameters.

    Keyword Args:
        S_min_idx (int): Starting position of the frame for ``S``, default is 0.
        T_min_idx (int): Starting position of the frame for ``T``, default is 0.
        S_max_idx (int): Ending position (non-inclusive) of the frame for ``S``,
            default is the length of ``S``.
        T_max_idx (int): Ending position (non-inclusive) of the frame for ``T``,
            default is the length of ``T``.

    Attributes:
        c_obj (cffi.cdata): points to the underlying ``align_problem`` struct.
        c_dp_table (cffi.cdata): points to the underlying ``double **`` DP table.
    """
    def __init__(self, S, T, params, **kw):
        align_type = kw.pop('align_type', GLOBAL)
        S_min_idx, T_min_idx = kw.pop('S_min_idx', 0), kw.pop('T_min_idx', 0)
        S_max_idx = kw.pop('S_max_idx', S.length)
        T_max_idx = kw.pop('T_max_idx', T.length)
        assert(S_max_idx <= S.length)
        assert(T_max_idx <= T.length)

        self.S, self.T, self.params, self.align_type = S, T, params, align_type
        self.c_obj = ffi.new('align_problem*', {
            'S': S.c_idxseq,
            'T': T.c_idxseq,
            'S_min_idx': S_min_idx,
            'S_max_idx': S_max_idx,
            'T_min_idx': T_min_idx,
            'T_max_idx': T_max_idx,
            'type': align_type,
            'params': params.c_obj
        })

    def __enter__(self):
        self.c_dp_table = lib.init_dp_table(self.c_obj)
        if self.c_dp_table == -1:
            raise('Got -1 from align.define().')
        self.c_dp_row_cnt = self.S_max_idx - self.S_min_idx + 1
        self.c_dp_col_cnt = self.T_max_idx - self.T_min_idx + 1
        return self

    def __exit__(self, *args):
        if self.c_dp_table not in [ffi.NULL, -1]:
            lib.free_dp_table(self.c_dp_table, self.c_dp_row_cnt, self.c_dp_col_cnt)

    def score(self, opseq):
        """Calculates the score for an arbitray opseq. Opseqs are allowed to be
        partial alignments (i.e finishing before reaching the end of frame).::

            P = AlignProblem(...)
            P.score('MMMSSISSD') #=> 23.50
        """
        return self.params.score(self.S, self.T, opseq, self.S_min_idx, self.T_min_idx)

    def solve(self, print_dp_table=False):
        """Populates the DP table and returns the optimal score.

        Args:
            print_dp_table(Optional[bool]): whether or not to print the
                fully calculated DP table, default is ``False``.
        Returns:
            float|NoneType: Optimal alignment score or ``None`` if no alignment
                found.
        """
        self.opt = lib.solve(self.c_dp_table, self.c_obj)
        if print_dp_table:
            mat = self.dp_table
            for i in range(len(mat)):
                print [round(f,2) for f in mat[i]]
        if self.opt == ffi.NULL:
            self.opt = None
            return None
        score = self.opt.choices[0].score
        return score

    def traceback(self):
        """Traces back any optimal alignment found via ``solve()``.

        Returns:
            align.Transcript
        """
        if self.opt is None:
            return None

        raw_transcript = lib.traceback(self.c_dp_table, self.c_obj, self.opt)
        if raw_transcript == ffi.NULL:
            return None
        tx = Transcript(raw_transcript=ffi.string(raw_transcript))
        lib.free(raw_transcript)
        return tx

    def __getattr__(self, name):
        """Allow attributes to access members of the underlying `align_problem`
        struct. Additionally provide a `dp_table` attribute which builds a
        python list of lists from the corresponding C data structure.
        """
        if name == 'dp_table':
            i_idx = range(self.S_max_idx - self.S_min_idx + 1)
            j_idx = range(self.T_max_idx - self.T_min_idx + 1)
            def score(i,j):
                if self.c_dp_table[i][j].num_choices != 0:
                    return self.c_dp_table[i][j].choices[0].score
                else:
                    return None
            return [[score(i,j) for j in j_idx] for i in i_idx]
        else:
            return super(AlignProblem, self).__getattr__(name)

class Transcript(object):
    """A wrapper for alignment transcripts. Solutions to the alignment problem are
    represented by transcript strings with the following format::

        (<Si,Tj>),<score>:<opseq>

    ``Si`` and ``Tj`` are integers specifying the positiong along each string
    where the alignment begins (relative to the corresponding frames). ``score``
    is the score of the transcript to 2 decimal places. And ``opseq`` is a
    sequence of edit "ops" defined as follows where insertion/deletions are
    meant to mean "from S to T"::

        M match
        S substitution
        I insert
        D delete


    Args:
        raw_transcript (Optional[str]): If provided all other arguments are
            ignored and instead this string is parsed to populate the
            attributes.
    """
    def __init__(self, idx_S=0, idx_T=0, score=0.0, opseq='', raw_transcript=None):
        if raw_transcript is None:
            self.idx_S, self.idx_T = idx_S, idx_T
            self.score, self.opseq = score, opseq
            return

        assert(re.match('\([0-9]+,[0-9]+\),[0-9-\.]+:[MISD]+', raw_transcript) is not None)
        infostr, opseq = raw_transcript.split(':', 1)
        indices, score = infostr.rsplit(',', 1)
        idx_S, idx_T = indices[1:-1].split(',') # skip the open/close parens
        self.idx_S, self.idx_T = int(idx_S), int(idx_T)
        self.opseq, self.score = opseq, float(score)

    def __repr__(self):
        return '(%d,%d),%.2f:%s' % (self.idx_S, self.idx_T, self.score, self.opseq)

    def pretty_print(self, S, T, f, width=120, margin=20, colors=True):
        """Pretty prints a transcript to f.

        Args:
            S (seq.Sequence): The "from" sequence.
            T (seq.Sequence): The "to" sequence.
            f (file): Open file handle to write the output to; for standard
                output use ``sys.stdout``.
            width (Optional[int]): Terminal width used for wrapping;
                default is 120.
            margin (Optional[length]): Length of leading and trailing sequences
                in S and T before and after the alignment; default is 20.
            colors (Optional[bool]): Whether or not to use colors in output;
                default is ``True``.
        """
        assert(S.alphabet.letter_length == T.alphabet.letter_length)
        assert(S.alphabet.letters == T.alphabet.letters)
        letlen = S.alphabet.letter_length
        idx_S, idx_T = self.idx_S, self.idx_T

        slines = tlines = []
        sline = tline = ''

        def print_lines(sline, tline, f):
            maxlen = max(len(sline), len(tline))
            sline, tline = sline.rjust(maxlen), tline.rjust(maxlen)
            f.write('%s\n%s\n' % (sline,tline))

        def new_line(sline, tline, _idx_S, _idx_T, f):
            print_lines(sline, tline, f)
            sline, tline = 'S[%d]: ' % _idx_S, 'T[%d]: ' % _idx_T
            return (max(len(sline), len(tline)), sline, tline)

        # The pre margin:
        pre_margin = min(margin, max(idx_S, idx_T))
        sline = 'S[%d]: ' % idx_S
        tline = 'T[%d]: ' % idx_T
        counter = max(len(sline), len(tline))
        for i in reversed(range(1, pre_margin)):
            if counter >= width:
                counter, sline, tline = new_line(sline, tline, idx_S+i, idx_T+i, f)
            sline += S[idx_S-i] if i <= idx_S else ' '
            tline += T[idx_T-i] if i <= idx_T else ' '
            counter += letlen

        gap = '-' * S.alphabet.letter_length
        # The alignment itself:
        for op in self.opseq:
            if counter >= width:
                counter, sline, tline = new_line(sline, tline, idx_S, idx_T, f)
            if op in 'MS':
                s, t = S[idx_S], T[idx_T]
                idx_S += 1
                idx_T += 1
            elif op == 'I':
                s, t = gap, T[idx_T]
                idx_T += 1
            elif op == 'D':
                s, t = S[idx_S], gap
                idx_S += 1
            else:
                raise ValueError('Invalid edit operation: %c' % op)
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
        post_margin = min(margin, max(S.length - idx_S, T.length - idx_T))
        for i in range(post_margin):
            if counter >= width:
                counter, sline, tline = new_line(
                    sline, tline, idx_S + i, idx_T + i, f)
            sline += S[idx_S+i] if idx_S + i < S.length else ' '
            tline += T[idx_T+i] if idx_T + i < T.length else ' '
            counter += letlen

        print_lines(sline, tline, f)
