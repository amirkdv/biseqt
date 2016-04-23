"""Provides various pairwise sequence alignment algorithms. The following class
constants are inherited from ``libalign`` and are used to indicate the type of
alignment problem:

- ``GLOBAL``: Global alignment problem, i.e Needeman-Wunsch.
- ``LOCAL``: Local alignment problem, i.e Smith-Waterman.
- ``OVERLAP``: Find a suffix-prefix alignment in any direction; this includes
  alignments where a prefix of either sequence matches a suffix of the other
  and alignments where one sequence is a substring of the other.
- ``START_ANCHORED``: Find a local alignment demanding that it begins at the
  start of frame of both sequences.
- ``END_ANCHORED``: Find a local alignment demanding that it ends at the end
  of frame of both sequences.
- ``START_ANCHORED_OVERLAP``: Find a suffix-prefix alignment demanding that it
  begins at the start of frame of both sequences.
- ``END_ANCHORED_OVERLAP``:  Find a suffix-prefix alignment demanding that it
  ends at the end of frame of both sequences.
"""
from math import log
import re
import sys
from termcolor import colored
from contextlib import contextmanager
from matplotlib import pyplot as plt
from . import ffi, lib, seq, CffiObject

# alignment modes
STD_MODE = lib.STD_MODE
BANDED_MODE = lib.BANDED_MODE

# standard alignment types:
GLOBAL = lib.GLOBAL
LOCAL = lib.LOCAL
START_ANCHORED = lib.START_ANCHORED
END_ANCHORED = lib.END_ANCHORED
OVERLAP = lib.OVERLAP
START_ANCHORED_OVERLAP = lib.START_ANCHORED_OVERLAP
END_ANCHORED_OVERLAP = lib.END_ANCHORED_OVERLAP

# banded alignment types:
B_GLOBAL = lib.GLOBAL
B_OVERLAP = lib.B_OVERLAP

def hp_tokenize(string):
    """Generates (yields) homopolymeric stretches of the given sequences in
    order in tuples of the form ``(char, num, pos)``. For example::

        hp_tokenize('AAACCG') #=> [('A', 3, 0), ('C', 2, 3), ('G', 1, 5)]
    """
    for match in re.finditer(r'(.)\1*', string):
        match, pos = match.group(0), match.start()
        yield match[0], len(match), pos


class AlignScores(CffiObject):
    """Wraps the C struct ``alnscores``, see ``pwlib.h``.

    Attributes:
        alphabet (Alphabet): alphabet used to represent the sequences.
        c_obj (cffi.cdata): points to the underlying C ``alnscores`` struct.
        subst_scores (List[List]): rebuilt from the C ``double**`` upon access.

    Additionally all struct members of ``alnscores`` can be read (and not
    written) as usual attributes.
    """
    def __init__(self, alphabet=None, subst_scores=[],
                 go_score=0, ge_score=0, content_dependent_gap_scores=None):
        self.alphabet = alphabet
        # each row in the substitution matrix must be "owned" by an object that
        # is kept alive.
        L = len(self.alphabet)
        self._c_subst_rows_ka = [
            ffi.new('double[]', subst_scores[i]) for i in range(L)
        ]
        self._c_subst_full_ka = ffi.new('double *[]', self._c_subst_rows_ka)
        kw = {
            'subst_scores': self._c_subst_full_ka,
            'gap_open_score': go_score,
            'gap_extend_score': ge_score,
        }
        if content_dependent_gap_scores:
            self._c_gaps_ka = ffi.new('double []', content_dependent_gap_scores)
            kw['content_dependent_gap_scores'] = self._c_gaps_ka;
        else:
            self._c_gaps_ka = ffi.NULL
        self.c_obj = ffi.new('alnscores*', kw)

    @classmethod
    def subst_scores_from_probs(cls, alphabet, **kw):
        """Converts a substitution probability matrix to a substitution score
        matrix using a null-hypothesis letters distribution. The scores are
        natural logs of odds ratios:

            :math:`S(a_i,a_j) = \log[(1-g)\Pr(a_j|a_i)] - \log[\Pr(a_j)]`

        where :math:`S(a_i,a_j)` is the substitution score of letter
        :math:`a_i` to letter :math:`a_j` and :math:`g` is the gap probability.

        Args:
            alphabet(seq.Alphabet): the underlying alphabet, needed since
                probabilities are in order of letter index in alphabet.

        Keyword Args:
            subst_probs(List[List[float]]): As in
                :func:`biseqt.seq.Sequence.mutate`.
            letter_dist(Optional[List[float]]): Probability distributions of
                each letter of the alphabet in the null (random) hypothesis,
                default is uniform.

        Returns:
            List[List[float]]: Substitution score matrix for given alphabet,
                as expected by :func:`AlignScores`.

        Note:
            Only a linear gap model is supported for translating probabilities
            to scores. This is because under an affine model where the
            probability of opening a gap differs from that of extending a gap,
            the substitution probabilities also become context-dependent (see
            where :math:`g` appears in the formula above) which is not
            supported by ``libalign``.
        """
        L = len(alphabet)
        subst_probs = kw['subst_probs']
        letter_dist = kw.get('letter_dist', [1.0/L for k in range(L)])
        gap_prob = kw.get('gap_prob', 0)
        subst_scores = [[0 for _ in range(L)] for _ in range(L)]
        for i in range(L):
            assert(abs(1 - sum([subst_probs[i][j] for j in range(L)])) < 0.001)
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
              single indel following a substitution/match or an indel of a
              different kind.
            * The gap extend probability :math:`g_e` is the probability of a
              single indel following an indel of the same kind.

        Note:
            In the above sense the score (log likelihood) of a gap of
            length :math:`n \\ge 1` is :math:`\log g_o + (n-1)\log g_e`.
            This differs by the 1 offset from textbook definitions of the
            affine gap penalty (and from what ``libalign`` expects). The two
            are equivalent since the above gap penalty function can be
            rewritten as :math:`\log {g_o \over g_e} + n \log g_e`.
            These are precisely the scores this function returns.

            Consequently, to ensure that the gap open score is not positive,
            we require that the gap open probability be less than the gap
            extend probability.

            To enforce a linear gap model, naturally, pass identical values for
            ``go_prob`` and ``ge_prob``. This translates to a gap open
            score of 0.

        """
        return log(go_prob/ge_prob), log(ge_prob)

    def __getattr__(self, name):
        if name == 'subst_scores':
            idx = range(len(self.alphabet))
            return [[self.c_obj.subst_scores[i][j] for j in idx] for i in idx]
        elif name == 'gap_scores':
            return (self.c_obj.gap_open_score, self.c_obj.gap_extend_score)
        else:
            return super(AlignScores, self).__getattr__(name)

    def score(self, S, T, opseq, S_min_idx=0, T_min_idx=0):
        """Calculates the score for an arbitray opseq over given sequences.
        Opseqs are allowed to be partial alignments (i.e finishing before
        reaching the end of frame)::

            C = AlignScores(...)
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
        for op, num, _ in hp_tokenize(opseq):
            if op in 'MS':
                for k in range(num):
                    S_let_idx = S.c_idxseq[i + k]
                    T_let_idx = T.c_idxseq[j + k]
                    score += self.subst_scores[S_let_idx][T_let_idx]
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


class AlignFrame(CffiObject):
    """FIXME
    Args:
        S (seq.Sequence): The "from" sequence.
        T (seq.Sequence): The "to" sequence.

    Keyword Args:
        S_min_idx (int): Starting position of frame for ``S``, default is 0.
        T_min_idx (int): Starting position of frame for ``T``, default is 0.
        S_max_idx (int): Ending position (non-inclusive) of the frame for
            ``S``, default is the length of ``S``.
        T_max_idx (int): Ending position (non-inclusive) of the frame for
            ``T``, default is the length of ``T``.

    """
    def __init__(self, S, T, **kw):
        assert(isinstance(S, seq.Sequence) and isinstance(T, seq.Sequence))
        S_min_idx = kw.pop('S_min_idx', 0)
        T_min_idx = kw.pop('T_min_idx', 0)
        S_max_idx = kw.pop('S_max_idx', S.length)
        T_max_idx = kw.pop('T_max_idx', T.length)
        assert(S_max_idx <= S.length and T_max_idx <= T.length)

        self.S, self.T, = S, T
        self.c_obj = ffi.new('alnframe*', {
            'S': S.c_idxseq,
            'T': T.c_idxseq,
            'S_range': (S_min_idx, S_max_idx),
            'T_range': (T_min_idx, T_max_idx),
        })


class AlignTable(CffiObject):
    """Wraps the C struct ``dptable`` and provides a context manager to
    solve and potentially traceback an alignment problem. The corresponding
    ``alnprob`` or ``banded_alnprob`` is contained in this class as well.
    Example::

        A = seq.Alphabet('ACGT')
        S, T = A.randseq(100), A.randseq(100)
        C = biseqt.AlignScores(
            ... # snip
        )
        with biseqt.AlignTable(S, T, C, align_type=biseqt.GLOBAL) as P:
            score = P.solve()
            transcript = P.traceback()

        transcript.pretty_print(S, T, sys.stdout)

    All arguments (keyword and not) become attributes with identical names.

    Args:
        frame (pw.AlignFrame): Alignment frame.
        scores (pw.AlignScores): Alignment parameters.

    Keyword Args:
        alntype: FIXME

    Attributes:
        dp_table (List[List[float]]): The dynamic programming table built upon
            access from the underlying C ``double **``.
        c_obj (cffi.cdata): points to the underlying ``alndef`` struct.
    """
    def __init__(self, frame, scores, **kw):
        alnmode = kw['alnmode']
        alntype = kw['alntype']
        assert(isinstance(frame, AlignFrame))
        assert(alnmode in [STD_MODE, BANDED_MODE])
        alnprob_args = {
            'frame': frame.c_obj,
            'scores': scores.c_obj,
            'mode': alnmode,
            'max_new_mins': kw.get('max_new_mins', -1),
        }
        if alnmode == STD_MODE:
            assert(alntype in [GLOBAL, LOCAL, START_ANCHORED, END_ANCHORED,
                OVERLAP, START_ANCHORED_OVERLAP, END_ANCHORED_OVERLAP])
            self.c_alnparams = ffi.new('std_alnparams*', {'type': alntype});
            alnprob_args['std_params'] = self.c_alnparams
        elif alnmode == BANDED_MODE:
            assert(alntype in [B_GLOBAL, B_OVERLAP])
            self.c_alnparams = ffi.new('banded_alnparams*', {
                'type': alntype, 'dmin': kw['dmin'], 'dmax': kw['dmax']
            });
            alnprob_args['banded_params'] = self.c_alnparams

        self.c_alnprob = ffi.new('alnprob*', alnprob_args)
        self.c_obj = ffi.new('dptable*', {
            'prob': self.c_alnprob,
            'cells': ffi.NULL,
            'num_rows': -1,
            'row_lens': ffi.NULL,
        })

    def __enter__(self):
        if lib.dptable_init(self.c_obj) == -1:
            raise Exception('Failed to initialize the DP table.')
        return self

    def __exit__(self, *args):
        lib.dptable_free(self.c_obj)

    def __getattr__(self, name):
        # TODO look up attributes from self.c_obj or self.aln_prob
        return super(AlignTable, self).__getattr__(name)

    def rasterplot(self, path, transcript, fullview=False):
        ts, ss, cs = [transcript.T_idx], [transcript.S_idx], ['k']
        nums = {'M': 0, 'S': 0, '-': 0}
        colormap = {'M': 'g', 'S': 'y', '-': 'r'}
        for op in transcript.opseq:
            ts += [ts[-1] + 1 if op in 'MSI' else ts[-1]]
            ss += [ss[-1] + 1 if op in 'MSD' else ss[-1]]
            cs += colormap[op if op in 'MS' else '-']
            nums[op if op in 'MS' else '-'] += 1
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1, 1, 1)
        ax.scatter(ts, ss, color=cs, s=40, alpha=0.7, edgecolors='none')
        ax.set_aspect('equal')
        ax.grid(True)
        if fullview:
            ax.set_xlim(0, max(ts))
            ax.set_ylim(0, max(ss))
        else:
            ax.set_xlim(min(ts), max(ts))
            ax.set_ylim(min(ss), max(ss))
        # include ticks for endpoints
        ax.set_xticks(list(ax.get_xticks())[1:-1] + [min(ts), max(ts)])
        ax.set_yticks(list(ax.get_yticks())[1:-1] + [min(ss), max(ss)])
        ax.set_xlabel('T')
        ax.set_ylabel('S')
        ax.invert_yaxis()
        ax.xaxis.set_tick_params(labeltop='on')
        ax.xaxis.set_label_position('top')
        # inset plot: op stats
        width = 0.07
        inset = fig.add_axes([0.9 - 4*width, .65, 4*width, 0.2], frameon=False)
        ind = [width*i for i in range(2)]
        ind = [width * i * 1.3 for i in range(3)]
        inset.bar(ind, [nums[i]*1./len(transcript.opseq) for i in 'MS-'], width, color=[colormap[i] for i in 'MS-'])
        inset.set_aspect('equal')
        inset.set_xticks([i+width/2. for i in ind])
        inset.set_xticklabels(['M', 'S', '-'])
        inset.set_yticks([i * .2 for i in range(1,6)])
        fig.savefig(path)

    def score(self, opseq):
        """Calculates the score for an arbitray opseq. Opseqs are allowed to be
        partial alignments (i.e finishing before reaching the end of frame).::

            P = AlignProblem(...)
            P.score('MMMSSISSD') #=> 23.50
        """
        return self.scores.score(
            self.frame.S, self.frame.T, opseq, self.frame.S_range.i, self.frame.T_range.i
        )

    def solve(self):
        """Wraps :c:func:`stdpw.solve() <solve()>`: populates the DP table
        and returns the optimal score.

        Returns:
            float|NoneType: Optimal alignment score or ``None`` if no alignment
                found.
        """
        self.opt = lib.dptable_solve(self.c_obj)
        if self.opt.i == -1 or self.opt.j == -1:
            self.opt = None
            return None
        score = self.c_obj.cells[self.opt.i][self.opt.j].choices[0].score
        return score

    def traceback(self):
        """Traces back any optimal alignment found via ``solve()``.

        Returns:
            biseqt.pw.Transcript: The transcript corresponding to the alignment.
        """
        if self.opt is None:
            return None

        transcript = lib.dptable_traceback(self.c_obj, self.opt)
        if transcript == ffi.NULL:
            return None
        return Transcript(c_obj=transcript)

class Transcript(CffiObject):
    """Wrapps alignment transcripts represented as C `transcript*`.
    All keyword arguments become attributes with identical names.

    Keyword Args:
        S_idx (int): The starting position in the "from" sequence.
        T_idx (int): The starting position in the "to" sequence.
        score (float): The score of the alignment.
        opseq (str): The sequence of edit "ops" defined as follows where
            insertion/deletions are meant to mean *from S to T*::

                M match
                S substitution
                I insert
                D delete

        c_obj (Optional[cffi.cdata]): If provided all other arguments are
            ignored and instead this is used as the underlying C
            ``transcript *``.
    """
    def __init__(self, **kw):
        if 'c_obj' in kw:
            self.c_obj = kw['c_obj']
            self.c_opseq = self.c_obj.opseq
        else:
            self.c_opseq = ffi.new('char[]', kw['opseq'])
            self.c_obj = ffi.new('transcript*', {
                'S_idx': kw['S_idx'],
                'T_idx': kw['T_idx'],
                'score': kw['score'],
                'opseq': self.c_opseq,
            })

    @classmethod
    def parse_transcript(cls, raw_transcript):
        """Parses a raw transcript in string form into a :class:`Transcript`
        object. The format of raw_transcript is
        ``(<S_idx,T_idx>),<score>:<opseq>``.

        Args:
            raw_transcript (str): The raw transcript.

        Returns:
            Transcript: The populated transcript object.
        """
        assert(
            re.match('\([0-9]+,[0-9]+\),[0-9-\.]+:[MISD]+', raw_transcript)
            is not None
        )
        infostr, opseq = raw_transcript.split(':', 1)
        indices, score = infostr.rsplit(',', 1)
        # skip the open/close parens.
        S_idx, T_idx = indices[1:-1].split(',')
        kw = {
            'S_idx': int(S_idx),
            'T_idx': int(T_idx),
            'score': float(score),
            'opseq': opseq,
        }
        return cls(**kw)

    def __getattr__(self, name):
        if name == 'opseq':
            length = lib.strlen(self.c_opseq)
            return ''.join([self.c_opseq[i] for i in range(length)])
        else:
            return super(Transcript, self).__getattr__(name)

    def __repr__(self):
        return '(%d,%d),%.2f:%s' \
            % (self.S_idx, self.T_idx, self.score, self.opseq)

    def pretty_print(self, S, T, f=sys.stdout, width=120, margin=20,
                     colors=True):
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
        assert(S.alphabet.letters == T.alphabet.letters)
        S_idx, T_idx = self.S_idx, self.T_idx
        letlen = len(S.alphabet.letters[0])

        slines = tlines = []
        sline = tline = ''

        def print_lines(sline, tline, f):
            maxlen = max(len(sline), len(tline))
            sline, tline = sline.rjust(maxlen), tline.rjust(maxlen)
            f.write('%s\n%s\n' % (sline, tline))

        def new_line(sline, tline, _S_idx, _T_idx, f):
            print_lines(sline, tline, f)
            sline, tline = 'S[%d]: ' % _S_idx, 'T[%d]: ' % _T_idx
            return (max(len(sline), len(tline)), sline, tline)

        # The pre margin:
        pre_margin = min(margin, max(S_idx, T_idx) * letlen)
        sline = 'S[%d]: ' % S_idx
        tline = 'T[%d]: ' % T_idx
        counter = max(len(sline), len(tline))
        for i in reversed(range(1, pre_margin)):
            if counter >= width:
                counter, sline, tline = new_line(
                    sline, tline, S_idx+i, T_idx+i, f
                )
            sline += S[S_idx-i] if i <= S_idx else ' ' * letlen
            tline += T[T_idx-i] if i <= T_idx else ' ' * letlen
            counter += letlen

        gap = '-' * letlen
        # The alignment itself:
        for op in self.opseq:
            if counter >= width:
                counter, sline, tline = new_line(sline, tline, S_idx, T_idx, f)
            if op in 'MS':
                s, t = S[S_idx], T[T_idx]
                S_idx += 1
                T_idx += 1
            elif op == 'I':
                s, t = gap, T[T_idx]
                T_idx += 1
            elif op == 'D':
                s, t = S[S_idx], gap
                S_idx += 1
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
        post_margin = min(
            margin,
            max(
                (S.length - S_idx) * letlen,
                (T.length - T_idx) * letlen
            )
        )
        for i in range(post_margin):
            if counter >= width:
                counter, sline, tline = new_line(
                    sline, tline, S_idx + i, T_idx + i, f)
            sline += S[S_idx+i] if S_idx + i < S.length else ' ' * letlen
            tline += T[T_idx+i] if T_idx + i < T.length else ' ' * letlen
            counter += letlen

        print_lines(sline, tline, f)

class Segment(object):
    """Wraps a C ``segment``: represents an aligned pair of substrings in
    two sequences.

    Attributes:
        S_id (int): The id of the "from" sequence as found in ``seq``.
        T_id (int): The id of the "to" sequence as found in ``seq``.
        tx (biseqt.pw.Transcript): The alignment transctipt.
    """
    def __init__(self, S_id, T_id, tx):
        self.S_id = S_id
        self.T_id = T_id
        self.tx = tx

    def __repr__(self):
        return 'Segment(S_id=%d,T_id=%d,tx=%s)' \
            % (self.S_id, self.T_id, self.tx)
