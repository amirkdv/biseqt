from math import log
import re
from termcolor import colored

from . import ffi, lib, seq, CffiObject
from .homopolymeric import hp_tokenize

global lib
"""Values for enum `align_type`.
"""
ALIGN_GLOBAL = lib.GLOBAL
ALIGN_LOCAL = lib.LOCAL
ALIGN_START_ANCHORED = lib.START_ANCHORED
ALIGN_END_ANCHORED = lib.END_ANCHORED
ALIGN_OVERLAP = lib.OVERLAP

class AlignParams(CffiObject):
    """Wraps the C struct align_params, see `libalign.h`

    Attributes:
        alphabet (Alphabet): alphabet used to represent the sequences.
        _c_subst_rows_ka (list[cffi.cdata]): has ownership of (keeps alive) C pointers
            to rows in the substitution matrix.
        _c_subst_full_ka (cffi.cdata): has ownership of (keeps alive) the C
            pointer to the full substitution matrix.
        c_obj (cffi.cdata): points to the underlying `align_params` struct.
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
    def subst_scores_from_probs(cls, subst_probs, alphabet, letter_dist=None):
        """Converts a substitution probability matrix to a substitution score
        matrix using a null-hypothesis letters distribution. The scores are
        natural logs of odds ratios.

        :param subst_probs(list[list]): as in seq.Sequence.mutate()
        :param alphabet(Alphabet): the underlying alphabet, needed since
            probabilities are in order of letter index in alphabet.
        :param letter_dist(list[float]): probability distributions of each
            letter of the alphabet in the null (random) hypothesis.

        :return subst_scores(list[list]): as expected in AlignParams.__init__().
        """
        L = alphabet.length
        if letter_dist is None:
            letter_dist = [1.0/L for k in range(L)]
        subst_scores = [[0 for _ in range(L)] for _ in range(L)]
        for i in range(L):
            assert(abs(1-sum([subst_probs[i][j] for j in range(L)])) < 0.001)
            for j in range(L):
                assert(subst_probs[i][j] > 0)
                assert(letter_dist[i] * letter_dist[j] != 0)
                subst_scores[i][j] = log(subst_probs[i][j]) \
                    - log(letter_dist[i]) - log(letter_dist[j])
        return subst_scores

    @classmethod
    def gap_scores_from_probs(cls, go_prob, ge_prob):
        """Converts gap open/extend probabilities to gap open/extend scores
        in an affine gap penalty scheme. If go_prob = 1 (which means gap
        extension does not require a separate "opening" event) we get a linear
        gap penalty scheme.
        """
        return log(go_prob), log(ge_prob)

    def __getattr__(self, name):
        """Allow attributes to access members of the underlying `align_params`
        struct. Additionally provides a `subst_scores` attribute which builds
        a python list of lists from the corresponding C data structure.
        """
        if name == 'subst_scores':
            idx = range(self.alphabet.length)
            return [[self.c_obj.subst_scores[i][j] for j in idx] for i in idx]
        else:
            return super(AlignParams, self).__getattr__(name)

class AlignProblem(CffiObject):
    """Wraps the C struct `align_problem', see `libalign.h`

    Attributes:
        S (Sequence): "from" sequence.
        T (Sequence): "to" sequence.
        params (AlignParams): parameters for the alignment problem.
        align_type (int): any of the ALIGN_* constants defined here (pointing to
            values of C enum `align_type`).
        c_obj (cffi.cdata): points to the underlying `align_problem` struct.
        c_dp_table (cffi.cdata): points to the underlying `double **` DP table.
    """
    def __init__(self, S=None, T=None, params=None, align_type=ALIGN_GLOBAL,
        S_min_idx=0, T_min_idx=0, S_max_idx=None, T_max_idx=None):

        # protect against index overflows
        S_max_idx = S.length if S_max_idx else S.length
        T_max_idx = T.length if T_max_idx else T.length
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
        global lib

    def score(self, opseq):
        """Calculates the score for an arbitray opseq. Opseqs are allowed to be
        partial alignments (i.e finishing before reaching the end of frame).

            P = AlignProblem(...)
            P.score('MMMSSISSD') #=> 23.50
        """
        subst_scores = self.params.subst_scores
        score = 0
        i, j = self.S_min_idx, self.T_min_idx
        for op,num in hp_tokenize(opseq):
            if op in 'MS':
                for k in range(num):
                    score += subst_scores[self.S.c_idxseq[i+k]][self.T.c_idxseq[j+k]]
                i, j = i + num, j + num
            elif op in 'ID':
                score += self.params.gap_open_score + self.params.gap_extend_score * num
                if op == 'I':
                    j = j + num
                else:
                    i = i + num
            else:
                raise ValueError('Invalid edit operation: %c' % op)
        return score

    def solve(self, print_dp_table=False):
        """Populates the DP table and traces back an (any) optimal alignment.

        :param print_dp_table(optional): whether or not (truthy) to print the
            fully calculated DP table.
        :returns: a transcript string with the specified format.
        """
        global lib
        #print 'S range: %d, %d' % (self.S_min_idx, self.S_max_idx)
        #print 'T range: %d, %d' % (self.T_min_idx, self.T_max_idx)
        self.c_dp_table = lib.define(self.c_obj)
        if self.c_dp_table == -1:
            raise('Got -1 from align.define().')

        self.opt = lib.solve(self.c_dp_table, self.c_obj)
        if print_dp_table:
            mat = self.dp_table
            for i in range(len(mat)):
                print [round(f,2) for f in mat[i]]
        if self.opt == ffi.NULL:
            print 'Err: No Alignment found (go=%.2f, ge=%.2f, max_div=%d)' % \
                (self.params.gap_open_score, self.params.gap_extend_score, self.params.max_diversion)
            self.opt = None
            return None, None
        score = self.opt.choices[0].score
        # TODO how do I do this?
        # lib.free(self.c_dp_table)
        return self.opt, score

    def traceback(self):
        if self.opt is None:
            return None

        rtranscript = lib.traceback(self.c_dp_table, self.c_obj, self.opt)
        if rtranscript == ffi.NULL:
            return None
        return Transcript(ffi.string(rtranscript))

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

    def __setattr__(self, name, value):
        """TODO
        """
        if name in ['S_min_idx', 'S_max_idx', 'T_min_idx', 'T_max_idx', 'type']:
            setattr(self.c_obj, name, value)
        else:
            return super(AlignProblem, self).__setattr__(name, value)


class Transcript(object):
    """A wrapper for alignment transcripts. Solutions to the alignment problem are
    represented by transcript strings with the following format:

        (<Si,Tj>),<score>:...

    Si and Tj are integers specifying the positiong along each string where
    the alignment begins. Score is the score of the transcript to 2 decimal
    places. What follows the ':' is a sequence of "ops" defined as follows:
        M match
        S substitution
        I insert
        D delete
    All op sequences begin with a B and insertion/deletions are meant to
    mean "from S to T".
    """
    def __init__(self, rtranscript=None):
        """Parses a raw transcript as given by `libalign.so` and populates class
        attributes. Raw transcripts are expected to have the following format:

            (<idx_S,idx_T>),<score>:<opseq>
        """
        if rtranscript is None:
            self.idx_S, self.idx_T, self.score, self.opseq = 0, 0, 0, ''
            return

        print rtranscript
        assert(re.match('\([0-9]+,[0-9]+\),[0-9-\.]+:[MISD]+', rtranscript) is not None)
        infostr, opseq = rtranscript.split(':', 1)
        indices, score = infostr.rsplit(',', 1)
        idx_S, idx_T = indices[1:-1].split(',') # skip the open/close parens
        self.idx_S, self.idx_T = int(idx_S), int(idx_T)
        self.opseq, self.score = opseq, float(score)

    def __repr__(self):
        return '(%d,%d),%.2f:%s' % (self.idx_S, self.idx_T, self.score, self.opseq)

    # TODO allow appending another transcript to the front or back

    def pretty_print(self, S, T, f, width=120, margin=20, colors=True):
        """Pretty prints a transcript to f.

        :param S(Sequence): "from" sequence.
        :param T(Sequence): "to" sequence.
        :param f: file handle to write the output to; for standard output use
            `sys.stdout`.
        :param width(optional): terminal width (int) used for wrapping; default
            is 120.
        :param margin(optional): length (int) of leading and trailing sequences
            in S and T before and after the alignment; default 20.
        :param colors(optional): whether or not (truthy) to use colors in
            output; default True.
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
        sline = 'S[%d]: ' % max(0, idx_S - pre_margin)
        tline = 'T[%d]: ' % max(0, idx_T - pre_margin)
        counter = max(len(sline), len(tline))
        for i in reversed(range(1, pre_margin)):
            if counter >= width:
                counter, sline, tline = new_line(sline, tline, idx_S+i, idx_T+i, f)
            sline += S[idx_S-i] if i <= idx_S else ' '
            tline += T[idx_T-i] if i <= idx_T else ' '
            counter += letlen

        # The alignment itself:
        for i,op in enumerate(self.opseq):
            if counter >= width:
                counter, sline, tline = new_line(sline, tline, idx_S, idx_T, f)
            if op in 'MS':
                s, t = S[idx_S], T[idx_T]
                idx_S += 1
                idx_T += 1
            elif op == 'I':
                s, t = '-', T[idx_T]
                idx_T += 1
            elif op == 'D':
                s, t = S[idx_S], '-'
                idx_S += 1
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
        post_margin = min(margin, max(S.length-idx_S, T.length-idx_T))
        for i in range(post_margin):
            if counter >= width:
                counter, sline, tline = new_line(
                    sline, tline, idx_S+i, idx_T+i, f)
            sline += S[idx_S+i] if idx_S + i < S.length else ' '
            tline += T[idx_T+i] if idx_T + i < T.length else ' '
            counter += letlen

        print_lines(sline, tline, f)
