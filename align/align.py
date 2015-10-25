from . import ffi, lib, utils, CffiObject
from .distillery import hp_tokenize

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
        gap_open_score=0, gap_extend_score=-1, max_diversion=10):
        self.alphabet = alphabet
        # each row in the subst matrix must be "owned" by an object that's kept
        # alive.
        self._c_subst_rows_ka = [ffi.new('double[]', subst_scores[i]) for i in range(self.alphabet.length)]
        self._c_subst_full_ka = ffi.new('double *[]', self._c_subst_rows_ka)
        self.c_obj = ffi.new('align_params*', {
            'alphabet': self.alphabet.c_obj,
            'subst_scores': self._c_subst_full_ka,
            'gap_open_score': gap_open_score,
            'gap_extend_score': gap_extend_score,
            'max_diversion': max_diversion,
        })

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
        S_max_idx = min(S.length, S_max_idx) if S_max_idx else S.length
        T_max_idx = min(T.length, T_max_idx) if T_max_idx else T.length

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
        self.c_dp_table = lib.define(self.c_obj)
        if self.c_dp_table == -1:
            raise('Got -1 from align.define().')

    def score(self, opseq):
        """Calculates the score for an arbitray opseq. Opseqs are allowed to be
        partial alignments (i.e finishing before reaching the end of frame).

            P = AlignProblem(...)
            P.score('BMMMSSISSD') #=> 23.50
        """
        assert(opseq[0] == 'B')
        opseq = opseq[1:]
        subst_scores = self.params.subst_scores
        score = 0
        i, j = self.S_min_idx, self.T_min_idx
        for op,num in hp_tokenize(opseq):
            if op == 'M':
                score += num * subst_scores[self.S.c_idxseq[i]][self.T.c_idxseq[j]]
                i, j = i+num, j+num
            elif op == 'S':
                for k in range(num):
                    score += subst_scores[self.S.c_idxseq[i+k]][self.T.c_idxseq[j+k]]
                i, j = i+num, j+num
            elif op in 'ID':
                score += self.gap_open_score + self.gap_extend_score * num
                if op == 'I':
                    j = j + num
                else:
                    j = i + num
            else:
                raise ValueError('Invalid edit operation: %c' % op)
        return score

    def solve(self, print_dp_table=False):
        """Populates the DP table and traces back an (any) optimal alignment.
        Solutions to the alignment problem are represented by transcript
        strings with the following format:

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

        :param print_dp_table(optional): whether or not (truthy) to print the
            fully calculated DP table.
        :returns: a transcript string with the specified format.
        """
        global lib
        if -1 == lib.solve(self.c_dp_table, self.c_obj):
            raise RuntimeError('Got -1 from align.solve()')
        if (print_dp_table):
            mat = self.dp_table
            for i in range(len(mat)):
                print mat[i]
        ret = lib.traceback(self.c_dp_table, self.c_obj)
        if ret == ffi.NULL:
            return 'Err: No Alignment found (go=%.2f, ge=%.2f, max_div=%d)' % \
                (self.params.gap_open_score, self.params.gap_extend_score, self.params.max_diversion)
        return ffi.string(ret)

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

