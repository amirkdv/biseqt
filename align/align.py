#!/usr/bin/env python
from . import ffi
from . import lib
from . import utils

global lib
"""Values for enum `align_type`.
"""
ALIGN_GLOBAL = lib.GLOBAL
ALIGN_LOCAL = lib.LOCAL
ALIGN_START_ANCHORED = lib.START_ANCHORED
ALIGN_END_ANCHORED = lib.END_ANCHORED
ALIGN_OVERLAP = lib.OVERLAP

class CffiObject(object):
    """Generic cffi wrapper for C structs, delegates all unknown attributes to
    the underlying C pointer. Subclasses must populate self.c_obj in their
    constructors with a pointer to their underlying C struct.

    Attributes:
        c_obj (cffi.cdata): points to the the underlying C pointer.
    """
    def __init__(self, c_type, **kw):
        self.c_obj = ffi.new('%s *' % c_type)

    def __getattr__(self, name):
        return getattr(self.c_obj, name)

class Sequence(CffiObject):
    """Wraps a C char[] and keeps its length. Placeholder for potential
    additions.

    Attributes:
        c_obj (cffi.cdata): points to the the underlying C char[].
    """
    def __init__(self, string):
        self.c_obj = ffi.new('char[]', string)
        self._len = len(string)

    def __repr__(self):
        return ''.join([self.c_obj[i] for i in range(self._len)])


class AlignParams(CffiObject):
    """Wraps the C struct align_params, see `libalign.h`

    Attributes:
        alphabet (str|list): a string or a list of strings with equal length
        _c_subst_rows_ka (cffi.cdata): has ownership of (keeps alive) C pointers
            to rows in the substitution matrix.
        _c_subst_full_ka (cffi.cdata): has ownership of (keeps alive) the C
            pointer to the full substitution matrix.
        _alph_letters_ka (cffi.cdata): has ownership of (keeps alive) the C
            pointers to strings, each a letter of the alphabet.
        c_obj (cffi.cdata): points to the underlying `align_params` struct.
        c_alph (cffi.cdata): points to the underlying `sequence_alphabet` struct.
    """
    def __init__(self, alphabet='ACGT', subst_scores=[],
        gap_open_score=0, gap_extend_score=-1, max_diversion=10):
        if isinstance(alphabet, str):
            alphabet = [c for c in alphabet]
            alph_letter_len = 1
        assert(len(set([len(s) for s in alphabet])) == 1)
        self.alph_len, self.alph_let_len = len(alphabet), len(alphabet[0])
        # each letter string in the alphabet must be "owned" by an object
        # that's kept alive.
        self._c_alph_letters_ka = [ffi.new('char[]', alphabet[i]) for i in range(self.alph_len)]
        self._c_alph_ka = ffi.new('char *[]', self._c_alph_letters_ka)
        self.c_alph = ffi.new('sequence_alphabet*', {
            'length': self.alph_len,
            'letter_length': self.alph_let_len,
            'letters': self._c_alph_ka
        })
        # each row in the subst matrix must be "owned" by an object that's kept
        # alive.
        self._c_subst_rows_ka = [ffi.new('double[]', subst_scores[i]) for i in range(self.alph_len)]
        self._c_subst_full_ka = ffi.new('double *[]', self._c_subst_rows_ka)
        self.c_obj = ffi.new('align_params*', {
            'alphabet': self.c_alph,
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
            idx = range(self.alph_len)
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

        # the maximum *possible* last indices for S and T
        S_Max = S._len/params.alphabet.letter_length
        T_Max = T._len/params.alphabet.letter_length
        S_max_idx = min(S_Max, S_max_idx) if S_max_idx else S_Max
        T_max_idx = min(T_Max, T_max_idx) if T_max_idx else T_Max

        self.S, self.T, self.params, self.align_type = S, T, params, align_type
        self.c_obj = ffi.new('align_problem*', {
            'S': S.c_obj,
            'T': T.c_obj,
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

