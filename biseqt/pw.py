# -*- coding: utf-8 -*-
"""This modules provides various pairwise sequence alignment algorithms
implemented in `pwlib <biseqt.pwlib.html>`_."""
import os
from cffi import FFI
from math import log
import re
try:
    from matplotlib import pyplot as plt
except ImportError:
    pass
from .sequence import Sequence, Alphabet


pwlib_so = os.path.join(os.path.dirname(__file__), 'pwlib', 'pwlib.so')
pwlib_h = os.path.join(os.path.dirname(__file__), 'pwlib', 'pwlib.h')

ffi = FFI()
lib = ffi.dlopen(pwlib_so)
with open(pwlib_h) as f:
    ffi.cdef(f.read())


# alignment modes
STD_MODE = lib.STD_MODE
"""Standard alignment type; time and memory complexity is quadratic in sequence
lengths."""
BANDED_MODE = lib.BANDED_MODE
"""Banded alignment type; time and memory complexity is linear in sequence
lengths with a constant proportional to band width. This mode is incompatible
with local alignments."""

# standard alignment types:
GLOBAL = lib.GLOBAL
"""Standard global alignment problem, i.e Needleman-Wunsch."""
LOCAL = lib.LOCAL
"""Standard local alignment problem, i.e Smith-Waterman."""
START_ANCHORED = lib.START_ANCHORED
"""Standard local alignment demanding that it begins at the start of frame of
both sequences."""
END_ANCHORED = lib.END_ANCHORED
"""Standard local alignment demanding that it ends at the end of frame of both
sequences."""
OVERLAP = lib.OVERLAP
"""Standard suffix-prefix alignment in any direction; this includes alignments
where a prefix of either sequence matches a suffix of the other and alignments
where one sequence is a substring of the other."""
START_ANCHORED_OVERLAP = lib.START_ANCHORED_OVERLAP
"""Standard suffix-prefix alignment demanding that it begins at the start of
frame of both sequences."""
END_ANCHORED_OVERLAP = lib.END_ANCHORED_OVERLAP
"""Standard suffix-prefix alignment demanding that it ends at the end of frame
of both sequences."""

STD_ALN_TYPES = [GLOBAL, LOCAL, START_ANCHORED, END_ANCHORED, OVERLAP,
                 START_ANCHORED_OVERLAP, END_ANCHORED_OVERLAP]

# banded alignment types:
B_GLOBAL = lib.GLOBAL
"""Banded global alignment problem; may not be well-defined (end points of the
table may not lie in band)."""
B_OVERLAP = lib.B_OVERLAP
"""Banded suffix-prefix alignment problem in either direction including
substring alignments."""

BANDED_ALN_TYPES = [B_GLOBAL, B_OVERLAP]


# FIXME not needed at all! Instead have a python class defined in stochastics
# that wraps gap and substitution scores. Convert to c object when solving.
class AlignScoresC(object):
    def __init__(self, alphabet=None, subst_scores=[],
                 go_score=0, ge_score=0, content_dependent_gap_scores=None):
        assert isinstance(alphabet, Alphabet)
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
            self._c_gaps_ka = ffi.new(
                'double []', content_dependent_gap_scores
            )
            kw['content_dependent_gap_scores'] = self._c_gaps_ka
        else:
            self._c_gaps_ka = ffi.NULL
        self.c_obj = ffi.new('alnscores*', kw)


# FIXME not needed at all! move all this C-conversion crap to AlignTableC or
# the currently-nonexistent "align()"
class AlignFrameC(object):
    def __init__(self, origin, mutant, origin_range=None, mutant_range=None):
        assert isinstance(origin, Sequence) and isinstance(mutant, Sequence)
        if origin_range is None:
            origin_range = (0, len(origin))
        if mutant_range is None:
            mutant_range = (0, len(mutant))
        assert 0 <= origin_range[0] <= origin_range[1] <= len(origin)
        assert 0 <= mutant_range[0] <= mutant_range[1] <= len(mutant)

        self.c_origin = ffi.new('int[]', origin.contents)
        self.c_mutant = ffi.new('int[]', mutant.contents)
        self.c_obj = ffi.new('alnframe*', {
            'origin': self.c_origin,
            'mutant': self.c_mutant,
            'origin_range': origin_range,
            'mutant_range': mutant_range,
        })


class AlignTableC(object):
    def __init__(self, frame, scores, alnmode=STD_MODE, alntype=GLOBAL,
                 max_new_mins=-1, diag_range=None):
        assert isinstance(frame, AlignFrameC)
        assert isinstance(scores, AlignScoresC)
        self.frame, self.scores = frame, scores
        assert alnmode in [STD_MODE, BANDED_MODE]
        alnprob_args = {
            'frame': frame.c_obj,
            'scores': scores.c_obj,
            'mode': alnmode,
            'max_new_mins': max_new_mins,
        }
        if alnmode == STD_MODE:
            assert alntype in STD_ALN_TYPES
            self.c_alnparams = ffi.new('std_alnparams*', {'type': alntype})
            alnprob_args['std_params'] = self.c_alnparams
        elif alnmode == BANDED_MODE:
            assert alntype in BANDED_ALN_TYPES
            assert diag_range is not None and len(diag_range) == 2
            self.c_alnparams = ffi.new('banded_alnparams*', {
                'type': alntype,
                'dmin': diag_range[0],
                'dmax': diag_range[1],
            })
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

    def score(self, opseq):
        # FIXME
        return self.scores.score(
            self.frame.origin, self.frame.mutant, opseq,
            self.frame.origin_range.i, self.frame.mutant_range.i
        )

    def solve(self):
        self.opt = lib.dptable_solve(self.c_obj)
        if self.opt.i == -1 or self.opt.j == -1:
            self.opt = None
            return None
        score = self.c_obj.cells[self.opt.i][self.opt.j].choices[0].score
        return score

    def traceback(self):
        if self.opt is None:
            return None

        alignment = lib.dptable_traceback(self.c_obj, self.opt)
        if alignment == ffi.NULL:
            return None
        return AlignmentC(c_obj=alignment)


# FIXME not needed at all! Instead have a python class with the score()
# function.
class AlignmentC(object):
    def __init__(self, c_obj=None, **kw):
        if c_obj is not None:
            self.c_obj = c_obj
            self.c_tx = c_obj.transcript
        else:
            self.c_tx = ffi.new('char[]', kw['transcript'])
            self.c_obj = ffi.new('alignment*', {
                'origin_idx': kw['origin_idx'],
                'mutant_idx': kw['mutant_idx'],
                'score': kw['score'],
                'transcript': self.c_tx,
            })

    def __getattr__(self, name):
        if name == 'transcript':
            return ffi.string(self.c_tx)
        elif name in ['origin_idx', 'mutant_idx', 'score']:
            return getattr(self.c_obj, name)
        else:
            raise AttributeError('Unknown attribute %s' % name)

    def __str__(self):
        return '(%d,%d),%.2f:%s' \
            % (self.origin_idx, self.mutant_idx, self.score, self.transcript)

    def __repr__(self):
        return 'Alignment(origin_idx=%d, mutant_idx=%d, score=%.2f, opseq="%s")' \
            % (self.origin_idx, self.mutant_idx, self.score, self.transcript)
