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
from . import sequence  # FIXME


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

# banded alignment types:
B_GLOBAL = lib.GLOBAL
"""Banded global alignment problem; may not be well-defined (end points of the
table may not lie in band)."""
B_OVERLAP = lib.B_OVERLAP
"""Banded suffix-prefix alignment problem in either direction including
substring alignments."""


class AlignScores(object):
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
            self._c_gaps_ka = ffi.new(
                'double []', content_dependent_gap_scores
            )
            kw['content_dependent_gap_scores'] = self._c_gaps_ka
        else:
            self._c_gaps_ka = ffi.NULL
        self.c_obj = ffi.new('alnscores*', kw)


class AlignFrame(object):
    def __init__(self, S, T, **kw):
        assert isinstance(S, sequence.Sequence)
        assert isinstance(T, sequence.Sequence)
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


class AlignTable(object):
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
                               OVERLAP, START_ANCHORED_OVERLAP,
                               END_ANCHORED_OVERLAP])
            self.c_alnparams = ffi.new('std_alnparams*', {'type': alntype})
            alnprob_args['std_params'] = self.c_alnparams
        elif alnmode == BANDED_MODE:
            assert(alntype in [B_GLOBAL, B_OVERLAP])
            self.c_alnparams = ffi.new('banded_alnparams*', {
                'type': alntype, 'dmin': kw['dmin'], 'dmax': kw['dmax']
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
            self.frame.S, self.frame.T, opseq,
            self.frame.S_range.i, self.frame.T_range.i
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
            Alignment
        """
        if self.opt is None:
            return None

        alignment = lib.dptable_traceback(self.c_obj, self.opt)
        if alignment == ffi.NULL:
            return None
        return Alignment(c_obj=alignment)


class Alignment(object):
    """Wraps alignment solutions represented as C ``alignment*``.
    All keyword arguments become attributes with identical names.

    Keyword Args:
        S_idx (int): The starting position in the "from" sequence.
        T_idx (int): The starting position in the "to" sequence.
        score (float): The score of the alignment.
        opseq (str): The sequence of edit "ops", see :class:`EditTranscript
            <biseqt.sequence.EditTranscript>` for format.
        c_obj (cffi.cdata): If provided all other arguments are
            ignored and instead this is used as the underlying C
            ``alignment *``.
    """
    # FIXME make opseq our EditTranscript?
    def __init__(self, **kw):
        if 'c_obj' in kw:
            self.c_obj = kw['c_obj']
            self.c_opseq = self.c_obj.opseq
        else:
            self.c_opseq = ffi.new('char[]', kw['opseq'])
            self.c_obj = ffi.new('alignment*', {
                'S_idx': kw['S_idx'],
                'T_idx': kw['T_idx'],
                'score': kw['score'],
                'opseq': self.c_opseq,
            })

    def __getattr__(self, name):
        if name == 'opseq':
            length = lib.strlen(self.c_opseq)
            return ''.join([self.c_opseq[i] for i in range(length)])
        else:
            return super(Alignment, self).__getattr__(name)

    # FIXME what is the use of this?
    def __str__(self):
        return '(%d,%d),%.2f:%s' \
            % (self.S_idx, self.T_idx, self.score, self.opseq)

    def __repr__(self):
        return 'Alignment(S_idx=%d, T_idx=%d, score=%.2f, opseq="%s")' \
            % (self.S_idx, self.T_idx, self.score, self.opseq)


# FIXME make this simpler and better
def rasterplot(path, alignment, S_name='S', T_name='T', fullview=False):
    if 'plt' not in globals():
        raise ImportError('matplotlib is required for rasterplots')
    ts, ss, cs = [alignment.T_idx], [alignment.S_idx], ['k']
    nums = {'M': 0, 'S': 0, '-': 0}
    colormap = {'M': 'g', 'S': 'y', '-': 'r'}
    for op in alignment.opseq:
        ts += [ts[-1] + 1 if op in 'MSI' else ts[-1]]
        ss += [ss[-1] + 1 if op in 'MSD' else ss[-1]]
        cs += colormap[op if op in 'MS' else '-']
        nums[op if op in 'MS' else '-'] += 1

    plt.rc('text', usetex=True)
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
    ax.set_xticks(list(ax.get_xticks())[1:-2] + [min(ts), max(ts)])
    ax.set_yticks(list(ax.get_yticks())[1:-1] + [min(ss), max(ss)])
    ax.set_ylabel(S_name)
    ax.xaxis.set_tick_params(labeltop='on')
    ax.xaxis.set_label_position('top')
    ax.invert_yaxis()
    ax.set_xlabel(T_name)
    # inset plot: op stats
    width = 0.07
    inset = fig.add_axes([0.9 - 4*width, .65, 4*width, 0.2], frameon=False)
    ind = [width * i for i in range(2)]
    ind = [width * i * 1.3 for i in range(3)]
    inset.bar(
        ind,
        [nums[i]*1./len(alignment.opseq) for i in 'MS-'],
        width,
        color=[colormap[i] for i in 'MS-']
    )
    inset.set_aspect('equal')
    inset.set_xticks([i+width/2. for i in ind])
    inset.set_xticklabels(['M', 'S', '-'])
    inset.set_yticks([i * .2 for i in range(1, 6)])
    fig.savefig(path)
