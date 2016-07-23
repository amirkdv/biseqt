# -*- coding: utf-8 -*-
"""This modules provides various pairwise sequence alignment algorithms
implemented in `pwlib <biseqt.pwlib.html>`_."""
import os
import termcolor
from cffi import FFI
from collections import namedtuple
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

        self.origin = origin
        self.mutant = mutant

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
        return Alignment(self.frame.origin, self.frame.mutant,
                         ffi.string(alignment.transcript),
                         score=alignment.score,
                         origin_start=alignment.origin_idx,
                         mutant_start=alignment.mutant_idx)


class Alignment(object):
    """Represents a pairwise alignment.

    Attributes:
        origin (sequence.Sequence): The original :class:`Sequence`.
        mutant (sequence.Sequence): The mutant :class:`Sequence`.
        alphabet (sequence.Alphabet): The shared alphabet of *origin* and
            *mutant*.
        transcript (str): The sequence of edit operations that mutates *origin*
            to *mutant*.
        origin_start (int): Starting position on the original sequence; default
            is 0.
        mutant_start (int): Starting position on the mutant sequence; default
            is 0.
        score (float): The score of the alignment; default is None.
    """
    def __init__(self, origin, mutant, transcript, score=None,
                 origin_start=0, mutant_start=0):
        assert isinstance(origin, Sequence) and isinstance(mutant, Sequence)
        assert origin.alphabet == mutant.alphabet
        self.alphabet = origin.alphabet
        assert all(c in 'MSID' for c in transcript)
        assert len(transcript) > 0
        origin_end = origin_start + self.projected_len(transcript, on='origin')
        mutant_end = mutant_start + self.projected_len(transcript, on='mutant')
        assert 0 <= origin_start and origin_end <= len(origin)
        assert 0 <= mutant_start and mutant_end <= len(mutant)
        self.transcript = str(transcript)
        self.origin, self.mutant = origin, mutant
        self.origin_start, self.mutant_start = origin_start, mutant_start
        self.score = score

    def __str__(self):
        return self.render_term(term_width=float('+inf'), margin=0, colored=0)

    @classmethod
    def projected_len(cls, transcript, on='origin'):
        assert on in ['origin', 'mutant']
        ops = 'MSD' if on == 'origin' else 'MSI'
        return sum(int(op in ops) for op in transcript)

    def render_term(self, term_width=120, margin=0, colored=True):
        """Renders a textual representation of the alignment.

        Keyword Args:
            term_width (int): Terminal width used for wrapping; default is 120
                and the smallest valid value is 30.
            margin (length): Length of leading and trailing substring to
                include in original and mutant sequences; default is 20.
            colored (bool): Whether or not to use ANSI color codes in output;
                default is True.

        Returns:
            str
        """
        assert term_width >= 30
        assert margin >= 0
        letlen = self.alphabet._letlen

        Carriage = namedtuple('carriage',
                              ['pos', 'o_idx', 'm_idx', 'o_line', 'm_line'])

        term_color = {'M': 'green', 'S': 'red'}
        term_on_color = {'I': 'on_red', 'D': 'on_red'}

        # In the rest: o_X and m_X mean X for origin and mutatnt, resp.

        # Creates an alignment line preamble, i.e a double line for origin and
        # sequence, given starting positions on each. The output is a tuple
        # (pos, o_line, m_line) where pos is the position in line after the
        # preamble.
        def start_line(o_idx, m_idx):
            kw = {
                'o_idx': o_idx,
                'm_idx': m_idx,
                'o_line': 'origin[%d]: ' % o_idx,
                'm_line': 'mutant[%d]: ' % m_idx,
            }
            pos = max(len(kw['o_line']), len(kw['m_line']))
            assert pos <= term_width, \
                'Alignment preamble does not fit in width %d' % term_width
            return Carriage(pos=pos, **kw)

        # returns a right adjusted double line given the two lines of an
        # alignment, i.e the origin and mutant versions.
        def carriage_flush(carriage):
            line_len = max(len(carriage.o_line), len(carriage.m_line))
            o_line = carriage.o_line.rjust(line_len)
            m_line = carriage.m_line.rjust(line_len)
            return '%s\n%s\n' % (o_line, m_line)

        def carriage_fwd(carriage, op=None):
            gap = '.' * letlen if op is None else '-' * letlen
            o_contents, m_contents = gap, gap
            if op is None:
                if carriage.o_idx >= 0 and carriage.o_idx < len(self.origin):
                    o_contents = self.alphabet[self.origin[carriage.o_idx]]
                if carriage.m_idx >= 0 and carriage.m_idx < len(self.mutant):
                    m_contents = self.alphabet[self.mutant[carriage.m_idx]]
            else:
                assert op in 'MSID'
                if op in 'MSD':
                    o_contents = self.alphabet[self.origin[carriage.o_idx]]
                if op in 'MSI':
                    m_contents = self.alphabet[self.mutant[carriage.m_idx]]

            length = len(o_contents)
            assert length == len(m_contents)

            colors = {'color': None, 'on_color': None}
            if colored and op in term_color:
                colors['color'] = term_color[op]
            if colored and op in term_on_color:
                colors['on_color'] = term_on_color[op]
            o_contents = termcolor.colored(o_contents, **colors)
            m_contents = termcolor.colored(m_contents, **colors)

            output = ''
            if carriage.pos >= term_width:
                output += carriage_flush(carriage)
                carriage = start_line(carriage.o_idx, carriage.m_idx)

            return output, Carriage(
                pos=carriage.pos + length,
                o_idx=carriage.o_idx + int(op is None or op in 'MSD'),
                m_idx=carriage.m_idx + int(op is None or op in 'MSI'),
                o_line=carriage.o_line + o_contents,
                m_line=carriage.m_line + m_contents
            )

        # the arguments are the starting positions in the origin/mutant.
        def pre_margin(o_idx, m_idx):
            margin_len = min(margin, max(o_idx, m_idx) * letlen)
            carriage = start_line(o_idx - margin_len, m_idx - margin_len)
            output = ''
            # the pre-margin
            for i in range(margin_len):
                out, carriage = carriage_fwd(carriage, op=None)
                output += out
            return output, carriage

        # the arguments are the ending positions in the origin/mutant.
        def post_margin(carriage):
            output = ''
            margin_len = min(margin,
                             max((len(self.origin) - carriage.o_idx) * letlen,
                                 (len(self.mutant) - carriage.m_idx) * letlen))
            for i in range(margin_len):
                out, carriage = carriage_fwd(carriage, op=None)
                output += out
            return output + carriage_flush(carriage)

        output, carriage = pre_margin(self.origin_start, self.mutant_start)
        for op in self.transcript:
            out, carriage = carriage_fwd(carriage, op=op)
            output += out
        output += post_margin(carriage)

        # when output is not supposed to be cleared remove the spurious color
        # reset ANSI escape sequence that termcolor adds:
        if not colored:
            output = output.replace(termcolor.RESET, '')

        return output
