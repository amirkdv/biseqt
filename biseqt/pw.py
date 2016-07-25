# -*- coding: utf-8 -*-
"""This modules provides various pairwise sequence alignment algorithms
implemented in `pwlib <biseqt.pwlib.html>`_."""
import os
import termcolor
import re
from cffi import FFI
from collections import namedtuple

from .sequence import Sequence


lib = None
"""The loaded shared object for ``pwlib``. All functions defined in the header
file are accessible through this object. This object is automatically populated
upon loading this module (cf. :func:`setup_ffi`) and users never have to
manipulate it."""

ffi = None
"""The main FFI instance used throughout this module. This object is
automatically populated upon loading this module (cf. :func:`setup_ffi`) and
users never have to manipulate it."""


def setup_ffi():
    """Instantiates an FFI object as :attr:`ffi` and loads the shared object
    for pwlib into :attr:`lib`. This function is automatically called when
    this module loads.

    Note:
        CFFI has issues with loading macros as they are defined in a header
        file. For this reason, and since we don't use the macros in python code
        any line that begins with ``#define`` is ignored from the header file.
        This means multiline macros will not work.
    """
    pwlib_so = os.path.join(os.path.dirname(__file__), 'pwlib', 'pwlib.so')
    pwlib_h = os.path.join(os.path.dirname(__file__), 'pwlib', 'pwlib.h')
    global ffi, lib
    ffi = FFI()
    lib = ffi.dlopen(pwlib_so)
    with open(pwlib_h) as f:
        # ignore macro definitions in header file; we don't use them in python
        # and cffi is unhappy about parsing them:
        headers = '\n'.join(line for line in f.read().split('\n')
                            if not line.startswith('#define'))
        ffi.cdef(headers)

setup_ffi()

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

ALN_TYPES = {
    STD_MODE: [GLOBAL, LOCAL, START_ANCHORED, END_ANCHORED, OVERLAP,
               START_ANCHORED_OVERLAP, END_ANCHORED_OVERLAP],
    BANDED_MODE: [B_GLOBAL, B_OVERLAP],
}


class Aligner(object):
    """Provides a context that solves a pairwise alignment problem.  Memory is
    allocated upon entering the context and is freed upon leaving it.  All
    alignment calculations (:func:`solve` and :func:`traceback`) are explicitly
    invoked by the caller:

        >>> A = biseqt.sequence.Alphabet('ACGT')
        >>> S = A.parse('AAAA')
        >>> T = A.parse('AGA')
        >>> aligner = biseqt.pw.Aligner(S, T)
        >>> with aligner:
        ...    print 'score is', aligner.solve()
        ...    print aligner.traceback()
        score is 2.0
        origin[0]: AAAA
        mutant[0]: AGA-

    Args:
        origin (sequence.Sequence): The original ("from") sequence.
        mutant (sequence.Sequence): The mutant ("to") sequence.

    Keyword Args:
        origin_range (tuple): The original ("from") sequence; cf.
            :c:member:`alnframe::origin_range`.
        mutant_range (tuple): The mutant ("to") sequence; cf.
            :c:member:`alnframe::mutant_range`.
        alnmode (int): One of the :attr:`STD_MODE` or :attr:`BANDED_MODE`,
            default is ``STD_MODE``; cf. :c:member:`alnprob::mode`.
        alntype (int): One of the allowed alingment types for the given
            *alnmode*, see :attr:`ALN_TYPES`; default is ``GLOBAL``; cf.
            :c:type:`std_alnparams` and :c:type:`banded_alnparams`.
        subst_scores (list): The overriding definition of the substitution
            score matrix; cf. :c:member:`alnscores::subst_scores`. Default is
            None in which case the score matrix is populated based on match and
            mismatch scores.
        match_score (float): If ``subst_scores`` is not given, this parameter
            is used to populate the diagonal entries of the substitution score
            matrix; default is 1.
        mismatch_score (float): If ``subst_scores`` is not given, this
            parameter is used to populate the off-diagonal entries of the
            substitution score matrix; default is 0.
        go_score (float): The gap open score; cf.
            :c:member:`alnscores::gap_open_score`. Default is 0.
        ge_score (float): The gap extend score; cf.
            :c:member:`alnscores::gap_extend_score`. Default is 0.
        max_new_mins (int): Maximum number of tolerated new minima encountered
            in the running score of an alignment; cf.
            :c:member:`alnprob::max_new_mins`. Default is -1 in which case
            no such constraint is imposed.
        diag_range (tuple): If in :attr:`BANDED_MODE` this argument specifies
            the upper and lower limit on diagonals of the dynamic programming
            table to be populated; cf. :c:type:`banded_alnparams`.
    """
    def __init__(self, origin, mutant, **kw):
        self.alnmode = kw.get('alnmode', STD_MODE)
        self.alntype = kw.get('alntype', GLOBAL)
        assert self.alnmode in [STD_MODE, BANDED_MODE]
        assert self.alntype in ALN_TYPES[self.alnmode]

        # set origin, mutant, and alphabet
        assert isinstance(origin, Sequence) and isinstance(mutant, Sequence)
        assert origin.alphabet == mutant.alphabet
        self.origin, self.mutant = origin, mutant
        self.alphabet = origin.alphabet

        # set origin_range and mutant_range
        origin_range = kw.get('origin_range', (0, len(self.origin)))
        mutant_range = kw.get('mutant_range', (0, len(self.mutant)))
        assert 0 <= origin_range[0] <= origin_range[1] <= len(self.origin)
        assert 0 <= mutant_range[0] <= mutant_range[1] <= len(self.mutant)
        self.origin_range, self.mutant_range = origin_range, mutant_range

        # set alignment scores
        self.go_score = kw.get('go_score', 0)
        self.ge_score = kw.get('ge_score', 0)
        L = len(self.alphabet)
        subst_scores = kw.get('subst_scores', None)
        if subst_scores is None:
            mismatch = kw.get('mismatch_score', 0)
            match = kw.get('match_score', 1)
            subst_scores = [
                [match if i == j else mismatch for i in range(L)]
                for j in range(L)
            ]
        assert isinstance(subst_scores, list) and len(subst_scores) == L
        self.subst_scores = subst_scores

        self.max_new_mins = kw.get('max_new_mins', -1)
        self.diag_range = kw.get('diag_range', None)

        # create all the C data structures
        self.c_subst_scores_rows = [ffi.new('double[]', self.subst_scores[i])
                                    for i in range(L)]
        self.c_subst_scores = ffi.new('double *[]', self.c_subst_scores_rows)
        self.c_alnscores = ffi.new('alnscores*', {
            'subst_scores': self.c_subst_scores,
            'gap_open_score': self.go_score,
            'gap_extend_score': self.ge_score,
        })
        self.c_origin = ffi.new('int[]', self.origin.contents)
        self.c_mutant = ffi.new('int[]', self.mutant.contents)
        self.c_alnframe = ffi.new('alnframe*', {
            'origin': self.c_origin,
            'mutant': self.c_mutant,
            'origin_range': self.origin_range,
            'mutant_range': self.mutant_range,
        })

        if self.alnmode == STD_MODE:
            self.c_alnparams = ffi.new('std_alnparams*',
                                       {'type': self.alntype})
        elif self.alnmode == BANDED_MODE:
            self.min_diag, self.max_diag = kw['diag_range']
            assert -len(mutant) <= self.min_diag <= \
                self.max_diag <= len(origin)
            self.c_alnparams = ffi.new('banded_alnparams*', {
                'type': self.alntype,
                'dmin': self.min_diag,
                'dmax': self.max_diag,
            })
        self.c_alnprob = ffi.new('alnprob*', {
            'frame': self.c_alnframe,
            'scores': self.c_alnscores,
            'mode': self.alnmode,
            'max_new_mins': self.max_new_mins,
            'std_params' if self.alnmode == STD_MODE else 'banded_params':
                self.c_alnparams,
        })
        self.c_dptable = ffi.new('dptable*', {
            'prob': self.c_alnprob,
            'cells': ffi.NULL,
            'num_rows': -1,
            'row_lens': ffi.NULL,
        })

    def __enter__(self):
        """Allocates memory for the dynamic programming table and initializes
        all cells."""
        if lib.dptable_init(self.c_dptable) == -1:
            raise Exception('Failed to initialize the DP table.')
        return self

    def __exit__(self, *args):
        """Frees the allocated memory for the dynamic programming table."""
        lib.dptable_free(self.c_dptable)

    def solve(self):
        """Populates the regions of interest in the dynamic programming table
        and reports the optimal score; if any. This function must be called
        within the context, cf. :func:`__enter__`, :func:`__exit__`.

        Returns:
            score (float): The score of the optimal alignment or None if none
                found.
        """
        self.opt = lib.dptable_solve(self.c_dptable)
        if self.opt.i == -1 or self.opt.j == -1:
            self.opt = None
            return None
        return self.c_dptable.cells[self.opt.i][self.opt.j].choices[0].score

    def traceback(self):
        """Traces back the optimal alignment identified by :func:`solve`. This
        function has to be called within the context and after :func:`solve`.
        Otherwise no alignment would be found.

        Returns:
            Alignment: The optimal alignment or None if none found.
        """
        if self.opt is None:
            return None

        alignment = lib.dptable_traceback(self.c_dptable, self.opt)
        assert alignment != ffi.NULL
        transcript = ffi.string(alignment.transcript)
        if not transcript:
            return None
        return Alignment(self.origin, self.mutant, transcript,
                         score=alignment.score,
                         origin_start=alignment.origin_idx,
                         mutant_start=alignment.mutant_idx)

    def calculate_score(self, alignment):
        """Scores a given alignment for :attr:`origin` and :attr:`mutant`.
        Args:
            alignment (Alignment): The alignment to be evaluated.

        Returns:
            float:
                The score of the alignment based on :attr:`subst_scores`,
                :attr:`go_score`, and :attr:`ge_score`.
        """
        return alignment.calculate_score(self.subst_scores, self.go_score,
                                         self.ge_score)


class Alignment(object):
    """Represents a pairwise alignment.

    Attributes:
        origin (sequence.Sequence): The original ("from") sequence.
        mutant (sequence.Sequence): The mutant ("to") sequence.
        alphabet (sequence.Alphabet): The shared alphabet of *origin* and
            *mutant*.
        transcript (str): The sequence of edit operations that transforms
            *origin* to *mutant*. The alphabet for edit operations is ``M`` for
            match, ``S`` for substitution (mismatch), and ``I`` and ``D`` for
            insertion and deletion.
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

    def __eq__(self, other):
        assert isinstance(other, Alignment)
        return other.origin == self.origin and \
            other.mutant == self.mutant and \
            other.transcript == self.transcript and \
            other.origin_start == self.origin_start and \
            other.mutant_start == self.mutant_start

    @classmethod
    def projected_len(cls, transcript, on='origin'):
        """Calculates the projected length of a given transcript on either of
        the involved sequences. For instance:

            >>> biseqt.Alignment.projected_len('MSI', on='origin')
            2
            >>> biseqt.Alignment.projected_len('MSI', on='mutant')
            3

        Args:
            transcript (str): A sequence of edit operations, cf.
                :attr:`Alignment.transcript`.

        Keyword Args:
            on (str): Either of ``origin`` or ``mutant``.

        Returns:
            int: The projected length of the edit transcript.
        """
        assert on in ['origin', 'mutant']
        ops = 'MSD' if on == 'origin' else 'MSI'
        return sum(int(op in ops) for op in transcript)

    def calculate_score(self, subst_scores, go_score, ge_score):
        """Scores a this alignment according to given scoring scheme.

        Args:
            subst_scores (list): The substitution score matrix, cf.
                :attr:`Aligner.subst_scores`.
            go_score (float): The gap open score; cf. :attr:`Aligner.go_score`.
            ge_score (float): The gap extend score; cf.
                :attr:`Aligner.ge_score`.

        Returns:
            float:
                The score of the alignment for :attr:`origin` and
                :attr:`mutant` based on given scores.
        """
        score = 0.
        i, j = self.origin_start, self.mutant_start

        def tokens():
            for match in re.finditer(r'(.)\1*', self.transcript):
                match = match.group(0)
                yield match[0], len(match)

        for op, num in tokens():
            if op in 'MS':
                score += sum(
                    subst_scores[self.origin[i + k]][self.mutant[j + k]]
                    for k in range(num)
                )
                i, j = i + num, j + num
            else:
                assert op in 'ID'
                score += go_score + ge_score * num
                if op == 'I':
                    j = j + num
                else:
                    i = i + num
        return score

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
