# -*- coding: utf-8 -*-
import termcolor
from collections import namedtuple
# FIXME add code snippets
# FIXME figure out how to run code snippets
"""
.. code-block:: python

    >>> from biseqt.sequence import Sequence, Alphabet
"""


class Alphabet(object):
    """A sequence alphabet.

    Attributes:
        _letters (tuple):
            The letters in the alphabet. All ``getitem`` operations (i.e
            indexing and slicing) are delegated to this tuple. This attribute
            should be considered read-only.

        _letlen (int):
            The length of the letters in the alphabet when represented as a
            string. This attribute should be considered read-only.
    """
    def __init__(self, letters):
        """
        Args:
            letters (iterable):
                The elements of this iterable must be hashable, i.e can be
                keys of a dictionary, and must respond to :func:`len`.
                Typically, they are single character strings.
        """
        self._letters = tuple(letters)
        self._letlen = len(self._letters[0])
        assert all(len(l) == self._letlen for l in self._letters), \
            'All alphabet letters must have the same length'
        self._idx_by_letter = {l: idx for idx, l in enumerate(self._letters)}

    def to_idx(self, letters):
        """Translates provided letters to the integer sequence corresponding
        to the index of each letter in this alphabet.

        Args:
            letters (iterable):
                The original sequence whose elements are letters of this
                alphabet.

        Returns:
            tuple
        """
        return tuple(self._idx_by_letter[l] for l in letters)

    def parse(self, string):
        """Given a string representation of a sequence returns a corresponding
        :class:`Sequence` object.
        """
        assert len(string) % self._letlen == 0, 'String representation ' + \
            'of sequence must be a multiple of the alphabet letter length'
        contents = []
        idx = 0
        while idx < len(string):
            contents.append(string[idx:idx + self._letlen])
            idx += self._letlen
        return Sequence(self, self.to_idx(contents))

    def __len__(self):
        return len(self._letters)

    def __eq__(self, other):
        assert isinstance(other, Alphabet), \
            'Only alphabets can be compared with alphabets'
        return self._letters == other._letters

    def __getitem__(self, key):
        return self._letters.__getitem__(key)

    def __repr__(self):
        return 'Alphabet([%s])' % \
            ','.join('"%s"' % self[idx] for idx in range(len(self)))


class Sequence(object):
    """An immutable sequence of letters from some :class:`Alphabet` which
    behaves mostly like a tuple.

    Attributes:
        alphabet (Alphabet): The :class:`Alphabet` of the sequence.
        contents (tuple): The contents of the sequence represented as tuple of
            integers of the same length where each letter is represented by
            its position in the alphabet.
    """
    def __init__(self, alphabet, contents=()):
        """Initializes the sequence object: translates all letters to integers
        corresponding to the position of each letter in the alphabet.

        Args:
            alphabet (Alphabet):
                The :class:`Alphabet` of the sequence.
            contents (iterable):
                The contents of the sequence as an iterable, each element of
                which is the integer representation of a letter from the
                :class:`Alphabet`; default is an empty sequence. If the
                alphabet letter length is one, this argument can be a string.
        """
        assert isinstance(alphabet, Alphabet)
        self.alphabet = alphabet

        assert all(isinstance(c, int) and c < len(alphabet) for c in contents)
        self.contents = tuple(contents)

    def __str__(self):
        return ''.join(self.alphabet[idx] for idx in self.contents)

    def __repr__(self):
        return 'Sequence(%s, %s)' % (repr(self.alphabet), repr(self.contents))
        return ''.join([self.__getitem__(i) for i in range(len(self))])

    def __len__(self):
        return len(self.contents)

    def __nonzero__(self):
        return True if self.contents else False

    def __getitem__(self, key):
        if isinstance(key, int):
            return self.contents[key]
        else:
            return Sequence(self.alphabet, self.contents.__getitem__(key))

    def __eq__(self, other):
        return self.alphabet == other.alphabet and \
               self.contents == other.contents

    def __add__(self, other):
        if isinstance(other, Sequence):
            assert self.alphabet == other.alphabet
            contents = other.contents
        else:
            other = self.alphabet.to_idx(other)
        return Sequence(self.alphabet, self.contents + contents)


class EditTranscript(str):
    """Represents a edit transcript without reference to the original or mutant
    sequence."""
    def __new__(cls, content):
        obj = str.__new__(cls, content.upper())
        assert all(c in 'MSID' for c in obj)
        return obj

    def __repr__(self):
        return 'EditTranscript("%s")' % self

    def __getitem__(self, key):
        if isinstance(key, int):
            return str(self)[key]
        else:
            return EditTranscript(str(self)[key])

    def __add__(self, other):
        return EditTranscript(str(self) + str(other))

    def render_term(self, origin, mutant, origin_start=0, mutant_start=0,
                    term_width=120, margin=0, colored=True):
        """Pretty prints a transcript to f.

        Args:
            origin (Sequence): The original :class:`Sequence`.
            mutant (Sequence): The mutant :class:`Sequence`.

        Keyword Args:
            origin_start (int): Starting position on the original sequence;
                default is 0.
            mutant_start (int): Starting position on the mutant sequence;
                default is 0.
            term_width (int): Terminal width used for wrapping; default is 120
                and the smallest valid value is 30.
            margin (length): Length of leading and trailing substring to
                include in original and mutant sequences; default is 20.
            colored (bool): Whether or not to use colors in output; default is
                True.
        """
        assert origin_start >= 0
        assert mutant_start >= 0
        assert origin_start < len(origin)
        assert mutant_start < len(mutant)
        assert term_width >= 30
        assert margin >= 0
        assert origin.alphabet == mutant.alphabet
        alphabet = origin.alphabet
        letlen = alphabet._letlen

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
                if carriage.o_idx >= 0 and carriage.o_idx < len(origin):
                    o_contents = alphabet[origin[carriage.o_idx]]
                if carriage.m_idx >= 0 and carriage.m_idx < len(mutant):
                    m_contents = alphabet[mutant[carriage.m_idx]]
            else:
                assert op in 'MSID'
                if op in 'MSD':
                    o_contents = alphabet[origin[carriage.o_idx]]
                if op in 'MSI':
                    m_contents = alphabet[mutant[carriage.m_idx]]

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
                             max((len(origin) - carriage.o_idx) * letlen,
                                 (len(mutant) - carriage.m_idx) * letlen))
            for i in range(margin_len):
                out, carriage = carriage_fwd(carriage, op=None)
                output += out
            return output + carriage_flush(carriage)

        output, carriage = pre_margin(origin_start, mutant_start)
        for op in self:
            out, carriage = carriage_fwd(carriage, op=op)
            output += out
        output += post_margin(carriage)

        # when output is not supposed to be cleared remove the spurious color
        # reset ANSI escape sequence that termcolor adds:
        if not colored:
            output = output.replace(termcolor.RESET, '')

        return output
