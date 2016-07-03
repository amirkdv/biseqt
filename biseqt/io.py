# -*- coding: utf-8 -*-
import termcolor
import textwrap
from collections import namedtuple

from .sequence import NamedSequence


def read_fasta(f, alphabet, num=-1):
    """Lazy generator for content-identifiable :class:`NamedSequence
    <biseqt.sequence.NamedSequence>` objects loaded from a FASTA source.

    Args:
        f (file): A readable open file or equivalent (e.g.
            :class:`StringIO.StringIO`). The file position is only modified by
            reading; namely, sequences are read starting from the current
            position of the file and after this function returns the file
            position points at the last character of the last yielded sequence.
        alphabet (sequence.Alphabet): The alphabet of the sequences.
        num (int): The number of sequences to read from the file; default is -1
            in which case all sequences are read until end-of-file is reached.

    Yields:
        tuple:
            A :class:`NamedSequence <biseqt.sequence.NamedSequence>` whose
            :attr:`name <biseqt.sequence.NamedSequence.name>` is the FASTA
            record name, and the starting position of the sequence in ``f`` as
            an integer.
    """
    observed_names = []

    def _parse(seq, name):
        if cur_seq and cur_name:
            assert cur_name not in observed_names, \
                'Duplicate sequence name: %s' % cur_name
            observed_names.append(cur_name)
            return alphabet.parse(cur_seq, name=cur_name)

    count = cur_pos = 0
    cur_name = cur_seq = ''
    # NOTE `for raw_line in f` uses a read-ahead buffer which makes `f.tell()`
    # useless for remembering where a sequence begins.
    # cf. https://docs.python.org/2/library/stdtypes.html#file.next
    for raw_line in iter(f.readline, ''):
        line = raw_line.strip()
        if not line:
            continue
        if line[0] == '>':
            seq = _parse(cur_seq, cur_name)
            if seq:
                yield seq, cur_pos
                count += 1
                if num > 0 and count >= num:
                    return
            cur_pos = f.tell() - len(raw_line)
            cur_name = line[1:].strip()
            cur_seq = ''
        else:
            cur_seq += line
    seq = _parse(cur_seq, cur_name)
    if seq:
        yield seq, cur_pos


def write_fasta(f, seqs, width=80):
    """Writes the given sequence in FASTA format.

    Args:
        f (file): A writable open file handle or anything that responds to
            ``write()``.
        seqs (iterable): An iterable of :class:`NamedSequence
            <biseqt.sequence.NamedSequence>` objects.
    """
    observed_names = []
    for seq in seqs:
        assert isinstance(seq, NamedSequence), \
            'Can only write named sequences'
        name = seq.name if seq.name else seq.content_id

        assert name not in observed_names, 'Duplicate sequence name: %s' % name
        observed_names.append(name)

        contents = str(seq)
        if width is not None:
            contents = '\n'.join(textwrap.wrap(contents, width))
        f.write('>%s\n%s\n' % (name, contents))
        f.flush()


def pw_render_term(tx, origin, mutant, origin_start=0, mutant_start=0,
                   term_width=120, margin=0, colored=True):
    """Renders an edit transcript for display on a terminal.

    Args:
        tx (sequence.EditTranscript): The sequence of edit operations.
        origin (sequence.Sequence): The original :class:`Sequence`.
        mutant (sequence.Sequence): The mutant :class:`Sequence`.

    Keyword Args:
        origin_start (int): Starting position on the original sequence;
            default is 0.
        mutant_start (int): Starting position on the mutant sequence;
            default is 0.
        term_width (int): Terminal width used for wrapping; default is 120
            and the smallest valid value is 30.
        margin (length): Length of leading and trailing substring to
            include in original and mutant sequences; default is 20.
        colored (bool): Whether or not to use ANSI color codes in output;
            default is True.

    Returns:
        str
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
    for op in tx:
        out, carriage = carriage_fwd(carriage, op=op)
        output += out
    output += post_margin(carriage)

    # when output is not supposed to be cleared remove the spurious color
    # reset ANSI escape sequence that termcolor adds:
    if not colored:
        output = output.replace(termcolor.RESET, '')

    return output
