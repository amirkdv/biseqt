# -*- coding: utf-8 -*-
import textwrap
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
