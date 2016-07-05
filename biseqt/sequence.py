# -*- coding: utf-8 -*-
from hashlib import sha1
from itertools import chain


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

    def letter_to_idx(self, letters):
        """Translates provided letters to the integer sequence corresponding
        to the index of each letter in this alphabet.

        Args:
            letters (iterable): The letters to be translated to integer
                indices. Each element retrieved through iteration should be
                an element in :attr:`_letters`.

        Returns:
            tuple
        """
        return tuple(self._idx_by_letter[l] for l in letters)

    def parse(self, string, name=None):
        """Given a string representation of a sequence returns a corresponding
        :class:`Sequence` object.

        Args:
            string (str): The raw sequence represented as a string.
            name (str): The name for the sequence; default is None in which
                case a :class:`NamedSequence` will be returned.

        Returns:
            Sequence
        """
        assert isinstance(string, str), 'Raw sequence must be in string form'
        assert len(string) % self._letlen == 0, 'String representation ' + \
            'of sequence must be a multiple of the alphabet letter length'
        contents = []
        idx = 0
        while idx < len(string):
            contents.append(string[idx:idx + self._letlen])
            idx += self._letlen
        contents = self.letter_to_idx(contents)
        if name is None:
            return Sequence(self, contents)
        else:
            return NamedSequence(self, contents, name=name)

    def transform(self, seq, mappings={}):
        """Transforms the given sequence to another sequence in the same
        alphabet according to provided letter-to-letter mappings.

        Args:
            seq (Sequence): The original sequence.
            mappings (list|dict): If a dictionary is given, each entry
                represents a translation rule from the key to the value. If a
                list is given, each entry must have two elements and is taken
                to represent a bidirectional translation rule between those two
                elements. Each element in either a dictionary or a list can
                either be a letter in string format or an integer representing
                the position of the letter.
        Returns:
            Sequence

        For example, to get the complement of a DNA sequence::

            >>> from biseqt.sequence import Alphabet, complement
            >>> A = Alphabet('ACGT')
            >>> S = A.parse('AGGGT')
            >>> print A.transform(S, mappings=['AT', 'CG'])
            'TCCCA'

        whereas to get the same effect with a dictionary::

            >>> mappings = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            >>> print A.complement(S, mappings)
            'TCCCA'

        """
        if mappings is None:
            mappings = {}

        if isinstance(mappings, list):
            assert all(len(m) == 2 for m in mappings)
            mappings = dict(chain.from_iterable(
                [(rule[0], rule[1]), (rule[1], rule[0])] for rule in mappings
            ))

        # for any letter (as int) c, pair_of[c] determines what its mapped to.
        pair_of = range(len(self))  # do not modify by default

        for key, val in mappings.items():
            if not isinstance(key, int):
                key = self._idx_by_letter[key]
            if not isinstance(val, int):
                val = self._idx_by_letter[val]
            pair_of[key] = val

        assert all(isinstance(idx, int) for idx in pair_of)
        return Sequence(self, tuple(pair_of[c] for c in seq))

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

    def reverse(self):
        """Returns another sequence whose contents are the reverse of this
        sequence in order.

        Returns:
            Sequence
        """
        return Sequence(self.alphabet, tuple(reversed(self.contents)))

    def transform(self, mappings={}):
        """Wraps :func:`Alphabet.transform` for convenience."""
        return self.alphabet.transform(self, mappings=mappings)

    def to_named(self, name):
        """Names this sequence, i.e a :class:`NamedSequence` is returned with
        identical raw contents.

        Args:
            name (str): The name to give the sequence.

        Returns:
            NamedSequence
        """
        return NamedSequence(self.alphabet, self.contents, name=name)

    def __str__(self):
        return ''.join(self.alphabet[idx] for idx in self.contents)

    def __repr__(self):
        return 'Sequence(%s, contents=%s)' % \
            (repr(self.alphabet), repr(self.contents))

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
            other = self.alphabet.letter_to_idx(other)
        return Sequence(self.alphabet, self.contents + contents)


class NamedSequence(Sequence):
    """A named version of :class:`Sequence`. The main differences are the
    additional attributes which allow for faster equality comparison for large
    sequences.

    Attributes:
        name (str): The name of the sequence which need not be unique.
        content_id (str): The SHA1 of the sequence contents in integer form.
    """
    def __init__(self, alphabet, contents=(), name='', content_id=None):
        super(NamedSequence, self).__init__(alphabet, contents)
        self.content_id = sha1(str(self)).hexdigest()
        if content_id is not None:
            assert self.content_id == content_id, \
                'Provided content identifier does not match sequence contents'
        self.name = name

    def reverse(self, name=None):
        """Wraps :func:`Sequence.reverse` to make sure a named sequence is
        returned.

        Args:
            name(str): The name to give to the new sequence. Default is None
                in which case ``(reversed)`` is preprended to the original
                sequence name.

        Returns:
            NamedSequence
        """
        rev = super(NamedSequence, self).reverse()
        if name is None:
            name = '(reversed) ' + self.name
        return NamedSequence(rev.alphabet, rev.contents, name=name)

    def transform(self, mappings={}, name=None):
        """Wraps :func:`Sequence.transform` to make sure a named sequence is
        returned.

        Args:
            mappings (dict|list): As in :func:`Sequence.transform`.
            name (str): The name to give to the new sequence. Default is None
                in which case ``(transformed)`` is prependended to the
                original sequence name.

        Returns:
            NamedSequence
        """
        seq = super(NamedSequence, self).transform(mappings=mappings)
        if name is None:
            name = '(transformed) ' + self.name
        return NamedSequence(seq.alphabet, seq.contents, name=name)

    def __eq__(self, other):
        return isinstance(other, NamedSequence) and \
               self.content_id == other.content_id and \
               self.name == other.name

    def __repr__(self):
        return 'NamedSequence(%s, name=%s, contents=%s, content_id=%s)' % (
            repr(self.alphabet),
            repr(self.name),
            repr(self.contents),
            repr(self.content_id),
        )


class EditTranscript(str):
    """Represents an edit transcript without reference to the original or
    mutant sequences."""
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
