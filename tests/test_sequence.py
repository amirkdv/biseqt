# -*- coding: utf-8 -*-
import pytest
from biseqt.sequence import Alphabet, Sequence, NamedSequence


def test_alphabet():
    # letters of an alphabet must have identical lengths
    with pytest.raises(AssertionError):
        A = Alphabet(['00', '01', '1'])
    A = Alphabet(['00', '01', '10', '11'])
    assert A[0] == '00', 'letters should be available by index'
    # letters of an alphabet are immutable
    with pytest.raises(Exception):
        A[0] = '--'
    # strings as alphabets
    A = Alphabet('ACGT')
    assert A[0] == 'A', 'An alphabet can be created from a single string'


def test_alphabet_magic():
    A = Alphabet('ACGT')

    assert len(A) == 4, 'length of an alphabet is the number of its letters'
    assert A == Alphabet(['A', 'C', 'G', 'T']), 'equal if they same letters'
    assert A != Alphabet('AGCT'), 'not equal if order of letters differ'
    assert A == eval(repr(A)), 'repr() should provide eval-able string'


def test_sequence():
    A = Alphabet('HT')
    S = Sequence(A, (0, 1, 0, 1))

    assert tuple(S) == (0, 1, 0, 1), 'Sequences must support iteration'
    assert S.alphabet == A
    assert str(S) == 'HTHT', 'str() should provide readable representation'
    assert S == eval(repr(S)), 'repr() should provide eval-able string'


def test_sequence_magic():
    A = Alphabet('HT')
    contents = [0, 1, 0, 1]
    S = Sequence(A, contents)

    assert str(S) == 'HTHT'
    assert len(S) == 4, 'len() should work'
    assert S == Sequence(A, contents), 'equals if same contents and alphabet'
    assert S and not Sequence(A, []), 'truthy iff not empty'

    assert S[0] == 0, 'indexing by int should give an int'
    assert isinstance(S[0:1], Sequence) and str(S[0:1]) == 'H', \
        'indexing by a slice should give another sequence object'

    assert S + A.parse('TT') == A.parse('HTHTTT'), 'add by appending'
    assert S + 'TT' == A.parse('HTHTTT'), 'add by appending raw sequences'


def test_sequence_transforms():
    A = Alphabet(['00', '01', '11'])
    S = A.parse('0001')
    assert S.transform(mappings={'00': '01'}) == A.parse('0101'), \
        'dict mappings are unidirectional'
    assert S.transform(mappings=[('00', '01')]) == A.parse('0100'), \
        'list mappings are bidirectional'
    assert S.transform(mappings=[(0, 1)]) == A.parse('0100'), \
        'list mappings are bidirectional'

    S = A.parse('0011')
    assert S.transform(mappings={'00': '01'}) == A.parse('0111'), \
        'unmapped letters remain untouched'
    assert S.reverse() == A.parse('1100'), 'reverse() works'


def test_sequence_parsing():
    A = Alphabet(['00', '01', '10', '11'])
    with pytest.raises(AssertionError):
        A.parse('000')

    S = A.parse('001011')
    assert len(S) == 3 and S == Sequence(A, [0, 2, 3]), \
        'alphabets with > 1 long letters should be able to parse strings'


def test_named_sequence():
    A = Alphabet('ACGT')
    S = A.parse('AACT', name='foo')
    assert isinstance(S, NamedSequence)
    assert eval(repr(S)) == S, 'repr() should provide eval-able string'
    assert S.name == 'foo'
    assert S.content_id == A.parse(str(S), name='bar').content_id, \
        'content id should only depend on the contents of the sequence'
    assert S == A.parse(str(S), name='foo'), \
        'equality should work'


def test_named_sequence_tranforms():
    A = Alphabet('ACGT')
    S = A.parse('AACT', name='foo')
    assert S.reverse(name='bar') == A.parse('TCAA', name='bar'), \
        'reverse of named sequences should be a named sequence'
    complement = S.transform(mappings=['AT', 'CG'], name='bar')
    assert complement == A.parse('TTGA', name='bar'), \
        'result of transforming a named sequence is a named sequence'

    assert 'transformed' in S.transform(mappings=['AT', 'CG']).name
