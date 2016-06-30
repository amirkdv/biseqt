# -*- coding: utf-8 -*-
import pytest
from biseqt.sequence import Alphabet, Sequence, EditTranscript


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


def test_sequence_parsing():
    A = Alphabet(['00', '01', '10', '11'])
    with pytest.raises(AssertionError):
        A.parse('000')

    S = A.parse('001011')
    assert len(S) == 3 and S == Sequence(A, [0, 2, 3]), \
        'alphabets with > 1 long letters should be able to parse strings'


def test_transcript():
    with pytest.raises(AssertionError):
        EditTranscript('T')

    tx = EditTranscript('MM')
    assert str(tx) == 'MM', 'str() gives the raw opseq'
    assert len(tx) == 2, 'len() works'
    assert eval(repr(tx)) == tx, 'repr() should provide eval-able string'
    assert tx == EditTranscript('MM'), 'equal opseq means equal transcripts'

    assert tx + EditTranscript('S') == EditTranscript('MMS'), \
        'transcript + transcript gives a transcript'
    assert tx + 'S' == EditTranscript('MMS'), \
        'transcript + string gives transcript'

    assert tx[0] == 'M', 'indexing by int should give a string'
    assert tx[:1] == EditTranscript('M'), \
        'indexing by slice should give another transcript'


def test_transcript_render_basic():
    A = Alphabet('ACGT')
    S = A.parse('AACT')
    tx = EditTranscript('M' * len(S))
    origin = S + S
    mutant = S + S
    render_args = origin, mutant
    render_kw = {'origin_start': len(S), 'colored': False}

    assert tx.render_term(S, S, colored=False).count('\033') == 0, \
        'colored output should allow being turned off'

    assert tx.render_term(S, S, colored=True).count('\033') > 0, \
        'colored output should allow being turned on'

    # validate input
    with pytest.raises(AssertionError):
        tx.render_term(S, S, margin=-1)
    with pytest.raises(AssertionError):
        tx.render_term(S, S, term_width=5)
    with pytest.raises(AssertionError):
        tx.render_term(S, S, origin_start=-1)
    with pytest.raises(AssertionError):
        tx.render_term(S, S, mutant_start=4)

    no_margin = tx.render_term(*render_args, margin=0, **render_kw)
    assert '[%d]' % len(S) in no_margin, 'margin should allow being turned off'

    with_margin = tx.render_term(*render_args, margin=1, **render_kw)
    assert '[%d]' % (len(S) - 1) in with_margin, \
        'margin should allow being turned on'

    # shouldn't choke on too large margins
    full_margin = tx.render_term(*render_args, margin=30, **render_kw)
    assert str(S) + '.' * len(S) in full_margin, 'overhanging margins work'
    assert len(set(len(l) for l in full_margin.rstrip().split('\n'))) == 1, \
        'both lines of the output should have the same length'

    # deletion:
    #   AACT
    #   AG-T
    tx = EditTranscript('MSDM')
    with_del = tx.render_term(S + S, A.parse('AGT'), **render_kw)
    assert 'AG-T' in with_del, 'deletions are represented by - in mutant'
    lines = with_del.rstrip().split('\n')
    assert lines[0].index('C') == lines[1].index('-'), \
        'deleted content and - should be aligned'

    # insertion:
    #   AAC-T
    #   AACGT
    tx = EditTranscript('MMMIM')
    with_ins = tx.render_term(S + S, A.parse('AACGT'), **render_kw)
    assert 'AAC-T' in with_ins, 'insertions are represented by - in origin'
    lines = with_ins.rstrip().split('\n')
    assert lines[0].index('-') == lines[1].index('G'), \
        'inserted content and - should be aligned'


def test_transcript_render_width():
    A = Alphabet('ACGT')
    N = 100
    long_tx = EditTranscript('M' * N)
    origin = mutant = A.parse('A' * (2 * N))
    term_width = N / 2
    render = long_tx.render_term(origin, mutant, margin=2*N, colored=False,
                                 term_width=term_width, origin_start=N)
    line_lens = [len(l) for l in render.rstrip().split('\n')]
    assert all(length <= term_width for length in line_lens), \
        'terminal width should be adjustable'
    assert any(length == term_width for length in line_lens), \
        'terminal width should be fully used'
    assert len(set(line_lens)) <= 2, \
        'alignments longer than terminal width should work'


def test_transcript_render_longlet():
    A = Alphabet(['00', '11'])
    origin = A.parse('0011')
    mutant = A.parse('11')
    tx = EditTranscript('DM')
    assert '--11' in tx.render_term(origin, mutant, colored=False), \
        'alphabets with > 1 long letters should be rendered properly'
