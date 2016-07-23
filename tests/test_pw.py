# -*- coding: utf-8 -*-
import pytest
from biseqt.sequence import Alphabet
from biseqt.pw import Alignment


def test_projected_aln_len():
    assert Alignment.projected_len('MMM', on='origin') == 3
    assert Alignment.projected_len('MMM', on='mutant') == 3
    assert Alignment.projected_len('SMS', on='origin') == 3
    assert Alignment.projected_len('SMS', on='mutant') == 3
    assert Alignment.projected_len('DMS', on='origin') == 3
    assert Alignment.projected_len('DMS', on='mutant') == 2
    assert Alignment.projected_len('IMS', on='origin') == 2
    assert Alignment.projected_len('IMS', on='mutant') == 3


def test_alignment_basic():
    A = Alphabet('ACGT')
    S = A.parse('AACCGGTT')
    with pytest.raises(AssertionError):
        Alignment(S, S, 'MSSST')  # illegal character
    with pytest.raises(AssertionError):
        Alignment(S, S, 'M', origin_start=len(S))  # illegal starting point
    with pytest.raises(AssertionError):
        Alignment(S, S, 'MM', origin_start=len(S)-1)  # transcript too long


def test_pw_render_basic():
    A = Alphabet('ACGT')
    S = A.parse('AACT')
    aln = Alignment(S, S, 'M' * len(S))
    assert aln.render_term(colored=False).count('\033') == 0, \
        'colored output should allow being turned off'
    assert aln.render_term(colored=True).count('\033') > 0, \
        'colored output should allow being turned on'
    # validate input
    with pytest.raises(AssertionError):
        aln.render_term(margin=-1)
    with pytest.raises(AssertionError):
        aln.render_term(term_width=5)

    aln = Alignment(S + S, S + S, 'M' * len(S), origin_start=len(S))
    no_margin = aln.render_term(margin=0, colored=False)
    assert '[%d]' % len(S) in no_margin, 'margin should allow being turned off'

    with_margin = aln.render_term(margin=1, colored=False)
    assert '[%d]' % (len(S) - 1) in with_margin, \
        'margin should allow being turned on'

    # shouldn't choke on too large margins
    full_margin = aln.render_term(margin=30, colored=False)
    assert str(S) + '.' * len(S) in full_margin, 'overhanging margins work'
    assert len(set(len(l) for l in full_margin.rstrip().split('\n'))) == 1, \
        'both lines of the output should have the same length'

    # deletion:
    #   AACT
    #   AG-T
    aln = Alignment(S + S, A.parse('AGT'), 'MSDM', origin_start=len(S))
    with_del = aln.render_term(colored=False)
    assert 'AG-T' in with_del, 'deletions are represented by - in mutant'
    lines = with_del.rstrip().split('\n')
    assert lines[0].index('C') == lines[1].index('-'), \
        'deleted content and - should be aligned'
    # shouldn't crash when printing deletions with color
    aln.render_term(colored=True)

    # insertion:
    #   AAC-T
    #   AACGT
    aln = Alignment(S + S, A.parse('AACGT'), 'MMMIM', origin_start=len(S))
    with_ins = aln.render_term(colored=False)
    assert 'AAC-T' in with_ins, 'insertions are represented by - in origin'
    lines = with_ins.rstrip().split('\n')
    assert lines[0].index('-') == lines[1].index('G'), \
        'inserted content and - should be aligned'
    # shouldn't crash when printing with color
    with_ins = aln.render_term(colored=True)


def test_pw_render_width():
    A = Alphabet('ACGT')
    N = 100
    S = A.parse('A' * (2 * N))
    tx, term_width = 'M' * N, N/2
    aln = Alignment(S, S, tx, origin_start=N)
    render = aln.render_term(margin=2*N, colored=False, term_width=term_width)
    line_lens = [len(l) for l in render.rstrip().split('\n')]
    assert all(length <= term_width for length in line_lens), \
        'terminal width should be adjustable'
    assert any(length == term_width for length in line_lens), \
        'terminal width should be fully used'
    assert len(set(line_lens)) <= 2, \
        'alignments longer than terminal width should work'


def test_pw_render_longlet():
    A = Alphabet(['00', '11'])
    aln = Alignment(A.parse('0011'), A.parse('11'), 'DM')
    assert '--11' in aln.render_term(colored=False), \
        'alphabets with > 1 long letters should be rendered properly'
