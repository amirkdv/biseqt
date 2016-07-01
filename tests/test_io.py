# -*- coding: utf-8 -*-
import pytest
from tempfile import NamedTemporaryFile

from biseqt.sequence import Alphabet, NamedSequence, EditTranscript
from biseqt.io import pw_render_term, load_fasta


def test_load_fasta():
    A = Alphabet('ACGT')
    with NamedTemporaryFile() as tmp:
        tmp.write('> name1\nAAA\n\nTTT')
        tmp.flush()
        tmp.seek(0)
        seqs = [s for s in load_fasta(tmp, A)]
        assert len(seqs) == 1
        assert seqs[0] == A.parse('AAATTT', name='name1')
        assert isinstance(seqs[0], NamedSequence)

    with NamedTemporaryFile() as tmp:
        tmp.write('>name1\nAAA\n\n>name2\nTTT\n')
        tmp.flush()
        tmp.seek(0)

        seqs = [s for s in load_fasta(tmp, A)]
        assert len(seqs) == 2
        assert seqs[1] == A.parse('TTT', name='name2')

    with NamedTemporaryFile() as tmp:
        tmp.write('>name2\nAAA\n> name2\nTTT\n')
        tmp.flush()
        tmp.seek(0)
        with pytest.raises(AssertionError):
            [s for s in load_fasta(tmp, A)]  # duplicate names not allowed


def test_pw_render_basic():
    A = Alphabet('ACGT')
    S = A.parse('AACT')
    tx = EditTranscript('M' * len(S))
    origin = S + S
    mutant = S + S
    render_args = origin, mutant
    render_kw = {'origin_start': len(S), 'colored': False}

    assert pw_render_term(tx, S, S, colored=False).count('\033') == 0, \
        'colored output should allow being turned off'

    assert pw_render_term(tx, S, S, colored=True).count('\033') > 0, \
        'colored output should allow being turned on'

    # validate input
    with pytest.raises(AssertionError):
        pw_render_term(tx, S, S, margin=-1)
    with pytest.raises(AssertionError):
        pw_render_term(tx, S, S, term_width=5)
    with pytest.raises(AssertionError):
        pw_render_term(tx, S, S, origin_start=-1)
    with pytest.raises(AssertionError):
        pw_render_term(tx, S, S, mutant_start=4)

    no_margin = pw_render_term(tx, *render_args, margin=0, **render_kw)
    assert '[%d]' % len(S) in no_margin, 'margin should allow being turned off'

    with_margin = pw_render_term(tx, *render_args, margin=1, **render_kw)
    assert '[%d]' % (len(S) - 1) in with_margin, \
        'margin should allow being turned on'

    # shouldn't choke on too large margins
    full_margin = pw_render_term(tx, *render_args, margin=30, **render_kw)
    assert str(S) + '.' * len(S) in full_margin, 'overhanging margins work'
    assert len(set(len(l) for l in full_margin.rstrip().split('\n'))) == 1, \
        'both lines of the output should have the same length'

    # deletion:
    #   AACT
    #   AG-T
    tx = EditTranscript('MSDM')
    with_del = pw_render_term(tx, S + S, A.parse('AGT'), **render_kw)
    assert 'AG-T' in with_del, 'deletions are represented by - in mutant'
    lines = with_del.rstrip().split('\n')
    assert lines[0].index('C') == lines[1].index('-'), \
        'deleted content and - should be aligned'

    # insertion:
    #   AAC-T
    #   AACGT
    tx = EditTranscript('MMMIM')
    with_ins = pw_render_term(tx, S + S, A.parse('AACGT'), **render_kw)
    assert 'AAC-T' in with_ins, 'insertions are represented by - in origin'
    lines = with_ins.rstrip().split('\n')
    assert lines[0].index('-') == lines[1].index('G'), \
        'inserted content and - should be aligned'


def test_pw_render_width():
    A = Alphabet('ACGT')
    N = 100
    long_tx = EditTranscript('M' * N)
    origin = mutant = A.parse('A' * (2 * N))
    term_width = N / 2
    render = pw_render_term(long_tx, origin, mutant, margin=2*N, colored=False,
                            term_width=term_width, origin_start=N)
    line_lens = [len(l) for l in render.rstrip().split('\n')]
    assert all(length <= term_width for length in line_lens), \
        'terminal width should be adjustable'
    assert any(length == term_width for length in line_lens), \
        'terminal width should be fully used'
    assert len(set(line_lens)) <= 2, \
        'alignments longer than terminal width should work'


def test_pw_render_longlet():
    A = Alphabet(['00', '11'])
    origin = A.parse('0011')
    mutant = A.parse('11')
    tx = EditTranscript('DM')
    assert '--11' in pw_render_term(tx, origin, mutant, colored=False), \
        'alphabets with > 1 long letters should be rendered properly'
