# -*- coding: utf-8 -*-
import pytest
from tempfile import NamedTemporaryFile
from StringIO import StringIO

from biseqt.sequence import Alphabet, NamedSequence, EditTranscript
from biseqt.io import pw_render_term, read_fasta, write_fasta


def test_read_fasta_basic():
    A = Alphabet('ACGT')
    with NamedTemporaryFile() as f:
        f.write('> name1\nAAA\n\nTTT')
        f.flush()
        f.seek(0)
        recs = [r for r in read_fasta(f, A)]
        assert len(recs) == 1, 'should work when reading from file'
        assert recs[0][0] == A.parse('AAATTT', name='name1'), \
            'should properly parse what is in the file'
        assert isinstance(recs[0][0], NamedSequence), \
            'should return NamedSequence objects'
        assert recs[0][1] == 0, 'should report the right file positions'

    # duplicate names not allowed
    with pytest.raises(AssertionError):
        [r for r in read_fasta(StringIO('>name\nAAA\n> name\nTTT\n'), A)]


def test_read_fasta_advanced():
    A = Alphabet('ACGT')
    fasta = '>name0\nAAA\n\n>name1\nTTT\n\n>name2\nCCC\n>name3\nGGG'
    f = StringIO(fasta)

    recs = [r for r in read_fasta(f, A)]
    assert len(recs) == fasta.count('>'), 'should work with StringIO'

    # check loading from position
    for idx in range(4):
        assert recs[idx][1] == fasta.index('>name%d' % idx), \
            'should report the right file positions'
        f.seek(recs[idx][1])
        assert next(read_fasta(f, A, num=1))[0] == recs[idx][0], \
            'should be able to read single sequences from known positions'

    f.seek(recs[1][1])
    assert [r for r in read_fasta(f, A, num=2)] == recs[1:3], \
        'should be able to read known number of sequences from known positions'


def test_write_fasta():
    A = Alphabet('ACGT')
    S = A.parse('AAA', name='foo')
    T = A.parse('TTT', name='bar')

    with NamedTemporaryFile() as f:
        write_fasta(f, [S, T])
        f.seek(0)
        assert [s for s, _ in read_fasta(f, A)] == [S, T], \
            'read_fasta(write_fasta()) should be identity'

    f = StringIO('')
    write_fasta(f, [S, T])
    f.seek(0)
    assert f.read() == '>foo\nAAA\n>bar\nTTT\n', 'should work on StringIO'

    f = StringIO('')
    # duplicate names not allowed
    with pytest.raises(AssertionError):
        write_fasta(f, [S, S])

    f = StringIO('')
    S = A.parse('AAATTT', name='foo')
    write_fasta(f, [S], width=3)  # should take 3 lines
    f.seek(0)
    assert sum(1 for _ in f) == 3, 'FASTA width should be modifiable'


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
    origin, mutant = S + S, A.parse('AGT')
    with_del = pw_render_term(tx, origin, mutant, **render_kw)
    assert 'AG-T' in with_del, 'deletions are represented by - in mutant'
    lines = with_del.rstrip().split('\n')
    assert lines[0].index('C') == lines[1].index('-'), \
        'deleted content and - should be aligned'
    # shouldn't crash when printing with color
    pw_render_term(tx, origin, mutant, origin_start=len(S), colored=True)

    # insertion:
    #   AAC-T
    #   AACGT
    tx = EditTranscript('MMMIM')
    origin, mutant = S + S, A.parse('AACGT')
    with_ins = pw_render_term(tx, origin, mutant, **render_kw)
    assert 'AAC-T' in with_ins, 'insertions are represented by - in origin'
    lines = with_ins.rstrip().split('\n')
    assert lines[0].index('-') == lines[1].index('G'), \
        'inserted content and - should be aligned'
    # shouldn't crash when printing with color
    pw_render_term(tx, origin, mutant, origin_start=len(S), colored=True)


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
