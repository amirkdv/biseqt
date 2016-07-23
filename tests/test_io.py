# -*- coding: utf-8 -*-
import pytest
from tempfile import NamedTemporaryFile
from StringIO import StringIO

from biseqt.sequence import Alphabet, NamedSequence
from biseqt.io import read_fasta, write_fasta


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
