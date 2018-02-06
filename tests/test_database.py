# -*- coding: utf-8 -*-
import pytest
from tempfile import NamedTemporaryFile
from StringIO import StringIO
from mock import patch, MagicMock

from biseqt.sequence import Alphabet, Sequence
from biseqt.database import DB, NamedSequence, Record, read_fasta


def test_named_sequence():
    A = Alphabet('ACGT')
    S = NamedSequence(A.parse('AACT'), name='foo')
    assert isinstance(S, NamedSequence)
    T = eval(repr(S), {'Alphabet': Alphabet, 'Sequence': Sequence,
                       'NamedSequence': NamedSequence})
    assert T == S, 'repr() should provide eval-able string'
    assert S.name == 'foo'
    T = NamedSequence(A.parse(str(S)), name='bar')
    assert S.content_id == T.content_id, \
        'content id should only depend on the contents of the sequence'
    assert S == NamedSequence(A.parse(str(S)), name='foo'), \
        'equality should work'


def test_named_sequence_tranforms():
    A = Alphabet('ACGT')
    S = NamedSequence(A.parse('AACT'), name='foo')
    T = NamedSequence(A.parse('TCAA'), name='bar')
    assert S.reverse(name='bar') == T, \
        'reverse of named sequences should be a named sequence'
    T = NamedSequence(A.parse('TTGA'), name='bar')
    assert S.transform(mappings=['AT', 'CG'], name='bar') == T, \
        'result of transforming a named sequence is a named sequence'

    assert 'transformed' in S.transform(mappings=['AT', 'CG']).name


def test_database_basic():
    A = Alphabet('ACGT')
    db = DB(':memory:', A)
    db.initialize()
    db.initialize()  # should be able to call it twice
    with db.connection() as conn:
        # a sequence table should be created
        conn.cursor().execute('SELECT * FROM sequence LIMIT 1;')

    with pytest.raises(AssertionError):
        DB('/cannot/possibly/exist/directory/', A)


def test_database_insert():
    A = Alphabet('ACGT')
    S = NamedSequence(A.parse('AACT'), name='foo')
    db = DB(':memory:', A)
    db.initialize()
    attrs = {'key': 'value'}
    rec = db.insert(S, source_file='source.fa', source_pos=10, attrs=attrs)
    assert isinstance(rec.id, int)
    assert rec.content_id == S.content_id
    assert rec.source_pos == 10
    assert rec.source_file == 'source.fa'
    assert 'key' in rec.attrs and rec.attrs['key'] == 'value', \
        'attributes must be populated correctly'
    with db.connection() as conn:
        cursor = conn.cursor()
        cursor.execute('SELECT content_id FROM sequence WHERE id = ?',
                       (rec.id,))
        # NOTE for some reason if we just say next(cursor) ==  ...
        # the cursor remains open after the context is over (which should
        # not happen as per docs). This leads to BusyError further down.
        assert cursor.fetchall() == [(S.content_id,)], \
            'content identifier is properly populated'

    # add a second sequence
    T = NamedSequence(A.parse('GCTG'), name='bar')
    new_rec = db.insert(T)
    assert new_rec.id != rec.id, 'new ids are assigned to new sequences'
    with db.connection() as conn:
        cursor = conn.cursor()
        cursor.execute('SELECT content_id FROM sequence WHERE id = ?',
                       (new_rec.id,))
        assert next(cursor) == (T.content_id,), \
            'correct id must be populated'


def test_database_overwrite():
    A = Alphabet('ACGT')
    S = NamedSequence(A.parse('AACT'), name='foo')
    db = DB(':memory:', A)
    db.initialize()
    db.insert(S, source_file='old_source.fa')
    db.insert(S, source_file='new_source.fa')
    with db.connection() as conn:
        cursor = conn.cursor()
        cursor.execute(
            'SELECT source_file FROM sequence WHERE content_id = ?',
            (S.content_id,)
        )
        res = [x[0] for x in cursor]
        assert len(res) == 1 and res[0] == 'old_source.fa', \
            'Sequences with observed content id should be ignored'


def test_database_find():
    A = Alphabet('ACGT')
    S = NamedSequence(A.parse('AACT'), name='foo')
    T = NamedSequence(A.parse('GGCT'), name='bar')
    db = DB(':memory:', A)
    db.initialize()
    db.insert(S)
    db.insert(T)

    sql_condition = "attrs LIKE '%s'" % '%"name": "bar"%'
    found = [rec for rec in db.find(sql_condition=sql_condition)]
    assert len(found) == 1 and found[0].content_id == T.content_id, \
        'find() should work with sql_condition'

    def condition(rec): return rec.attrs['name'] == 'foo'

    found = [rec for rec in db.find(condition=condition)]
    assert len(found) == 1 and found[0].content_id == S.content_id, \
        'find() should work with callable condition'


def test_database_populate_fasta():
    A = Alphabet('ACGT')
    S = NamedSequence(A.parse('AACT'), name='S')
    T = NamedSequence(A.parse('GCAT'), name='T')

    db = DB(':memory:', A)
    db.initialize()

    fasta = StringIO()
    fasta.name = '/x.fasta'

    fasta.write('\n'.join(x.to_fasta() for x in [S, T]))
    fasta.seek(0)
    inserted = db.load_fasta(fasta, rc=False)
    assert len(inserted) == 2
    assert all(isinstance(r, Record) for r in inserted)
    assert all(rec.source_file == fasta.name for rec in inserted), \
        'source file of sequence records must be set'
    assert [db.load_from_record(rec, fasta) for rec in inserted] == [S, T], \
        'should be able to retrieve sequences by position in source'

    with patch('biseqt.database.open', create=True) as open_mock:
        open_mock.return_value = MagicMock(spec=file, wraps=fasta)
        assert db.load_from_record(inserted[0]) == S, \
            'load_from_record should work without an open file handle'
        # the same cannot be repeated for T because load_from_record calls
        # close on the fake file. If additional tests are needed wrap "fasta"
        # twice: MagicMock(spec=file, wraps=StringIO(fasta.getvalue()))


def test_database_populate_fasta_rc():
    A = Alphabet('ACGT')
    S = NamedSequence(A.parse('AACT'), name='S')
    T = NamedSequence(A.parse('GCAT'), name='T')

    db = DB(':memory:', A)
    db.initialize()
    fasta = StringIO()
    fasta.write('\n'.join(x.to_fasta() for x in [S, T]))
    fasta.seek(0)
    inserted = db.load_fasta(fasta, rc=True)

    assert len(inserted) == 4
    assert [r.attrs['rc_of'] for r in inserted if 'rc_of' in r.attrs] \
        == [S.content_id, T.content_id], \
        'reverse complements should know what their origin is'

    def cond_T_rc(r): return r.attrs.get('rc_of', None) == T.content_id

    found_T_rc = next(db.find(condition=cond_T_rc))
    T_rc = T.reverse().transform(['AT', 'CG'], name='(rc) ' + T.name)
    assert db.load_from_record(found_T_rc, fasta) == T_rc, \
        'reverse complements should load properly from a record'


def test_database_events():
    A = Alphabet('ACGT')
    S = NamedSequence(A.parse('AACT'), name='S')

    # NOTE python 2 does not support non-local, non-global variables, put it in
    # the function object.
    test_database_events.callback_called = 0

    def callback(self, *args):
        test_database_events.callback_called += 1

    db = DB(':memory:', A)
    db.add_event_listener('db-initialized', callback)
    db.add_event_listener('sequence-inserted', callback)
    db.initialize()
    assert test_database_events.callback_called == 1, \
        'event callbacks for "initialize" should be executed'

    db.insert(S)
    assert test_database_events.callback_called == 2, \
        'event callbacks for "insert-sequence" should be executed'


def test_read_fasta_basic():
    A = Alphabet('ACGT')
    with NamedTemporaryFile() as f:
        f.write('> name1\nAAA\n\nTTT')
        f.flush()
        f.seek(0)
        recs = [r for r in read_fasta(f, A)]
        assert len(recs) == 1, 'should work when reading from file'
        assert recs[0][0] == NamedSequence(A.parse('AAATTT'), name='name1'), \
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


def test_to_fasta():
    A = Alphabet('ACGT')
    S = NamedSequence(A.parse('AAA'), name='foo')
    T = NamedSequence(A.parse('TTT'), name='bar')

    with NamedTemporaryFile() as f:
        f.write('\n'.join(x.to_fasta() for x in [S, T]))
        f.seek(0)
        assert [s for s, _ in read_fasta(f, A)] == [S, T], \
            'read_fasta should be the opposite of to_fasta'

    S = NamedSequence(A.parse('AAATTT'), name='foo')
    f = StringIO('')
    f.write(S.to_fasta(width=3))
    f.seek(0)
    assert sum(1 for _ in f) == 3, 'FASTA width should be modifiable'
