# -*- coding: utf-8 -*-
from StringIO import StringIO

from biseqt.io import write_fasta
from biseqt.sequence import Alphabet
from biseqt.database import DB, Record


def test_database_basic():
    A = Alphabet('ACGT')
    db = DB(':memory:', A)
    db.initialize()
    db.initialize()  # should be able to call it twice
    with db.connection() as conn:
        # a sequence table should be created
        conn.cursor().execute('SELECT * FROM sequence LIMIT 1;')


def test_database_insert():
    A = Alphabet('ACGT')
    S = A.parse('AACT', name='foo')
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
    T = A.parse('GCTG', name='bar')
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
    S = A.parse('AACT', name='foo')
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
    S = A.parse('AACT', name='foo')
    T = A.parse('GGCT', name='bar')
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
    S = A.parse('AACT', name='S')
    T = A.parse('GCAT', name='T')

    db = DB(':memory:', A)
    db.initialize()

    fasta = StringIO()
    fasta.name = '/x.fasta'

    write_fasta(fasta, [S, T])
    fasta.seek(0)
    inserted = db.load_fasta(fasta, rc=False)
    assert len(inserted) == 2
    assert all(isinstance(r, Record) for r in inserted)
    assert all(rec.source_file == fasta.name for rec in inserted), \
        'source file of sequence records must be set'
    assert [db.load_from_record(rec, fasta) for rec in inserted] == [S, T], \
        'should be able to retrieve sequences by position in source'


def test_database_populate_fasta_rc():
    A = Alphabet('ACGT')
    S = A.parse('AACT', name='S')
    T = A.parse('GCAT', name='T')

    db = DB(':memory:', A)
    db.initialize()
    fasta = StringIO()
    write_fasta(fasta, [S, T])
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
    S = A.parse('AACT', name='S')

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
