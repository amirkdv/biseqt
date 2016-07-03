# -*- coding: utf-8 -*-
import pytest
from tempfile import NamedTemporaryFile
from sqlite3 import IntegrityError

from biseqt.io import read_fasta, write_fasta
from biseqt.sequence import Alphabet
from biseqt.database import DB, create_record


def test_database_basic():
    A = Alphabet('ACGT')
    with NamedTemporaryFile() as tmp:
        db = DB(tmp.name, A)
        assert db.path == tmp.name
        with db.connect() as conn:
            # a sequence table should be created
            conn.cursor().execute('SELECT * FROM sequence LIMIT 1;')

        # initialization script must be idempotent
        db = DB(tmp.name, A)


def test_database_record():
    A = Alphabet('ACGT')
    with NamedTemporaryFile() as tmp:
        db = DB(tmp.name, A)
        rec = create_record(A.parse('AACT', name='foo'),
                            source_file='source.fa', source_pos=0,
                            attrs={'key': 'value'})
        db.populate([rec])
        retrieved = next(db.find())
        assert isinstance(retrieved[0], int), 'an integer id must be set'
        assert rec[1:] == retrieved[1:], \
            'aside from id we should get what we put in'

        with pytest.raises(IntegrityError):
            db.populate([rec])

        rec = create_record(A.parse('AGCT', name='bar'),
                            source_file='source.fa', source_pos=10)
        db.populate([rec])
        sql_condition = "attrs LIKE '%s'" % '%"name": "bar"%'
        retrieved = [r for r in db.find(sql_condition=sql_condition)]
        assert len(retrieved) == 1 and rec[1:] == retrieved[0][1:], \
            'find() should work with sql_condition'

        def condition(r): return r.attrs['name'] == 'bar'

        retrieved = [r for r in db.find(condition=condition)]
        assert len(retrieved) == 1 and rec[1:] == retrieved[0][1:], \
            'find() should work with callable condition'


def test_database_populate_fasta():
    A = Alphabet('ACGT')
    S = A.parse('AACT', name='S')
    T = A.parse('GCAT', name='T')

    with NamedTemporaryFile() as tmp_db:
        db = DB(tmp_db.name, A)
        with NamedTemporaryFile() as tmp_fa:
            write_fasta(tmp_fa, [S, T])
            tmp_fa.seek(0)
            db.populate_from_fasta(tmp_fa, rc=False)
            retrieved = [x for x in db.find()]
            assert len(retrieved) == 2
            assert all(rec.source_file == tmp_fa.name for rec in retrieved), \
                'source file of sequence records must be set'
            retrieved_T = retrieved[1]
            tmp_fa.seek(retrieved_T.source_pos)
            assert next(read_fasta(tmp_fa, A, num=1))[0] == T, \
                'should be able to retrieve sequences by position in source'

    # with rc
    with NamedTemporaryFile() as tmp_db:
        db = DB(tmp_db.name, A)
        with NamedTemporaryFile() as tmp_fa:
            write_fasta(tmp_fa, [S, T])
            tmp_fa.seek(0)
            db.populate_from_fasta(tmp_fa, rc=True)

            def cond_T(r): return r.attrs['name'] == 'T'

            retrieved_T = next(db.find(condition=cond_T))
            tmp_fa.seek(retrieved_T.source_pos)
            assert next(read_fasta(tmp_fa, A, num=1))[0] == T, \
                'should be able to retrieve sequences by position in source'

            def cond_Trc(r): return r.attrs.get('rc_of', None) == T.content_id

            retrieved_Trc = next(db.find(condition=cond_Trc))
            tmp_fa.seek(retrieved_Trc.source_pos)
            assert next(read_fasta(tmp_fa, A, num=1))[0] == T, \
                'should be able to retrieve sequences by position in source'
            Trc = T.reverse().transform(['AT', 'CG'])
            assert Trc.content_id == retrieved_Trc.content_id, \
                'the stored content_id of reverse complement should be correct'
