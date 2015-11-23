#!/usr/bin/env python

import os.path
import sys
import tempfile
from Bio import SeqIO, SeqRecord
import subprocess
from uuid import uuid4

HERE = os.path.dirname(__file__)
BLAST_DB = os.environ['DB']
READS_IN = os.environ['READS']
READS_OUT = sys.argv[1]
NUM_READS = 300
# NUM_READS = 1
annotated = []

count = 0
for rec in SeqIO.parse(READS_IN, 'fasta'):
    if count == NUM_READS:
        break
    with tempfile.NamedTemporaryFile(suffix='.fasta') as tmp:
        tmp.write('>tmp\n' + str(rec.seq))
        args = ['blastn', '-db', BLAST_DB, '-query', tmp.name, '-outfmt', '7 qstart sstart']
        proc = subprocess.Popen(args, stdout=subprocess.PIPE)
        out, err = proc.communicate()

        # print out
        # count += 1
        # continue

        tokens = [l for l in out.split('\n') if l and l[0] != '#']
        if not tokens:
            continue
        count += 1
        qstart, start = [int(t) for t in tokens[0].split()]
        seqid = 'R%s_P%d' % (str(uuid4())[:8], start - qstart)
        print seqid
        annotated += [SeqRecord.SeqRecord(rec.seq, id=seqid, description=rec.id)]

SeqIO.write(annotated, READS_OUT, 'fasta')
