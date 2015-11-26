#!/usr/bin/env python

import os.path
import sys
import tempfile
from Bio import SeqIO, SeqRecord
import subprocess
from uuid import uuid4
from align import ProgressIndicator

HERE = os.path.dirname(__file__)
BLAST_DB = os.environ['DB']
READS_IN = os.environ['READS']
READS_OUT = sys.argv[1]
NUM_READS = int(os.environ.get('NUM_READS', -1))
annotated = []
indicator = ProgressIndicator('Finding correct positition of reads', NUM_READS, percentage=False)
indicator.start()
mismatch_stats = (0, 0) # number of observed numbers, the average
gap_stats = (0,0)
for rec in SeqIO.parse(READS_IN, 'fasta'):
    if NUM_READS > 0 and indicator.progress_cnt >= NUM_READS:
        break
    with tempfile.NamedTemporaryFile(suffix='.fasta') as tmp:
        tmp.write('>tmp\n' + str(rec.seq))
        args = ['blastn', '-db', BLAST_DB, '-query', tmp.name, '-outfmt', '7 qstart sstart length mismatch gaps']
        proc = subprocess.Popen(args, stdout=subprocess.PIPE)
        out, _ = proc.communicate()
        tokens = [l for l in out.split('\n') if l and l[0] != '#']
        if not tokens:
            continue
        indicator.progress()
        qstart, start, length, mismatch, gaps = [int(t) for t in tokens[0].split()]
        mismatch_stats = (
            mismatch_stats[0] + 1,
            (mismatch_stats[0]*mismatch_stats[1] + mismatch*1.0/length)/(mismatch_stats[0]+1)
        )
        gap_stats = (
            gap_stats[0] + 1,
            (gap_stats[0]*gap_stats[1] + gaps*1.0/length
            )/(gap_stats[0]+1)
        )
        seqid = 'R%s_P%d' % (str(uuid4())[:8], start - qstart)
        annotated += [SeqRecord.SeqRecord(rec.seq, id=seqid, description=rec.id)]

indicator.finish()
print 'Average substitution rate: %.2f%%' % (mismatch_stats[1] * 100)
print 'Average indel rate: %.2f%%' % (gap_stats[1] * 100)
SeqIO.write(annotated, READS_OUT, 'fasta')
