#!/usr/bin/env python

import os.path
import sys
import tempfile
from Bio import SeqIO, SeqRecord
import subprocess
from hashlib import sha1
from align import ProgressIndicator

BWA_IDX = os.environ['DB']
READS_IN = os.environ['READS']
READS_OUT = sys.argv[1]
NUM_READS = int(os.environ.get('NUM_READS', -1))
annotated = []
indicator = ProgressIndicator('Annotating first %d reads' % NUM_READS, NUM_READS, percentage=False)
indicator.start()
mismatch_stats = (0, 0) # number of observed numbers, the average
gap_stats = (0,0)

for rec in SeqIO.parse(READS_IN, 'fasta'):
    if NUM_READS > 0 and indicator.progress_cnt >= NUM_READS:
        break
    with tempfile.NamedTemporaryFile(suffix='.fasta') as tmp:
        tmp.write('>tmp\n' + str(rec.seq))
        tmp.flush() # flush to disk, a subprocess will be accessing it

        args = ['bwa', 'mem', BWA_IDX, tmp.name]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, _ = proc.communicate()

        indicator.progress()
        out = out.strip().split('\n')[1:]
        _, flag, _, pos, quality = out[0].split()[:5]
        # FIXME document quality of reads in their names
        flag, pos = int(flag), int(pos)
        if flag == 4:
            skip += 1
            # Segment unmapped
            continue
        assert(flag in [0,16])
        annotated += [SeqRecord.SeqRecord(
            rec.seq if flag == 0 else rec.seq.reverse_complement(),
            id='R%s_P%d' % (sha1(str(rec.seq)).hexdigest()[:8], pos),
            description=rec.id)
        ]

indicator.finish()
SeqIO.write(annotated, READS_OUT, 'fasta')
