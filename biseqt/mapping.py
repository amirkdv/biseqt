from sys import stderr
from subprocess import Popen
from subprocess import PIPE
from tempfile import NamedTemporaryFile
from collections import namedtuple
from Bio import SeqIO
from re import finditer
from hashlib import sha1
import sqlite3
import sys
import os

from . import seq, words, lib, ProgressIndicator

# FIXME docs
# rc is whether the read is reverse complemented to match the reference
Mapping = namedtuple('Mapping', ['ref', 'rc', 'ref_from', 'ref_to'])

def save_mappings(path, mappings):
    with open(path, 'w') as f:
        f.write('{\n')
        for read, mapping in mappings:
            f.write('%s: %s,\n' % (repr(read), repr(mapping)))
        f.write('\n}\n')

def _identify(path):
    recs = SeqIO.parse(path, 'fasta')
    return {rec.id: sha1(str(rec.seq)).hexdigest() for rec in recs}

def parse_mappings(ref, reads, sam_output, blasr=0):
    proc = Popen(['grep', '^>', reads], stdout=PIPE)
    num_reads = len(proc.communicate()[0].strip().split('\n'))
    indicator = ProgressIndicator(
        'Mapping %d reads' % num_reads, num_reads, percentage=False
    )
    ids = {
        'reads': _identify(reads),
        'refs': _identify(ref),
    }

    indicator.start()
    with open(sam_output) as f:
        for l in f:
            if l[0] == '@':
                continue

            # cf. SAM spec: https://samtools.github.io/hts-specs/SAMv1.pdf
            read, flag, ref, start, quality, cigar = l.strip().split()[:6]
            # FIXME why is blasr changing the ID in such a specfic way?
            if blasr:
                read = read.rsplit('/', 1)[0]
            start, flag = int(start), int(flag)
            if flag == 0x4:
                # Segment unmapped
                continue
            if flag in [0x800, 0x800 + 0x10, 0x100, 0x100 + 0x10]:
                # Secondary/Supplemantary alignment
                continue
            indicator.progress()
            # flag 0x10 means the query (the "segment") has been reveresed
            assert(flag in [0, 0x10])
            len_on_ref = 0
            for match in finditer('\d+(S|M|D|I)', cigar):
                op = match.group()
                assert(op[-1] in 'SMDI')
                length = int(op[:-1])
                # soft clipping
                if op[-1] == 'S' and (match.start() == 0 or match.end() == len(cigar)):
                    # soft clipping at the beginning of query or at the end of
                    # rc'd query.
                    if (match.start() == 0 and flag == 0) or (match.end() == len(cigar) and flag == 0x10):
                        start -= length
                    else:
                        len_on_ref += length
                elif op[-1] in 'MI':
                    len_on_ref += length
                else:
                    assert(op[-1] == 'D')

            end = start + len_on_ref
            yield ids['reads'][read], Mapping(
                ref=ids['refs'][ref],
                rc=1 if flag == 0x10 else 0,
                ref_from=int(start),
                ref_to=end,
            )

        indicator.finish()

