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

from . import seq, words, lib, ProgressIndicator, overlap

class ReadMapper(object):
    def __init__(self, **kw):
        self.read_src = kw['read_src']
        self.ref_src = kw['ref_src']

    def identify(self, path):
        recs = SeqIO.parse(path, 'fasta')
        return {rec.id: sha1(str(rec.seq)).hexdigest() for rec in recs}

    def num_reads(self):
        proc = Popen(['grep', '^>', self.read_src], stdout=PIPE)
        return len(proc.communicate()[0].strip().split('\n'))

    def indicator(self):
        n = self.num_reads()
        indicator = ProgressIndicator(
            'Mapping %d reads' % n, n, percentage=False
        )
        return indicator

    def mappings(self):
        raise NotImplementedError

# FIXME docs
Mapping = namedtuple('Mapping', ['ref', 'strand', 'ref_from', 'ref_to'])

def save_mappings(path, mappings):
    with open(path, 'w') as f:
        f.write('{\n')
        for read, mapping in mappings:
            f.write('%s: %s,\n' % (repr(read), repr(mapping)))
        f.write('\n}\n')

class BwaReadMapper(ReadMapper):
    def index(self):
        sys.stderr.write('Creating BWA index for %s.\n' % self.ref_src)
        proc = Popen(['bwa', 'index', self.ref_src], stdout=PIPE, stderr=PIPE)
        _, err = proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError('`bwa index` exited with code %d' % proc.returncode)
        stderr.write('|  %s\n' % '\n|  '.join(err.strip().split('\n')))

    def ids(self):
        return {
            'reads': self.identify(self.read_src),
            'refs': self.identify(self.ref_src),
        }

    def mappings(self):
        self.index()
        ids = self.ids()
        indicator = self.indicator()
        indicator.start()
        for rec in SeqIO.parse(self.read_src, 'fasta'):
            indicator.progress()
            bwa_mem_opts = [
                '-A', '1',
                '-B', '2',
                '-O', '3',
                '-E', '1',
            ]
            lastz_opts = [
                '--ambiguous=n',
                '--format=softsam',
                '--gap=3,1',
                '--match=1,2',
                '--chain', # NOTE
            ]
            read_len = len(rec.seq)
            with NamedTemporaryFile(suffix='.fasta') as tmp:
                tmp.write('>%s\n%s\n' % (rec.id, str(rec.seq)))
                tmp.flush() # flush to disk, a subprocess will be accessing it

                cmd = ['/home/amir/lastz-distrib/bin/lastz'] + lastz_opts + [self.ref_src, tmp.name]
                #cmd = ['bwa', 'mem'] + bwa_mem_opts + [self.ref_src, tmp.name]
                proc = Popen(cmd, stdout=PIPE, stderr=PIPE)
                out, err = proc.communicate()
                #print '\n', out
                if proc.returncode != 0:
                    raise RuntimeError('`bwa mem` exited with code %d' % proc.returncode)

            # FIXME why does this get silently ignored?
            # cf https://www.python.org/dev/peps/pep-0479/
            try:
                aln = (l for l in out.strip().split('\n') if l[0] != '@').next()
            except StopIteration:
                continue
            # cf. SAM spec: https://samtools.github.io/hts-specs/SAMv1.pdf
            read, flag, ref, start, quality, cigar = aln.strip().split()[:6]
            start, flag = int(start), int(flag)
            if flag == 4:
                # Segment unmapped
                continue
            assert(flag in [0,16])
            direction = 1 if flag == 0 else -1
            len_on_ref = 0
            for match in finditer('\d+(S|M|D|I)', cigar):
                op = match.group()
                length = int(op[:-1])
                # soft clipping
                # FIXME double check and simplify the if soup:
                if op[-1] == 'S':
                    if match.start() == 0:
                        if direction == 1:
                            start -= length
                        else:
                            len_on_ref += length
                    elif match.end() == read_len:
                        if direction == 1:
                            len_on_ref += length
                        else:
                            start -= length
                elif op[-1] in 'MI':
                    len_on_ref += length
                else:
                    assert(op[-1] == 'D')

            end = start + direction * len_on_ref
            yield ids['reads'][read], Mapping(
                ref=ids['refs'][ref],
                strand='+' if direction == 1 else '-',
                ref_from=int(start),
                ref_to=end,
            )

        indicator.finish()

class WordReadMapper(ReadMapper):
    def save_mapping(self, mapping):
        with sqlite3.connect(self.seqdb.db) as conn:
            q = """
                INSERT INTO mappings (%s) VALUES (?, ?, ?, ?)
            """ % ', '.join(mapping._fields)
            conn.cursor().execute(q, tuple(mapping))

    def populate(self, **kw):
        raise NotImplementedError

    def map_read(self, readid):
        raise NotImplementedError

    def map_all(self):
        reads = self.seqdb.seqids(seq_type=seq.READ)
        num_reads = len(reads)
        indicator = ProgressIndicator(
            'Mapping %d reads' % num_reads, num_reads, percentage=False
        )
        indicator.start()
        for seqid in reads:
            yield self.map_read(seqid)
            indicator.progress()
        indicator.finish()

