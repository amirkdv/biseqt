# -*- coding: utf-8 -*-
import pysam
from itertools import combinations

from biseqt.util import ProgressIndicator
from biseqt.pw import Aligner, BANDED_MODE, B_OVERLAP
from biseqt.database import DB, Record
from biseqt.kmers import KmerIndex
from biseqt.seeds import SeedIndex

# FIXME docs and tests


class Read(object):
    def __init__(self, seed_index, record, rc_record):
        assert isinstance(seed_index, SeedIndex)
        self.seed_index = seed_index
        assert isinstance(record, Record) and isinstance(rc_record, Record)
        assert rc_record.attrs.get('rc_of', None) == record.content_id
        self.record, self.rc_record = record, rc_record
        self.ids = [self.record.id, self.rc_record.id]
        self.seq = seed_index.db.load_from_record(self.record)
        self.rc_seq = seed_index.db.load_from_record(self.rc_record)

    def map(self, targets, min_band_score=20, **aligner_kw):
        if not isinstance(targets, list):
            targets = [targets]
        default = None, None, None
        assert all(isinstance(target, Record) for target in targets)
        targetids = [target.id for target in targets]
        (ours, match), band, score = self.seed_index.highest_scoring_band(
            self.ids, targetids, min_band_score=min_band_score)

        if ours is None:
            return default
        assert ours in self.ids and match in targetids
        target = [target for target in targets if target.id == match][0]
        rc = ours == self.rc_record.id
        seq = self.seq if not rc else self.rc_seq
        target_seq = self.seed_index.db.load_from_record(target)
        aligner_kw.update(alnmode=BANDED_MODE, alntype=B_OVERLAP,
                          diag_range=band)
        with Aligner(seq, target_seq, **aligner_kw) as aligner:
            score = aligner.solve()
            if score is None:
                return default
            alignment = aligner.traceback()
            if alignment is not None:
                ours = self.record if not rc else self.rc_record
                # FIXME maybe return True False for rc only, users already
                # have this object
                return ours, target, alignment


class ReadMapper(object):
    def __init__(self, alphabet, wordlen, db_path):
        self.db = DB(db_path, alphabet)
        self.kmer_index = KmerIndex(self.db, wordlen)
        self.seed_index = SeedIndex(self.kmer_index)
        self.bands_indexed = False

    def log(self, *args, **kwargs):
        self.db.log(*args, **kwargs)

    def initialize(self, reads_fa, refs_fa=None, num_reads=-1):
        self.db.initialize()
        with open(reads_fa) as f:
            self.db.load_fasta(f, num=num_reads, rc=True)
        if refs_fa is not None:
            with open(refs_fa) as f:
                self.db.load_fasta(f, rc=False)

    def index_bands(self, **kw):
        self.kmer_index.score_kmers()
        self.seed_index.score_diagonals(**kw)
        self.bands_indexed = True

    def load_reads(self):
        recs_by_content_id = {r.content_id: r for r in list(self.db.find())}
        reads = []
        for record in recs_by_content_id.values():
            if 'rc_of' in record.attrs:
                pair = (recs_by_content_id[record.attrs['rc_of']], record)
                reads.append(Read(self.seed_index, *pair))
        return sorted(reads, key=lambda read: read.record.id)

    def load_refs(self):
        recs_by_content_id = {r.content_id: r for r in list(self.db.find())}
        for record in recs_by_content_id.values():
            if 'rc_of' in record.attrs:
                recs_by_content_id.pop(record.attrs['rc_of'])
                recs_by_content_id.pop(record.content_id)
        return recs_by_content_id.values()

    def map_all_to_all(self, min_band_score, **aligner_kw):
        assert self.bands_indexed, 'Bands must be indexed first'
        self.log('Mapping all reads against each other')
        reads = self.load_reads()  # NOTE comes in sorted order of id
        indic = ProgressIndicator(num_total=len(reads))
        indic.start()
        for read in reads:
            indic.progress()
            # NOTE only compare to reads after us
            others = (r.record for r in reads if r.record.id > read.record.id)
            for other in others:
                rec, target_rec, aln = read.map(other,
                                                min_band_score=min_band_score,
                                                **aligner_kw)
                if rec is None:
                    continue
                yield rec, target_rec, aln
        indic.finish()

    def map_all_to_refs(self, min_band_score, **aligner_kw):
        # FIXME it would be nice to only calculate bands for read v. ref not
        # all pairwise of reads too
        assert self.bands_indexed, 'Bands must be indexed first'
        self.log('Mapping all reads against reference sequences')
        reads, refs = self.load_reads(), self.load_refs()
        indic = ProgressIndicator(num_total=len(reads))
        indic.start()
        for read in reads:
            indic.progress()
            rec, target_rec, aln = read.map(refs,
                                            min_band_score=min_band_score,
                                            **aligner_kw)
            if rec is not None:
                yield rec, target_rec, aln
        indic.finish()

    def mappings_from_sam(self, sampath):
        """Loads mappings from a SAM mapping file and translates sequence
        names to integer identifiers as stored by :class:`biseqt.database.DB`.

        Args:
            db (database.DB): The sequence database where ids are looked up.
            sampath (str): The path to SAM mappings file.

        Yields:
            tuple:
                A 3-tuple containing the read record, the reference name and
                the ``pysam.calignedsegment.AlignedSegment`` mapping it to the
                reference.
        """
        self.log('Loading SAM mappings from %s.' % sampath)
        reads_by_name = {r.record.attrs['name']: r for r in self.load_reads()}
        samfile = pysam.AlignmentFile(sampath)
        for mapping in samfile.fetch():
            qname, rname = mapping.query_name, mapping.reference_name
            # NOTE this is because BLASR does a weird thing with sequence names
            qname = qname.rsplit('/', 1)[0]
            if qname not in reads_by_name:
                continue
            yield reads_by_name[qname], rname, mapping

    def overlaps_from_sam_mappings(self, sampath, min_overlap=-1):
        """Finds all pairs of overlapping sequences based on their mappings to
        a reference.

        Args:
            sampath (str): The path to SAM mappings file.
            min_overlap (int): The minimum required length for overlaps to be
                reported; default is -1 in which case no overlap is excluded.

        Yields:
            tuple:
                A tuple of sequence integer ids (in increasing order) that are
                deemed as overlapping based on SAM mappings.
        """
        self.log('Finding overlaps from SAM mappings.')
        mappings = {read.record.id: (read, ref, mapping)
                    for read, ref, mapping in self.mappings_from_sam(sampath)}
        seqids = sorted(mappings.keys())
        for id0, id1 in combinations(seqids, 2):
            (r0, ref0, map0), (r1, ref1, map1) = mappings[id0], mappings[id1]
            if ref0 != ref1:
                continue
            # TODO ignoring query_alignment_start and query_alignment_end
            overlap_len = min(map0.reference_end, map1.reference_end) - \
                max(map0.reference_start, map1.reference_start)
            if overlap_len <= 0 or overlap_len < min_overlap:
                continue
            # FIXME the second thing we yield is not reported by our own
            # map_all_to_all.
            if map0.is_reverse == map1.is_reverse:
                yield r0.record, r1.record
                yield r0.rc_record, r1.rc_record
            else:
                yield r0.record, r1.rc_record
                yield r0.rc_record, r1.record
