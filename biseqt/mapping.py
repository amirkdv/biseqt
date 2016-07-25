# -*- coding: utf-8 -*-
import pysam
from itertools import combinations

from biseqt.util import ProgressIndicator
from biseqt.pw import Aligner, BANDED_MODE, B_OVERLAP

def _all_records(db):
    """Internal helper that loads all records from the database."""
    db.log('Loading all sequence records from database')
    records = list(db.find())
    return records


# TODO report the reference too, for now it's only one
def mappings_from_sam(db, sampath):
    """
    Yields:
        tuple:
            A tuple of sequence names (as found in original source) and the
            ``pysam.calignedsegment.AlignedSegment`` mapping it to the
            reference.
    """
    records_by_name = {rec.attrs['name']: rec for rec in _all_records(db)}
    samfile = pysam.AlignmentFile(sampath)
    db.log('Loading SAM mappings from %s.' % sampath)
    for mapping in samfile.fetch():
        name = mapping.query_name.rsplit('/', 1)[0]
        if name not in records_by_name:
            continue
        yield records_by_name[name], mapping


def overlaps_from_sam_mappings(db, sampath, min_overlap=-1):
    # dict from Record to mapping
    mappings = {rec.id: mapping
                for rec, mapping in mappings_from_sam(db, sampath)}
    db.log('Finding overlaps from SAM mappings.')
    seqids = sorted(mappings.keys())
    for id0, id1 in combinations(seqids, 2):
        map0, map1 = mappings[id0], mappings[id1]
        if map0.is_reverse != map1.is_reverse:
            continue
        # TODO ignoring query_alignment_start and query_alignment_end
        overlap_len = min(map0.reference_end, map1.reference_end) - \
                      max(map0.reference_start, map1.reference_start)
        if overlap_len <= 0 or overlap_len < min_overlap:
            continue
        yield id0, id1


# TODO what about rc?
def overlaps_from_seed_index(seed_index, seqids=None, targets=None,
    min_diagonal_score=20, **aligner_kw):
    records_by_id = {rec.id: rec for rec in _all_records(seed_index.db)}
    aligner_kw.update(alnmode=BANDED_MODE, alntype=B_OVERLAP)

    def _load_seq(_id):
        return seed_index.db.load_from_record(records_by_id[_id])

    if seqids is None:
        seqids = records_by_id.keys()
    elif isinstance(seqids, int):
        seqids = [seqids]
    assert isinstance(seqids, list)
    if targets is None:
        assert len(seqids) > 1
        all_to_all = True
        targets = seqids
    else:
        all_to_all = False
        assert not set(seqids).intersection(set(targets))

    indic = ProgressIndicator(num_total=len(seqids))
    indic.start()

    for seqid in seqids:
        indic.progress()
        seq = _load_seq(seqid)
        max_score = float('-inf')
        bands = {}
        for target in targets:
            if seqid == target:
                continue
            diag_range, score = seed_index.highest_scoring_band(seqid, target)
            if diag_range is None or score < min_diagonal_score:
                continue
            if not all_to_all and score < max_score:
                continue
            bands[target] = diag_range
        for target in bands:
            target_seq = _load_seq(target)
            with Aligner(seq, target_seq, diag_range=bands[target],
                         **aligner_kw) as aligner:
                score = aligner.solve()
                if score is None:
                    continue
                alignment = aligner.traceback()
                if alignment is not None:
                    yield seqid, target, alignment
    indic.finish()
