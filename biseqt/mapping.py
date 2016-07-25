# -*- coding: utf-8 -*-
from biseqt.pw import Aligner, BANDED_MODE, B_OVERLAP


class Mapper(object):
    def __init__(self, seed_index, min_diagonal_score=20):
        self.seed_index = seed_index
        self.db = seed_index.db
        self.min_diagonal_score = min_diagonal_score

    # FIXME what about rc?
    def map(self, seqids, targets=[], to_all=False, **aligner_kw):
        aligner_kw.update(alnmode=BANDED_MODE, alntype=B_OVERLAP)
        def load_seq_by_id(_id):
            record = self.db.find(sql_condition='id = %d' % _id).next()
            return self.db.load_from_record(record)

        mappings = {}

        if not targets:
            targets = seqids
        for seqid in seqids:
            mappings[seqid] = {}
            seq = load_seq_by_id(seqid)
            # FIXME fix hack in highest_scoring_band
            bands = {to: self.seed_index.highest_scoring_band(seqid, to)
                     for to in targets}
            bands = {to: bands[to] for to in bands
                     if bands[to] != (None, None) and
                         bands[to][1] >= self.min_diagonal_score}
            if not to_all:
                best_target = max(bands, lambda to: bands[to][1])
                bands = {best_target:bands[best_target]}
            for to in bands:
                diag_range = bands[to][0]
                target_seq = load_seq_by_id(to)
                with Aligner(seq, target_seq, diag_range=diag_range,
                             **aligner_kw) as aligner:
                    aligner.solve()
                    alignment = aligner.traceback()
                    if alignment is not None:
                        mappings[seqid][to] = alignment
        return mappings
