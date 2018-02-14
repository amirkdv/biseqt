#!/usr/bin/env python
import pysam
import re
from util import log
import numpy as np


def opseq_from_cigar_string(cigarstring):
    cigar_ops = [x for x in re.split('(I|M|D)', cigarstring) if x]
    assert len(cigar_ops) % 2 == 0, 'incomplete cigar string'
    ops = []
    for i in range(len(cigar_ops) / 2):
        ops += [cigar_ops[2 * i + 1]] * int(cigar_ops[2 * i])
    assert all(op in 'IMD' for op in ops)
    return ''.join(ops)


def estimate_gap_probs(opseq, radius):
    assert isinstance(radius, int)
    assert len(opseq) > 2 * radius
    gaps = np.zeros(len(opseq))
    n_M = opseq[:2 * radius].count('M')
    for i in range(radius, len(opseq) - radius):
        if opseq[i - radius] == 'M':
            n_M -= 1
        if opseq[i + radius] == 'M':
            n_M += 1
        gaps[i] = 1 - n_M / (2. * radius)
    return gaps


def sam_to_opseqs(sam_path, opseq_path):
    samfile = pysam.AlignmentFile(sam_path)
    log('loading SAM records from %s.' % sam_path)
    with open(opseq_path, 'w') as f:
        for rec in samfile.fetch():
            opseq = opseq_from_cigar_string(rec.cigarstring)
            assert opseq.count('I') + opseq.count('M') == \
                rec.query_alignment_length
            f.write(opseq)

        f.write('\n')


def sample_opseq(opseq, K, gap, n_samples, resolution=1e-3):
    radius = K / 2
    gaps = estimate_gap_probs(opseq, radius)
    cands = [idx for idx in range(len(gaps))
             if abs(round(gaps[idx], 2) - gap) < resolution]
    assert len(cands) > n_samples, 'not enough samples found'
    samples = np.random.choice(cands, size=n_samples)
    for pos in samples:
        yield opseq[pos - radius:pos + radius]
