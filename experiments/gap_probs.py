#!/usr/bin/env python
from matplotlib import pyplot as plt
import numpy as np
from util import estimate_gap_probs_in_opseq, savefig


# TODO clean this up and make it an experiment in band_radius.py
if __name__ == '__main__':
    opseqs_path = 'data/leishmenia/blasr_opseqs.fa'

    n = 20000

    with open(opseqs_path) as f:
        blasr_opseq = f.read().strip()
    start = np.random.randint(0, len(blasr_opseq) - n)
    blasr_opseq = blasr_opseq[start:start + n]
    blasr_gap_hat = 1 - 1. * blasr_opseq.count('M') / len(blasr_opseq)

    g = blasr_gap_hat
    sim_opseq = np.random.choice(list('MID'), p=[1-g, g / 2, g / 2], size=n)
    sim_opseq = ''.join(sim_opseq)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    for K, alpha, ls, lw in zip([500, 2000], [.3, .9], ['-', '--'], [1, 2]):
        radius = K / 2
        blasr_gaps = estimate_gap_probs_in_opseq(blasr_opseq, radius)
        ax.plot(range(len(blasr_gaps)), blasr_gaps, ls=ls, c='g',
                label='bio. data (K = %d)' % K, alpha=alpha)

        sim_gaps = estimate_gap_probs_in_opseq(sim_opseq, radius)
        ax.plot(range(len(sim_gaps)), sim_gaps, ls=ls, c='k',
                label='sim. data (K = %d)' % K, alpha=alpha)

    ax.axhline(y=blasr_gap_hat, c='k', alpha=.2, lw=3)
    ax.grid(True)
    ax.legend(fontsize=11)
    ax.set_xlim(K, n - K)
    ax.set_ylim(0, .5)
    savefig(fig, 'gap probabilities.png')
