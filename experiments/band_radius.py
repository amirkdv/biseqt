#!/usr/bin/env python
import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt
from bisect import bisect_left
from scipy.special import erf
import re

from util import log, color_code, plot_with_sd, with_dumpfile, savefig
from util import sample_opseq, estimate_gap_probs_in_opseq
from util import load_fasta

from biseqt.blot import band_radius


def time_in_band(K, g, r):
    A = r / (np.sqrt(2 * g * K))
    return erf(A) \
        + A * (2 / np.sqrt(np.pi) * np.exp(-A ** 2) - 4 * A * (1 - erf(A)))


def sample_edit_sequences(K, g, n_samples, bio=False):
    if bio:
        pws_path = 'data/actn2/actn2-7vet-pws.fa'
        opseq = ''
        with open(pws_path) as f:
            for seq, name, _ in load_fasta(f):
                opseq += seq
        opseq = re.sub('D{10,}', 'D', opseq)
        opseq = re.sub('I{10,}', 'I', opseq)
        return sample_opseq(opseq, K, n_samples, g)
    else:
        return (''.join(np.random.choice(list('MID'),
                        p=[1-g, g / 2, g / 2], size=K))
                for _ in range(n_samples))


@with_dumpfile
def sim_time_in_band(K, gs, rs, n_samples, **kw):
    sim_data = {
        'in_band': {'simulated': np.zeros((len(gs), len(rs), n_samples)),
                    'biological': np.zeros((len(gs), len(rs), n_samples))},
        'gs': gs,
        'rs': rs,
        'K': K,
    }
    for g_idx, g in enumerate(gs):
        d0 = K
        for key, bio in zip(['simulated', 'biological'], [False, True]):
            log('sampling homologies for K = %d, g = %.2f (%s data)' %
                (K, g, key))
            samples = sample_edit_sequences(K, g, n_samples, bio=bio)
            for sample_idx, opseq in enumerate(samples):
                time_at_d_ = np.zeros(2*K)
                i, j = 0, 0
                for op in opseq:
                    d = i - j
                    time_at_d_[d + d0] += 1
                    if op in 'DM':
                        i += 1
                    if op in 'IM':
                        j += 1
                cum_time_at_d_ = np.cumsum(time_at_d_)
                for r_idx, r in enumerate(rs):
                    in_band = cum_time_at_d_[d0 + r] - cum_time_at_d_[d0 - r]
                    prop = in_band / K
                    assert prop <= 1
                    sim_data['in_band'][key][g_idx][r_idx][sample_idx] = prop
    return sim_data


def plot_time_in_band(sim_data, cutoff_epsilon, path=None):
    assert path
    gs = sim_data['gs']
    rs = sim_data['rs']
    K = sim_data['K']

    fig = plt.figure(figsize=(10, 6))

    grids = gridspec.GridSpec(2, 2, width_ratios=[1.5, 1])

    ax_mod = fig.add_subplot(grids[:, 0])  # model
    ax_sim = fig.add_subplot(grids[0, 1])  # simulation
    ax_bio = fig.add_subplot(grids[1, 1])  # biological data

    colors = color_code(gs, cmap='PuOr')

    for g_idx, (color, g) in enumerate(zip(colors, gs)):
        r_lim = rs[bisect_left(
            sim_data['in_band']['biological'][g_idx].mean(axis=1),
            1 - .5 * cutoff_epsilon
        )]
        r_cutoff = band_radius(K, g, 1 - cutoff_epsilon)

        vs = [erf(r / (2 * np.sqrt(g * K))) for r in rs]
        us = [time_in_band(K, g, r) for r in rs]

        kw = {'color': color, 'lw': 1.5, 'alpha': .8}
        ax_mod.plot(rs, vs, label='$g = %.2f$' % g, **kw)  # simplified model
        ax_mod.plot(rs, us, ls='--', **kw)                 # full correct model
        ax_mod.axvline(r_cutoff, color=color, lw=5, alpha=.3)

        ax_mod.grid(True)
        ax_mod.set_xlabel('diagonal band radius')
        ax_mod.set_ylabel('proportion of time in band')
        ax_mod.legend(loc='lower right', fontsize=12)
        ax_mod.set_xlim(0, r_lim)
        ax_mod.set_ylim(0, 1.2)

        for key, ax in zip(['simulated', 'biological'], [ax_sim, ax_bio]):
            res = sim_data['in_band'][key][g_idx, :, :]

            r_eff = rs[bisect_left(vs, 1 - cutoff_epsilon)]
            in_band_eff = sim_data['in_band'][key][g_idx, r_eff, :].mean()
            label = 'out of band: \%%%.2f' % (100 * (1 - in_band_eff))

            plot_with_sd(ax, rs, res, axis=1, y_max=1, color=color, lw=1.5,
                         label=label)
            ax.axvline(r_cutoff, color=color, lw=5, alpha=.3)
            ax.set_xlim(0, r_lim)
            ax.set_ylim(0, 1.2)
            ax.grid(True)
            ax.set_title('%s data' % key, fontsize=10)
            ax.set_xlabel('diagonal band radius', fontsize=8)
            ax.set_ylabel('proportion of time in band', fontsize=8)
            ax.legend(loc='lower right', fontsize=8)

    # n_samples = sim_data['in_band']['simulated'].shape[2]
    # fig.suptitle('$K = %d$, no. samples = $%d$' % (K, n_samples))
    # fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    savefig(fig, path, comment='time in band simulation')


def exp_time_in_band():
    """Plots proportion of alignment time (i.e. number of edit operations)
    spent within diagonal band for:

    * Simulated edit operations,
    * Biological edit operations obtained from real aligned genomes,
    * Full probabilistic model (expected amount of time spent in band),
    * Approximated probabilistic model (probability of last step in band).

    .. figure::
        https://www.dropbox.com/s/93cg9t4oh0h2guh/band_radius.png?raw=1
       :target:
        https://www.dropbox.com/s/93cg9t4oh0h2guh/band_radius.png?raw=1
       :alt: lightbox

       Predicted proportion of time spent in diagonal band (left) in
       alignments of total length K=500 according to full (dashed)
       and approximate (solid) probabilistic model and for two values of
       gap probabilities. Horizontal lines indicate the calculated band
       radius for %.01 accuracy. Empericial observations (n=500 samples)
       for proportion of time spent in band based on simulated edit sequences
       (top right) and real biological data (bottom right) are shown, and in
       each case the empirical percentage of time spent outside of band is
       reported.  Biological alignments are chosen from projected pairwise
       alignments of *Actinin 2* for 7 vertebrates obtained from UCSC genome
       browser. Each sample biological edit sequence is chosen such that it has
       the desired length and gap probability; gaps of length more than 10
       nucleotides are removed to avoid nonstationary gap probabilities.
    """
    K = 500
    gs = [.05, .15]
    n_samples = 500
    rs = range(0, 400)
    cutoff_epsilon = 1e-4  # for vertical lines showing calculated cutoff

    dumpfile = 'band_radius.txt'
    plot_path = 'band_radius.png'
    sim_data = sim_time_in_band(K, gs, rs, n_samples, dumpfile=dumpfile)
    plot_time_in_band(sim_data, cutoff_epsilon, path=plot_path)


def exp_gap_probs():
    opseqs_path = 'data/acta2/acta2-7vet-opseqs.fa'

    with open(opseqs_path) as f:
        opseq = f.read().strip()

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(1, 1, 1)

    K = 2000
    sample_len = 2 * K
    gap_hat = 0
    while gap_hat < .05:
        start = np.random.randint(0, len(opseq) - sample_len)
        sample = opseq[start:start + sample_len]
        gap_hat = 1. * (sample.count('I') + sample.count('D')) / sample_len

    g = gap_hat
    sim_sample = ''.join(np.random.choice(list('MID'),
                                          p=[1 - g, g / 2, g / 2],
                                          size=sample_len))

    radius = 200
    gaps = estimate_gap_probs_in_opseq(sample, radius)
    ax.plot(range(len(gaps)), gaps, c='g', lw=1, alpha=.9,
            label='biological')

    sim_gaps = estimate_gap_probs_in_opseq(sim_sample, radius)
    ax.plot(range(len(sim_gaps)), sim_gaps, c='k', lw=1,
            alpha=.9, label='simulated')

    ax.axhline(y=gap_hat, c='k', alpha=.2, lw=3)
    ax.grid(True)
    ax.legend(fontsize=11)
    ax.set_xlim(radius, sample_len - radius)
    ax.set_xlabel('position in edit sequence')
    ax.set_ylabel('estimated gap probability ($r=%d$)' % radius)
    ax.set_ylim(-.1, 1)
    savefig(fig, 'gap probabilities.png')


if __name__ == '__main__':
    exp_time_in_band()
    exp_gap_probs()
