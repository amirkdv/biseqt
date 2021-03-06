#!/usr/bin/env python
import numpy as np
from matplotlib import gridspec
from matplotlib import pyplot as plt
from bisect import bisect_left
from scipy.special import erf, erfcinv
from scipy.ndimage.filters import gaussian_filter1d
import re

from util import log, color_code, plot_with_sd, with_dumpfile, savefig
from util import sample_opseq, estimate_gap_probs_in_opseq
from util import load_fasta, opseq_path

from biseqt.seeds import SeedIndex
from biseqt.blot import band_radius


def time_in_band(K, gap, radius):
    alpha = radius / (np.sqrt(2 * gap * K))

    return erf(alpha) + alpha * (2 / np.sqrt(np.pi) * np.exp(-alpha ** 2) -
                                 2 * alpha * (1 - erf(alpha)))


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


def plot_time_in_band(sim_data, cutoff_epsilon):
    gs = sim_data['gs']
    rs = sim_data['rs']
    K = sim_data['K']
    n_samples = sim_data['in_band']['simulated'].shape[2]

    fig = plt.figure(figsize=(10, 8))

    grids = gridspec.GridSpec(3, 2, height_ratios=[1, 1, 1])

    ax_mod = fig.add_subplot(grids[:2, 0])  # model
    ax_sim = fig.add_subplot(grids[2, 0])  # simulated data
    ax_bio = fig.add_subplot(grids[2, 1])  # biological data
    ax_rw1 = fig.add_subplot(grids[0, 1])
    ax_rw2 = fig.add_subplot(grids[1, 1])

    colors = color_code(gs, cmap='PuOr')

    for g_idx, (color, g, ax_rw) in enumerate(zip(colors, gs,
                                                  [ax_rw1, ax_rw2])):
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

        ax_mod.set_title('model predictions', fontsize=10)
        ax_mod.grid(True)
        ax_mod.set_xlabel('diagonal band radius')
        ax_mod.set_ylabel('proportion of time in band')
        ax_mod.legend(loc='lower right', fontsize=12)
        ax_mod.set_xlim(0, r_lim)
        ax_mod.set_ylim(0, 1.2)

        kw = {'lw': 1, 'alpha': .2, 'color': color}
        for opseq in sample_edit_sequences(K, g, n_samples, bio=False):
            xs, ys = opseq_path(opseq, x0=0, y0=0)
            ds = [SeedIndex.to_diagonal_coordinates(x, y)[0]
                  for x, y in zip(xs, ys)]
            ax_rw.plot(range(K + 1), ds, **kw)

        for eps, lw, alpha, label in zip([1e-1, 1e-6], [2, 3], [.7, .3],
                                         ['$\epsilon=10^{-1}$',
                                          '$\epsilon=10^{-6}$']):

            def radius_(k): return erfcinv(eps) * np.sqrt(2 * g) * np.sqrt(k)

            plot_kw = {'c': 'k', 'alpha': alpha, 'lw': lw}
            radii = np.array([radius_(k) for k in range(K)])
            ax_rw.plot(range(K), radii, label=label, **plot_kw)
            ax_rw.plot(range(K), -radii, **plot_kw)

        ax_rw.legend(loc='upper left', fontsize=8)
        ax_rw.set_ylim(-50, 50)
        ax_rw.set_xlabel('number of alignment steps')
        ax_rw.set_ylabel('diagonal position $d$')

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

    savefig(fig, 'band_radius.png', comment='time in band')


def exp_time_in_band():
    """Shows theoretical and empirical results for proportion of alignment time
    (i.e. number of edit operations) spent within diagonal bands of varying
    radius.

    **Supported Claims**

    * The random walk model for capturing the probability distribution of
      diagonal positions is accurate.
    * Calculating the band radius such that the *endpoint* of an alignment is
      within diagonal band with probability at least :math:`1 - \epsilon` is a
      good conservative estimate for the more reasonable requirement that the
      *average proportion* of time spent in band is at least :math:`1 -
      \epsilon`.
    * Calculated band radii are accurate and reliable when applied to
      simulated as well as biological edit sequences.

    .. figure::
        https://www.dropbox.com/s/93cg9t4oh0h2guh/band_radius.png?raw=1
       :target:
        https://www.dropbox.com/s/93cg9t4oh0h2guh/band_radius.png?raw=1
       :alt: lightbox

       Predicted proportion of time spent in diagonal band (*top left*) in
       alignments of total length K=500 according to full (dashed) and
       approximate (solid) probabilistic model and for two values of gap
       probabilities. Vertical lines indicate the calculated band radius for
       %.01 accuracy.
       Simulated alignments (n=500 samples) are drawn in diagonal coordinates
       as a function of alignment time for gap probability 0.05 (*right top*)
       and 0.15 (*right middle*) with bounds obtained from
       :func:`biseqt.blot.band_radius` for two sensitivity values.
       Proportion of time spent in band as a function of band radius is shown
       for simulated sequences (*bottom left*) and real biological data
       (*bottom right*), and in each case the empirical percentage of time
       spent *outside* of band is reported; shaded regions indicate one
       standard deviation.  Biological alignments are chosen from projected
       pairwise alignments of *alpha-Actinin 2 (ACTN2)* for 7 vertebrates
       obtained from UCSC genome browser. Each sample biological edit sequence
       is chosen such that it has the desired length and gap probability; gaps
       of length more than 10 nucleotides are removed to avoid nonstationary
       gap probabilities.
    """
    K = 500
    gs = [.05, .15]
    n_samples = 500
    rs = range(0, 400)
    cutoff_epsilon = 1e-4  # for vertical lines showing calculated cutoff

    dumpfile = 'band_radius.txt'
    sim_data = sim_time_in_band(K, gs, rs, n_samples, dumpfile=dumpfile)
    plot_time_in_band(sim_data, cutoff_epsilon)


def exp_gap_probs():
    """Shows gap probability variability in simulated and biological edit
    sequences.

    **Supported Claims**

    * Gap probabilities are nonstationary in biological data.
    * Removing long indel stretches mitigates this nonstationarity for the most
      part and allows us to treat gap probability in biological alignments
      as stationary.

    .. figure::
        https://www.dropbox.com/s/vplg0x22ykc47am/gap_probabilities.png?raw=1
       :target:
        https://www.dropbox.com/s/vplg0x22ykc47am/gap_probabilities.png?raw=1
       :alt: lightbox

       Gap probability variability along alignments (length K=1200) in
       simulated (black) and biological edit sequences (green), before (left)
       and after (right) removing long indel stretches (deletion or insertion
       stretches of length at least 10) from biological samples; n=50 samples
       in both cases. Local gap probabilites are calculated in a window of
       radius 100 nt.  Simulated edit sequences (black) are generated using
       stationary gap probability (i.e constant along alignment) of 0.15.
       Biological edit sequences are sampled from projected pairwise alignments
       of *Actinin 2* for 7 vertebrates obtained from UCSC genome browser such
       that the overall gap probability (i.e over K=1200 operations) is 0.15.
    """
    pws_path = 'data/actn2/actn2-7vet-pws.fa'
    full_opseq = ''
    with open(pws_path) as f:
        for seq, name, _ in load_fasta(f):
            full_opseq += seq

    fig = plt.figure(figsize=(10, 4))
    ax = fig.add_subplot(1, 2, 1)
    ax_stationary = fig.add_subplot(1, 2, 2)

    gap = .15
    n_samples = 50
    K = 1200
    radius = 100

    for mode, ax in zip(['stationary', 'original'], [ax_stationary, ax]):
        opseq = full_opseq
        if mode == 'stationary':
            opseq = re.sub('D{10,}', 'D', opseq)
            opseq = re.sub('I{10,}', 'I', opseq)

        real_samples = sample_opseq(opseq, K, n_samples, gap)
        real_labeled = False
        for real_sample in real_samples:
            label = 'biological' if not real_labeled else None
            real_gaps = gaussian_filter1d(
                estimate_gap_probs_in_opseq(real_sample, radius), radius/5)
            ax.plot(range(K), real_gaps, c='g', lw=1, alpha=.2, label=label)
            real_labeled = True

        sim_labeled = False
        for _ in range(n_samples):
            sim_sample = ''.join(
                np.random.choice(
                    list('MID'), p=[1 - gap, gap / 2, gap / 2], size=K
                )
            )
            label = 'simulated' if not sim_labeled else None
            sim_gaps = gaussian_filter1d(
                estimate_gap_probs_in_opseq(sim_sample, radius), radius/5)
            ax.plot(range(K), sim_gaps, c='k', lw=1, alpha=.4, label=label)
            sim_labeled = True

        ax.axhline(y=gap, c='k', alpha=.2, lw=3)
        ax.grid(True)
        ax.legend(loc='upper right', fontsize=8)
        ax.set_xlim(radius, K - radius)
        ax.set_xlabel('position in edit sequence')
        ax.set_ylabel('gap probability')
        ax.set_ylim(-.05, .6)

    savefig(fig, 'gap_probabilities.png')


if __name__ == '__main__':
    exp_time_in_band()
    exp_gap_probs()
