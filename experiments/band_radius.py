#!/usr/bin/env python
import numpy as np
from matplotlib import pyplot as plt
from bisect import bisect_left
from scipy.special import erf, erfcinv
from util import color_code, with_dumpfile, savefig
plt.rc('text', usetex=True)


def time_in_band(K, g, r):
    A = r / (2 * np.sqrt(g * K))
    return erf(A) \
        + A * (2 / np.sqrt(np.pi) * np.exp(-A ** 2) - 4 * A * (1 - erf(A)))


def real_sensitivity(K, g, epsilon):
    r = erfcinv(epsilon) * 2 * np.sqrt(g * K)
    return time_in_band(K, g, r)


@with_dumpfile
def sim_time_in_band(K, gs, rs, n_samples, **kw):
    sim_data = {
        'in_band': np.zeros((len(gs), len(rs), n_samples)),
        'gs': gs,
        'rs': rs,
        'K': K,
    }
    for g_idx, g in enumerate(gs):
        d0 = K
        for sample_idx in range(n_samples):
            print '%d / %d' % (sample_idx + 1, n_samples)
            path = np.random.choice(['M', 'I', 'D'], p=[1-2*g, g, g], size=K)
            time_at_d_ = np.zeros(2*K)
            i, j = 0, 0
            for op in path:
                d = i - j
                time_at_d_[d + d0] += 1
                if op in 'DM':
                    i += 1
                if op in 'IM':
                    j += 1
            cum_time_at_d_ = np.cumsum(time_at_d_)
            for r_idx, r in enumerate(rs):
                in_band = cum_time_at_d_[d0 + r] - cum_time_at_d_[d0 - r]
                assert in_band <= K
                sim_data['in_band'][g_idx][r_idx][sample_idx] = in_band / K
    return sim_data


def plot_time_in_band(sim_data, cutoff_epsilon, path=None):
    assert path and dumpfile
    gs = sim_data['gs']
    rs = sim_data['rs']
    K = sim_data['K']
    n_samples = sim_data['in_band'].shape[2]

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1, 2, 1)
    ax_sim = fig.add_subplot(1, 2, 2)

    colors = color_code(gs)

    for g_idx, (color, g) in enumerate(zip(colors, gs)):
        vs = [erf(r / (2 * np.sqrt(g * K))) for r in rs]
        r_lim = rs[bisect_left(vs, 1 - .5 * cutoff_epsilon)]
        r_cutoff = rs[bisect_left(vs, 1 - cutoff_epsilon)]
        us = [time_in_band(K, g, r) for r in rs]
        kw = {'color': color, 'lw': 1.5, 'alpha': .8}
        ax.plot(rs, vs, label='$g = %.2f$' % g, **kw)  # simplified model
        ax.plot(rs, us, ls='--', **kw)                 # full correct model
        ax.axvline(r_cutoff, color=color, lw=5, alpha=.3)
        ax.grid(True)
        ax.set_xlabel('Band radius')
        ax.set_ylabel('proportion of time in band')
        ax.legend(loc='lower right', fontsize=12)
        ax.set_xlim(0, r_lim)
        ax.set_ylim(0, 1.2)

        res = sim_data['in_band'][g_idx, :, :]

        num_sds = 1
        means = res.mean(axis=1)
        sds = np.sqrt(res.var(axis=1))
        ys_l = means - num_sds * sds
        ys_h = np.minimum(1, means + num_sds * sds)
        ax_sim.plot(rs, means, lw=1.5, color=color, zorder=2)
        ax_sim.fill_between(rs, ys_l, ys_h, facecolor=color, edgecolor=color,
                            alpha=.2)
        ax_sim.set_xlabel('Band radius')
        ax_sim.set_ylabel('proportion of time in band')
        ax_sim.axvline(r_cutoff, color=color, lw=5, alpha=.3)
        ax_sim.set_xlim(0, r_lim)
        ax_sim.set_ylim(0, 1.2)
        ax_sim.grid(True)

    r_cutoff = rs[bisect_left(vs, 1 - cutoff_epsilon)]
    print 'effective mean time in band (sim.) = %f' % \
        (sim_data['in_band'][0, r_cutoff, :].mean())

    fig.suptitle('$K = %d$, \# samples = $%d$' % (K, n_samples))
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    savefig(fig, path, comment='time in band simulation')


def plot_real_vs_calculated_sensitivity(K, g, epsilons, path=None):
    assert path
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(1, 1, 1)
    # NOTE g and K don't effect the relationship!! use .1 and 1000 for example
    real_epsilons = [1. - real_sensitivity(K, g, epsilon)
                     for epsilon in epsilons]
    ax.plot(epsilons, real_epsilons, color='k', lw=1.5, alpha=.8)
    ax.grid(True)
    ax.set_xlim(0, 1e-3)
    ax.set_ylim(0, 1.1e-2)
    ax.set_xlabel('$\\hat{\\epsilon}$', fontsize=18)
    ax.set_ylabel('$\epsilon$', fontsize=18)
    savefig(fig, path, comment='approximate vs exact sensitivity')


def exp_time_in_band():
    K = 5000
    gs = [.05, .3]
    n_samples = 1000
    rs = range(0, 400)
    cutoff_epsilon = 1e-3  # for vertical lines showing calculated cutoff

    dumpfile = 'band_radius.txt'
    plot_path = 'band_radius.png'
    sim_data = sim_time_in_band(K, gs, rs, n_samples, dumpfile=dumpfile,
                                ignore_existing=False)
    plot_time_in_band(sim_data, cutoff_epsilon, path=plot_path)

    plot_path = 'approximated_sensitivity.png'
    epsilons = np.arange(2e-5, 1e-3, 2e-5)  # for comparing with "exact" model
    plot_real_vs_calculated_sensitivity(1000, .1, epsilons, path=plot_path)


if __name__ == '__main__':
    exp_time_in_band()
