#!/usr/bin/env python
from matplotlib import pyplot as plt
import sys
from time import time
import numpy as np
from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess
from biseqt.seeds import SeedIndex
from biseqt.blot import band_radius, H0_moments, H1_moments
from util import seq_pair, color_code
from util import plot_with_sd, with_dumpfile, log, savefig
from logging import WARN


@with_dumpfile
def sim_count_seeds(**kw):
    ns, n_samples = kw['ns'], kw['n_samples']
    gap, subst, wordlen = kw['gap'], kw['subst'], kw['wordlen']
    log('simulating # of seeds (%d samples), lengths = %s' %
        (n_samples, str(ns)))

    def _zero(): return np.zeros((len(ns), n_samples))

    A = Alphabet('ACGT')
    seed_index_kw = {
        'alphabet': A,
        'wordlen': wordlen,
        'path': kw.get('path', ':memory:'),
        'log_level': WARN,
    }
    sim_data = {
        'time': {'pos': _zero(), 'neg': _zero()},
        'n_seeds': {'pos': _zero(), 'neg': _zero()},
        'gap': gap,
        'match': (1 - gap) * (1 - subst),
        'ns': ns,
        'seed_index_kw': seed_index_kw,
    }
    M = MutationProcess(A, go_prob=gap, ge_prob=gap, subst_probs=subst)
    for n_idx, n in enumerate(ns):
        log('n = %d ' % n, newline=False)
        for idx in range(n_samples):
            sys.stderr.write('.')
            S_rel, T_rel = seq_pair(n, A, mutation_process=M)
            S_urel, T_urel = rand_seq(A, n), rand_seq(A, n)
            for key, (S, T) in zip(['pos', 'neg'],
                                   [(S_rel, T_rel), (S_urel, T_urel)]):
                t = time()
                seed_index = SeedIndex(S, T, **seed_index_kw)
                n_seeds = seed_index.seed_count()
                sim_data['time'][key][n_idx][idx] = time() - t
                sim_data['n_seeds'][key][n_idx][idx] = n_seeds
        sys.stderr.write('\n')
    return sim_data


def plot_count_seeds_moments(sim_data, K=None, suffix=''):
    ns, wordlen = sim_data['ns'], sim_data['seed_index_kw']['wordlen']
    match = sim_data['match']

    kw = {'marker': 'o', 'markersize': 3, 'lw': 1.5, 'alpha': .5}

    mus_H0, sds_H0 = [], []
    mus_H1, sds_H1 = [], []
    for n in ns:
        area = n ** 2
        mu_H0, sd_H0 = H0_moments(4, wordlen, area)
        mus_H0.append(mu_H0)
        sds_H0.append(sd_H0)

        mu_H1, sd_H1 = H1_moments(4, wordlen, area, n, match)
        mus_H1.append(mu_H1)
        sds_H1.append(sd_H1)

    fig = plt.figure(figsize=(14, 5))
    ax_t = fig.add_subplot(1, 3, 1)
    ax_mu = fig.add_subplot(1, 3, 2)
    ax_sd = fig.add_subplot(1, 3, 3)

    # time to find all seeds
    plot_with_sd(ax_t, ns, 1000 * sim_data['time']['neg'], axis=1, color='r',
                 label='unrelated', **kw)
    plot_with_sd(ax_t, ns, 1000 * sim_data['time']['pos'], axis=1, color='g',
                 label='related', **kw)

    # average no. of seeds
    pos = sim_data['n_seeds']['pos']
    neg = sim_data['n_seeds']['neg']
    ax_mu.plot(ns, neg.mean(axis=1), c='r', label='unrelated', **kw)
    ax_mu.plot(ns, mus_H0, color='r', ls='--', **kw)
    ax_mu.plot(ns, pos.mean(axis=1), c='g', label='related', **kw)
    ax_mu.plot(ns, mus_H1, color='g', ls='--', **kw)

    # std dev. of no of seeds
    ax_sd.plot(ns, np.sqrt(neg.var(axis=1)), c='r', label='unrelated', **kw)
    ax_sd.plot(ns, sds_H0, color='r', ls='--', **kw)
    ax_sd.plot(ns, np.sqrt(pos.var(axis=1)), c='g', label='related', **kw)
    ax_sd.plot(ns, sds_H1, color='g', ls='--', **kw)

    _ns = np.arange(min(ns), max(ns))
    ax_mu.plot(_ns, _ns, color='k', alpha=.9, lw=.5, ls='--')
    ax_t.plot(_ns, _ns, color='k', alpha=.9, lw=.5, ls='--')

    for ax in [ax_sd, ax_mu, ax_t]:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('sequence length', fontsize=10)
        ax.set_xticks(ns)
        ax.set_xticklabels(ns, fontsize=10)
        ax.legend(loc='best', fontsize=10)

    ax_sd.set_ylabel('standard deviation of no. of matching %d-mers' % wordlen)
    ax_mu.set_ylabel('average no. of matching %d-mers' % wordlen)
    ax_t.set_ylabel('time to find matching %d-mers (ms)' % wordlen)

    n_samples = pos.shape[1]
    fig.suptitle('no. samples = %d, match = %.2f' % (n_samples, match),
                 fontsize=8)
    savefig(fig, 'num_seeds_moments%s.png' % suffix)


def exp_count_seeds():
    ns = [200 * 2 ** i for i in range(8)]
    gap = .1
    subst = .1
    n_samples = 50
    wordlen = 8
    suffix = '[w=%d]' % wordlen
    dumpfile = 'num_seeds%s.txt' % suffix
    sim_data = sim_count_seeds(
        ns=ns, n_samples=n_samples, gap=gap, subst=subst,
        wordlen=wordlen, dumpfile=dumpfile, ignore_existing=False)
    plot_count_seeds_moments(sim_data, suffix=suffix)


@with_dumpfile
def sim_count_seeds_segment(**kw):
    Ks, g_radii, n_samples = kw['Ks'], kw['g_radii'], kw['n_samples']
    gap, subst, wordlen = kw['gap'], kw['subst'], kw['wordlen']

    def _zero(): return np.zeros((len(Ks), len(g_radii), n_samples))

    A = Alphabet('ACGT')
    seed_index_kw = {
        'alphabet': A,
        'wordlen': wordlen,
        'path': kw.get('path', ':memory:'),
        'log_level': WARN,
    }
    sim_data = {
        'n_seeds': {'pos': _zero(), 'neg': _zero()},
        'p_hat': {'pos': _zero(), 'neg': _zero()},
        'gap': gap,
        'match': (1 - gap) * (1 - subst),
        'Ks': Ks,
        'g_radii': g_radii,
        'seed_index_kw': seed_index_kw,
    }
    M = MutationProcess(A, go_prob=gap, ge_prob=gap, subst_probs=subst)
    for K_idx, K in enumerate(Ks):
        log('K = %d' % K, newline=False)
        for idx in range(n_samples):
            sys.stderr.write('.')
            S_rel, T_rel = seq_pair(K, A, mutation_process=M)
            S_urel, T_urel = rand_seq(A, K), rand_seq(A, K)
            for key, (S, T) in zip(['pos', 'neg'],
                                   [(S_rel, T_rel), (S_urel, T_urel)]):
                seed_index = SeedIndex(S, T, **seed_index_kw)
                for g_idx, g_max in enumerate(g_radii):
                    radius = band_radius(K, g_max, 1 - 1e-4)
                    d_band = (-radius, radius)
                    n_seeds = seed_index.seed_count(d_band=d_band)
                    sim_data['n_seeds'][key][K_idx, g_idx, idx] = n_seeds

                    area = 2 * radius * K
                    word_p_null = (1./len(A)) ** wordlen
                    word_p = (n_seeds - area * word_p_null) / K
                    try:
                        p_hat = np.exp(np.log(word_p) / wordlen)
                    except Warning:
                        # presumably this happened because word_p was too small
                        p_hat = 0
                    p_hat = min(p_hat, 1)
                    sim_data['p_hat'][key][K_idx, g_idx, idx] = p_hat
        sys.stderr.write('\n')
    return sim_data


def plot_count_seeds_segment(sim_data, suffix=''):
    Ks, g_radii = sim_data['Ks'], sim_data['g_radii']
    match = sim_data['match']

    kw = {'marker': 'o', 'markersize': 3, 'lw': 1, 'alpha': .6}

    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(1, 1, 1)

    colors = color_code(g_radii)
    for g_idx, (g_max, color) in enumerate(zip(g_radii, colors)):
        pos = sim_data['p_hat']['pos'][:, g_idx, :]
        neg = sim_data['p_hat']['neg'][:, g_idx, :]
        plot_with_sd(ax, Ks, neg, axis=1, color=color, ls='--', **kw)
        plot_with_sd(ax, Ks, pos, axis=1, color=color,
                     label='$g_{\max} = %.2f$' % g_max, **kw)

    ax.set_xlabel('homology length', fontsize=10)
    ax.set_xticks(Ks)
    ax.set_xticklabels(Ks, fontsize=6)
    ax.legend(loc='best', fontsize=6)
    ax.set_ylim(-.2, 1.1)
    ax.set_ylabel('estimated match probability')

    n_samples = pos.shape[1]
    fig.suptitle('no. samples = %d, match = %.2f' % (n_samples, match),
                 fontsize=8)
    savefig(fig, 'num_seeds_segment%s.png' % suffix)


def exp_count_seeds_segment():
    Ks = [100 * i for i in range(1, 21)]
    g_radii = [.05, .1, .2, .4]
    gap = .1
    subst = .15
    n_samples = 50
    wordlen = 6
    suffix = '[w=%d]' % wordlen
    dumpfile = 'num_seeds_segment%s.txt' % suffix
    sim_data = sim_count_seeds_segment(
        Ks=Ks, g_radii=g_radii, n_samples=n_samples,
        gap=gap, subst=subst, wordlen=wordlen,
        dumpfile=dumpfile, ignore_existing=False)
    plot_count_seeds_segment(sim_data, suffix=suffix)


if __name__ == '__main__':
    exp_count_seeds()
    exp_count_seeds_segment()
