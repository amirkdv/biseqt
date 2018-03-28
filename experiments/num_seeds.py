#!/usr/bin/env python
from matplotlib import pyplot as plt
from time import time
import numpy as np
from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess
from biseqt.seeds import SeedIndex
from biseqt.blot import band_radius, H0_moments, H1_moments
from util import seq_pair
from util import plot_with_sd, with_dumpfile, log, savefig
from logging import WARN


# FIXME keyword arguments are split between kw and else, merge them!
@with_dumpfile
def sim_count_seeds(ns, n_samples, gap=None, subst=None, banded=False, **kw):
    log('simulating # of seeds: %d samples, ns = %s' % (n_samples, str(ns)))

    def _zero(): return np.zeros((len(ns), n_samples))

    sim_data = {
        'time': {'pos': _zero(), 'neg': _zero()},
        'n_seeds': {
            'full': {'pos': _zero(), 'neg': _zero()},
            'banded': {'pos': _zero(), 'neg': _zero()},
        },
        'gap': gap,
        'match': (1 - gap) * (1 - subst),
        'ns': ns,
        'wordlen': kw['wordlen'],
    }
    A = Alphabet('ACGT')
    M = MutationProcess(A, go_prob=gap, ge_prob=gap, subst_probs=subst)
    for n_idx, n in enumerate(ns):
        radius = band_radius(n, gap, 1 - 1e-4)
        d_band = (-radius, radius)
        log('n = %d' % n)
        for idx in range(n_samples):
            S_rel, T_rel = seq_pair(n, A, mutation_process=M)
            S_urel, T_urel = rand_seq(A, n), rand_seq(A, n)
            for key, (S, T) in zip(['pos', 'neg'],
                                   [(S_rel, T_rel), (S_urel, T_urel)]):
                t = time()
                seed_index = SeedIndex(S, T, wordlen=kw['wordlen'],
                                       path=kw.get('path', ':memory:'),
                                       alphabet=A, log_level=WARN)
                # full table
                n_seeds = seed_index.seed_count(d_band=None)
                sim_data['time'][key][n_idx][idx] = 1000 * (time() - t)
                sim_data['n_seeds']['full'][key][n_idx][idx] = n_seeds

                # banded
                n_seeds = seed_index.seed_count(d_band=d_band)
                sim_data['n_seeds']['banded'][key][n_idx][idx] = n_seeds
    return sim_data


@with_dumpfile
def sim_count_seeds_fixed_K(ns, K, n_samples, gap=None, subst=None, **kw):
    log('simulating # of seeds: %d samples, ns = %s' % (n_samples, str(ns)))

    def _zero(): return np.zeros((len(ns), n_samples))

    sim_data = {
        'n_seeds': {  # only consider banded
            'banded': {'pos': _zero(), 'neg': _zero()},
        },
        'gap': gap,
        'match': (1 - gap) * (1 - subst),
        'ns': ns,
        'K': K,
        'wordlen': kw['wordlen'],
    }
    assert min(ns) >= K
    A = Alphabet('ACGT')
    M = MutationProcess(A, go_prob=gap, ge_prob=gap, subst_probs=subst)
    for n_idx, n in enumerate(ns):
        radius = band_radius(K, gap, 1 - 1e-4)
        d_band = (-radius, radius)
        log('n = %d' % n)
        for idx in range(n_samples):
            S_rel, T_rel = seq_pair(K, A, mutation_process=M)
            S_rel += rand_seq(A, n - K)
            T_rel += rand_seq(A, n - K)
            S_urel, T_urel = rand_seq(A, n), rand_seq(A, n)
            for key, (S, T) in zip(['pos', 'neg'],
                                   [(S_rel, T_rel), (S_urel, T_urel)]):
                seed_index = SeedIndex(S, T, wordlen=kw['wordlen'],
                                       path=kw.get('path', ':memory:'),
                                       alphabet=A, log_level=WARN)
                n_seeds = seed_index.seed_count(d_band=d_band)
                sim_data['n_seeds']['banded'][key][n_idx][idx] = n_seeds
    return sim_data


def plot_time_seeds(sim_data, suffix=''):
    fig_t = plt.figure(figsize=(6, 5))
    ax_t = fig_t.add_subplot(1, 1, 1)

    ns = sim_data['ns']
    wordlen = sim_data['wordlen']
    kw = {'n_sds': 1, 'marker': 'o', 'markersize': 5, 'axis': 1}

    plot_with_sd(ax_t, ns, sim_data['time']['pos'], color='g', alpha=.6,
                 label='homologous', **kw)
    plot_with_sd(ax_t, ns, sim_data['time']['neg'], color='r', alpha=.6,
                 label='non-homologous', **kw)
    ax_t.set_ylabel('Time to find all matching kmers (ms)', fontsize=10)

    ax_t.set_xlabel('(non)homologous sequence length', fontsize=10)
    ax_t.grid(True)
    ax_t.tick_params(axis='x', labelsize=10)
    ax_t.tick_params(axis='y', labelsize=10)
    ax_t.legend(fontsize=10)
    ax_t.set_xticks(ns)
    ax_t.set_xticklabels(ns, rotation=90, fontsize=6)

    n_samples = sim_data['time']['pos'].shape[1]
    fig_t.suptitle('wordlen = %d, no. samples = %d' % (wordlen, n_samples))
    savefig(fig_t, 'time_seeds%s.png' % suffix)


# only concerned with diagonal band as with sim_count_seeds_fixed_K
# if K is None it means in each experiment the whole sequences where
# homologous or not (i.e. as in sim_count_seeds), otherwise it's a fixed
# number (i.e. as in sim_count_seeds_fixed_K).
def plot_count_seeds_moments(sim_data, K=None, suffix=''):
    ns, wordlen = sim_data['ns'], sim_data['wordlen']
    gap, match = sim_data['gap'], sim_data['match']

    kw = {'marker': 'o', 'markersize': 3, 'lw': 1.5, 'alpha': .7}

    for key in ['full', 'banded']:
        if key not in sim_data['n_seeds']:
            # sometimes (e.g. fixed K, we don't bother with full table)
            continue
        mus_H0, sds_H0 = [], []
        mus_H1, sds_H1 = [], []
        for n in ns:
            K_ = n if K is None else K
            assert 0 < K_ <= n
            if key == 'banded':
                radius = band_radius(K_, gap, 1 - 1e-4)
                area = 2 * radius * n
            else:
                area = n ** 2
            mu_H0, sd_H0 = H0_moments(4, wordlen, area)
            mus_H0.append(mu_H0)
            sds_H0.append(sd_H0)

            mu_H1, sd_H1 = H1_moments(4, wordlen, area, K_, match)
            mus_H1.append(mu_H1)
            sds_H1.append(sd_H1)

        fig = plt.figure(figsize=(12, 5))
        ax_mu = fig.add_subplot(1, 2, 1)
        ax_sd = fig.add_subplot(1, 2, 2)

        pos = sim_data['n_seeds'][key]['pos']
        neg = sim_data['n_seeds'][key]['neg']
        ax_mu.plot(ns, neg.mean(axis=1), color='r', **kw)
        ax_mu.plot(ns, mus_H0, color='r', ls='--', **kw)
        ax_mu.plot(ns, pos.mean(axis=1), color='g', **kw)
        ax_mu.plot(ns, mus_H1, color='g', ls='--', **kw)

        ax_sd.plot(ns, np.sqrt(neg.var(axis=1)), color='r', **kw)
        ax_sd.plot(ns, sds_H0, color='r', ls='--', **kw)
        ax_sd.plot(ns, np.sqrt(pos.var(axis=1)), color='g', **kw)
        ax_sd.plot(ns, sds_H1, color='g', ls='--', **kw)

        for ax in [ax_sd, ax_mu]:
            ax.set_xlabel('(non)-homologous sequence length')
            ax.set_ylim(-.1 * ax.get_ylim()[1], None)
            ax.set_xticks(ns)
            ax.set_xticklabels(ns, rotation=90, fontsize=6)

        ax_sd.set_ylabel('standard deviation of no. of matching kmers')
        ax_mu.set_ylabel('expected no. of matching kmers')

        n_max = min(ax_mu.get_ylim()[1], ax_mu.get_xlim()[1])
        _ns = np.arange(n_max)
        ax_mu.plot(_ns, _ns, color='k', alpha=.9, lw=.5, ls='--')

        n_samples = pos.shape[1]
        fig.suptitle('wordlen = %d, no. samples = %d' % (wordlen, n_samples))
        savefig(fig, 'num_seeds_moments[%s]%s.png' % (key, suffix))


def exp_count_seeds():
    ns = [200 * 2 ** i for i in range(8)]
    gap = .1
    subst = .1
    n_samples = 50
    wordlen = 8
    suffix = '[w=%d]' % wordlen
    dumpfile = 'num_seeds%s.txt' % suffix
    sim_data = sim_count_seeds(ns, n_samples, gap=gap, subst=subst,
                               wordlen=wordlen,
                               dumpfile=dumpfile,
                               ignore_existing=False)
    plot_time_seeds(sim_data, suffix=suffix)
    plot_count_seeds_moments(sim_data, suffix=suffix)


def exp_count_seeds_fixed_K():
    K = 200
    ns = [K * 2 ** i for i in range(8)]
    gap = .1
    subst = .1
    n_samples = 50
    wordlen = 8
    suffix = '[fixed_K][w=%d]' % wordlen
    dumpfile = 'num_seeds%s.txt' % suffix
    sim_data = sim_count_seeds_fixed_K(
        ns, K, n_samples, gap=gap, subst=subst,
        wordlen=wordlen,
        dumpfile=dumpfile,
        ignore_existing=False)
    plot_count_seeds_moments(sim_data, K=K, suffix=suffix)


if __name__ == '__main__':
    exp_count_seeds()
    exp_count_seeds_fixed_K()
