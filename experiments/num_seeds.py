#!/usr/bin/env python
from matplotlib import pyplot as plt
from time import time
import numpy as np
from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess
from biseqt.seeds import SeedIndex
from biseqt.blot import band_radius, H0_moments, H1_moments
from util import seq_pair, seqs_bio, bio_opseqs, apply_opseq, load_fasta
from util import plot_with_sd, with_dumpfile, log, savefig
from logging import WARN


@with_dumpfile
def sim_count_seeds(ns, n_samples, gap=None, subst=None, banded=False, **kw):
    log('simulating # of seeds: %d samples, ns = %s' % (n_samples, str(ns)))

    def _zero(): return np.zeros((len(ns), n_samples))

    sim_data = {
        'time': {'pos': _zero(), 'neg': _zero()},
        'n_seeds': {'pos': _zero(), 'neg': _zero()},
        'gap': gap,
        'match': (1 - gap) * (1 - subst),
        'ns': ns,
        'wordlen': kw['wordlen'],
        'banded': banded,
    }
    A = Alphabet('ACGT')
    M = MutationProcess(A, go_prob=gap, ge_prob=gap, subst_probs=subst)
    for n_idx, n in enumerate(ns):
        # FIXME perform banded and non-banded experiments in one-go!! here and
        # in sim_count_seeds_bio; need to fix plot_count_seeds, and
        # plot_count_seeds_moments to produce two plots each (banded and not)
        if banded:
            radius = band_radius(n, gap, 1 - 1e-4)
            d_band = (-radius, radius)
        else:
            d_band = None
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
                n_seeds = seed_index.seed_count(d_band=d_band)
                sim_data['time'][key][n_idx][idx] = 1000 * (time() - t)
                sim_data['n_seeds'][key][n_idx][idx] = n_seeds
    return sim_data


@with_dumpfile
def sim_count_seeds_bio(ns, n_samples, gap=None, match=None, banded=False,
                        bio_source_seq=None, bio_source_opseq=None, **kw):
    log('simulating # of seeds: %d samples, ns = %s' % (n_samples, str(ns)))

    def _zero(): return np.zeros((len(ns), n_samples))

    sim_data = {
        'time': {'pos': _zero(), 'neg': _zero()},
        'n_seeds': {'pos': _zero(), 'neg': _zero()},
        'gap': gap,
        'match': match,
        'ns': ns,
        'wordlen': kw['wordlen'],
        'banded': banded,
    }
    A = Alphabet('ACGT')
    for n_idx, n in enumerate(ns):
        log('n = %d' % n)
        # NOTE just use gap for band radius purposes, use match for finding
        # opseqs
        opseqs = list(bio_opseqs(bio_source_opseq, n, n_samples, gap=None,
                                 match=match))
        if banded:
            radius = band_radius(n, gap, 1 - 1e-4)
            d_band = (-radius, radius)
        else:
            d_band = None
        for idx in range(n_samples):
            S_urel, T_urel = seqs_bio([n, n], bio_source_seq)
            S_rel = seqs_bio([n], bio_source_seq)[0]
            T_rel = apply_opseq(S_rel, opseqs[idx])
            for key, (S, T) in zip(['pos', 'neg'],
                                   [(S_rel, T_rel), (S_urel, T_urel)]):
                t = time()
                seed_index = SeedIndex(S, T, wordlen=kw['wordlen'],
                                       path=kw.get('path', ':memory:'),
                                       alphabet=A, log_level=WARN)
                n_seeds = seed_index.seed_count(d_band=d_band)
                sim_data['time'][key][n_idx][idx] = 1000 * (time() - t)
                sim_data['n_seeds'][key][n_idx][idx] = n_seeds
    return sim_data


# automatically banded
@with_dumpfile
def sim_count_seeds_fixed_K(ns, K, n_samples, gap=None, subst=None, **kw):
    log('simulating # of seeds: %d samples, ns = %s' % (n_samples, str(ns)))

    def _zero(): return np.zeros((len(ns), n_samples))

    sim_data = {
        'n_seeds': {'pos': _zero(), 'neg': _zero()},
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
                sim_data['n_seeds'][key][n_idx][idx] = n_seeds
    return sim_data


@with_dumpfile
def sim_count_seeds_fixed_K_bio(ns, K, n_samples, gap=None, match=None,
                                bio_source_seq=None, bio_source_opseq=None,
                                **kw):
    log('simulating # of seeds: %d samples, ns = %s' % (n_samples, str(ns)))

    def _zero(): return np.zeros((len(ns), n_samples))

    sim_data = {
        'n_seeds': {'pos': _zero(), 'neg': _zero()},
        'gap': gap,
        'match': match,
        'ns': ns,
        'K': K,
        'wordlen': kw['wordlen'],
    }
    assert min(ns) >= K
    A = Alphabet('ACGT')
    # NOTE just use gap for band radius purposes, use match for finding opseqs
    opseqs = list(bio_opseqs(bio_source_opseq, K, n_samples, gap=None,
                             match=match))
    for n_idx, n in enumerate(ns):
        log('n = %d' % n)
        radius = band_radius(K, gap, 1 - 1e-4)
        d_band = (-radius, radius)
        for idx in range(n_samples):
            S_urel, T_urel = seqs_bio([n, n], bio_source_seq)
            S_rel = seqs_bio([K], bio_source_seq)[0]
            T_rel = apply_opseq(S_rel, opseqs[idx])
            junk = seqs_bio([n - K, n - K], bio_source_seq)
            S_rel += junk[0]
            T_rel += junk[1]

            for key, (S, T) in zip(['pos', 'neg'],
                                   [(S_rel, T_rel), (S_urel, T_urel)]):
                seed_index = SeedIndex(S, T, wordlen=kw['wordlen'],
                                       path=kw.get('path', ':memory:'),
                                       alphabet=A, log_level=WARN)
                n_seeds = seed_index.seed_count(d_band=d_band)
                sim_data['n_seeds'][key][n_idx][idx] = n_seeds
    return sim_data


def plot_count_seeds(sim_data, suffix=''):
    log('plotting # of seed stats')
    fig = plt.figure(figsize=(11, 5))
    ax_n = fig.add_subplot(1, 2, 1)
    ax_t = fig.add_subplot(1, 2, 2)

    ns = sim_data['ns']
    wordlen = sim_data['wordlen']
    kw = {'n_sds': 1, 'marker': 'o', 'markersize': 5, 'axis': 1}
    plot_with_sd(ax_n, ns, sim_data['n_seeds']['pos'], color='g',
                 label='homologous', **kw)
    plot_with_sd(ax_n, ns, sim_data['n_seeds']['neg'], color='r',
                 label='nonhomologous', **kw)

    plot_with_sd(ax_t, ns, sim_data['time']['pos'], color='g', **kw)
    plot_with_sd(ax_t, ns, sim_data['time']['neg'], color='r', **kw)

    ax_n.set_ylabel('No. of matching kmers', fontsize=10)
    ax_t.set_ylabel('Time to find matching kmers (ms)', fontsize=10)
    for ax in [ax_t, ax_n]:
        ax.set_xlabel('(non)homologous sequence length', fontsize=10)
        ax.grid(True)
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)

    ax_n.legend(loc='upper left', fontsize=10)
    n_samples = sim_data['n_seeds']['pos'].shape[1]
    fig.suptitle('wordlen = %d, no. samples = %d' % (wordlen, n_samples))
    savefig(fig, 'num_seeds%s.png' % suffix)


def plot_count_seeds_fixed_K(sim_data, suffix=''):
    log('plotting # of seed stats')
    K, ns = sim_data['K'], sim_data['ns']
    wordlen = sim_data['wordlen']

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(1, 1, 1)
    kw = {'n_sds': 1, 'marker': 'o', 'markersize': 3, 'alpha': .7, 'axis': 1}
    plot_with_sd(ax, ns, sim_data['n_seeds']['pos'], color='g',
                 label='homologous', **kw)
    plot_with_sd(ax, ns, sim_data['n_seeds']['neg'], color='r',
                 label='non-homologous', **kw)
    ax.set_xlabel('full sequence length')
    ax.set_ylabel('No. of exactly matching kmers')

    ax.legend(loc='upper left', fontsize=10)
    n_samples = sim_data['n_seeds']['pos'].shape[1]
    fig.suptitle('K = %d, wordlen = %d, no. samples = %d' %
                 (K, wordlen, n_samples))
    savefig(fig, 'num_seeds%s.png' % suffix)


# only concerned with diagonal band as with sim_count_seeds_fixed_K
# if K is None it means in each experiment the whole sequences where
# homologous or not (i.e. as in sim_count_seeds{,_bio}), otherwise it's a fixed
# number (i.e. as in sim_count_seeds_fixed_K{,_bio}).
def plot_count_seeds_moments(sim_data, K=None, suffix=''):
    ns, gap, match = sim_data['ns'], sim_data['gap'], sim_data['match']
    wordlen, banded = sim_data['wordlen'], sim_data['banded']

    fig = plt.figure(figsize=(12, 5))
    ax_mu = fig.add_subplot(1, 2, 1)
    ax_sd = fig.add_subplot(1, 2, 2)

    kw = {'marker': 'o', 'markersize': 3, 'lw': 1.5, 'alpha': .7}

    mus_H0, sds_H0 = [], []
    mus_H1, sds_H1 = [], []
    for n in ns:
        K_ = n if K is None else K
        assert 0 < K_ <= n
        if banded:
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

    ax_mu.plot(ns, sim_data['n_seeds']['neg'].mean(axis=1), color='r', **kw)
    ax_mu.plot(ns, mus_H0, color='r', ls='--', **kw)
    ax_mu.plot(ns, sim_data['n_seeds']['pos'].mean(axis=1), color='g', **kw)
    ax_mu.plot(ns, mus_H1, color='g', ls='--', **kw)

    ax_sd.plot(ns, np.sqrt(sim_data['n_seeds']['neg'].var(axis=1)), color='r',
               **kw)
    ax_sd.plot(ns, sds_H0, color='r', ls='--', **kw)
    ax_sd.plot(ns, np.sqrt(sim_data['n_seeds']['pos'].var(axis=1)), color='g',
               **kw)
    ax_sd.plot(ns, sds_H1, color='g', ls='--', **kw)

    for ax in [ax_sd, ax_mu]:
        ax.grid(True)
        ax.set_xlabel('(non)-homologous sequence length')
    ax_sd.set_ylabel('standard deviation of no. of matching kmers')
    ax_mu.set_ylabel('expected no. of matching kmers')

    n_samples = sim_data['n_seeds']['pos'].shape[1]
    fig.suptitle('wordlen = %d, no. samples = %d' % (wordlen, n_samples))
    savefig(fig, 'num_seeds_moments%s.png' % suffix)


def exp_count_seeds():
    ns = [200 * 2 ** i for i in range(8)]
    gap = .1
    subst = .1
    n_samples = 50
    wordlen = 8
    dumpfile = 'num_seeds.txt'
    sim_data = sim_count_seeds(ns, n_samples, gap=gap, subst=subst,
                               wordlen=wordlen, banded=False,
                               dumpfile=dumpfile,
                               ignore_existing=False)
    plot_count_seeds(sim_data)
    plot_count_seeds_moments(sim_data)

    suffix = '[banded]'
    dumpfile = 'num_seeds%s.txt' % suffix
    sim_data = sim_count_seeds(ns, n_samples, gap=gap, subst=subst,
                               wordlen=wordlen, banded=True,
                               dumpfile=dumpfile,
                               ignore_existing=False)
    plot_count_seeds(sim_data, suffix=suffix)
    plot_count_seeds_moments(sim_data, suffix=suffix)

    # ===============
    # Biological Data
    # ===============
    match = (1 - gap) * (1 - subst)
    A = Alphabet('ACGT')
    with open('data/leishmenia/reference.fa') as f:
        bio_source_seq = ''.join(s for s, _, _ in load_fasta(f))
        bio_source_seq = bio_source_seq.upper()
        bio_source_seq = ''.join(x for x in bio_source_seq if x in 'ACGT')
    with open('data/leishmenia/blasr_opseqs.fa') as f:
        bio_source_opseq = f.read().strip()

    suffix = '[bio]'
    dumpfile = 'num_seeds%s.txt' % suffix
    sim_data_bio = sim_count_seeds_bio(ns, n_samples, gap=gap, match=match,
                                       wordlen=wordlen, banded=False,
                                       bio_source_seq=A.parse(bio_source_seq),
                                       bio_source_opseq=bio_source_opseq,
                                       dumpfile=dumpfile,
                                       ignore_existing=False)
    plot_count_seeds(sim_data_bio, suffix=suffix)
    plot_count_seeds_moments(sim_data_bio, suffix=suffix)

    suffix = '[banded][bio]'
    dumpfile = 'num_seeds%s.txt' % suffix
    sim_data_bio = sim_count_seeds_bio(ns, n_samples, gap=gap, match=match,
                                       wordlen=wordlen, banded=True,
                                       bio_source_seq=A.parse(bio_source_seq),
                                       bio_source_opseq=bio_source_opseq,
                                       dumpfile=dumpfile,
                                       ignore_existing=False)
    plot_count_seeds(sim_data_bio, suffix=suffix)
    plot_count_seeds_moments(sim_data_bio, suffix=suffix)


# FIXME estimating moments of fixed K doesn't work (we underestimate everything
# in all cases: mean and sd in pos and neg) for biological data. This is
# important for band discovery. I think the issue is:
# 1. kmer occurances in biological data is not uniform and hence pw_H1 is more
#    than what we think (multiplying it by 2.2 with w=8 seems to fix the neg
#    case, both mean and sd).
# 2. finding unrelated biological data is problematic! using 7vet data is a bad
#    idea because we must be sampling homologous chunks, for some reason
#    leishmenia chromosome 1 doesn't work either. But when I keep only the
#    homologous part biologic and the non-homologous part synthetic everything
#    is fine; of course the neg data is not biological anymore, but pos works
#    fine for mean and sd. The latter fact is basically the same as plotting
#    moments with K increasing because there we don't have the issue of
#    specifying non homologous biological data.
def exp_count_seeds_fixed_K():
    K = 200
    ns = [K * 2 ** i for i in range(8)]
    gap = .1
    subst = .1
    n_samples = 50
    wordlen = 8
    suffix = '[fixed_K]'
    dumpfile = 'num_seeds%s.txt' % suffix
    sim_data = sim_count_seeds_fixed_K(
        ns, K, n_samples, gap=gap, subst=subst,
        wordlen=wordlen,
        dumpfile=dumpfile,
        ignore_existing=False)
    plot_count_seeds_fixed_K(sim_data, suffix=suffix)
    plot_count_seeds_moments(sim_data, K=K, suffix=suffix)

    # ===============
    # Biological Data
    # ===============
    A = Alphabet('ACGT')
    with open('data/acta2/acta2-7vet.fa') as f:
        bio_source_seq = ''.join(s for s, _, _ in load_fasta(f))
        bio_source_seq = bio_source_seq.upper()
        bio_source_seq = ''.join(x for x in bio_source_seq if x in 'ACGT')
    with open('data/acta2/acta2_opseqs.fa') as f:
        bio_source_opseq = ''.join(s for s, _, _ in load_fasta(f))

    suffix = '[fixed_K][bio]'
    dumpfile = 'num_seeds%s.txt' % suffix
    match = (1 - gap) * (1 - subst)
    sim_data_bio = sim_count_seeds_fixed_K_bio(
        ns, K, n_samples, gap=gap, match=match,
        wordlen=wordlen,
        bio_source_seq=A.parse(bio_source_seq),
        bio_source_opseq=bio_source_opseq,
        dumpfile=dumpfile,
        ignore_existing=False)
    plot_count_seeds_fixed_K(sim_data_bio, suffix=suffix)
    plot_count_seeds_moments(sim_data_bio, K=K, suffix=suffix)


if __name__ == '__main__':
    exp_count_seeds()
    exp_count_seeds_fixed_K()
