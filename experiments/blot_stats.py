import numpy as np
import pickle
from matplotlib import pyplot as plt

import logging
from biseqt.blot import HomologyFinder
from biseqt.blot import band_radius
from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess

ALPHABET = Alphabet('ACGT')
LEISH_PATH = '../dense-band-alignment/biseqt-HEAD/leishmania/reference.fa'
load_real = False
if load_real:
    with open(LEISH_PATH) as f:
        LEISH_SEQ = ''.join([x.strip().upper() for x in f.readlines() if x[0] != '>'])

    LEISH_SEQ = [{'A': 0, 'C': 1, 'G': 2, 'T': 3}[s] for s in LEISH_SEQ if s in 'ACGT']


# FIXME remove duplicates bw here and complexity-analysis
def pickle_dump(path, obj, comment=None):
    with open(path, 'w') as f:
        pickle.dump(obj, f)
    print 'dumped %s to %s.' % (comment if comment else 'object', path)


def pickle_load(path):
    with open(path, 'rb') as f:
        return pickle.load(f)


def rand_seq_pair(N, related=False, real=False, subst=None, gap=None):
    if related:
        M = MutationProcess(ALPHABET, subst_probs=subst, ge_prob=gap, go_prob=gap)
    if not real:
        S = rand_seq(ALPHABET, N)
        if related:
            T, _ = M.mutate(S)
        else:
            T = rand_seq(ALPHABET, N)
    else:
        if related:
            assert N < len(LEISH_SEQ)
            nS = np.random.randint(0, len(LEISH_SEQ) - N)
            S = LEISH_SEQ[nS:nS+N]
            T, _ = M.mutate(S)
        else:
            assert N * 2 < len(LEISH_SEQ)
            nS = np.random.randint(0, len(LEISH_SEQ) - 2 * N)
            nT = np.random.randint(nS, len(LEISH_SEQ) - N)
            S = LEISH_SEQ[nS:nS+N]
            T = LEISH_SEQ[nT:nT+N]
    return S, T


# Add unrelated sequences to the end of the given sequences
def pad_seq_pair(S, T, pad_length):
    assert pad_length >= 0
    S_, T_ = rand_seq(ALPHABET, pad_length), rand_seq(ALPHABET, pad_length)
    return S + S_, T + T_

def adjust_pw_plot(ax, N0, N1):
    ax.set_ylim(-N0 * .1, N0 * 1.1)
    ax.set_xlim(-N1 * .1, N1 * 1.1)
    ax.set_xlabel('Pos. in T', fontsize=12)
    ax.set_ylabel('Pos. in S', fontsize=12)
    ax.set_aspect('equal')
    line_kw = {'c': 'k', 'alpha': .1, 'lw':'5'}
    ax.axvline(x=-5,  **line_kw)
    ax.axvline(x=N1+5, **line_kw)
    ax.axhline(y=-5,  **line_kw)
    ax.axhline(y=N0+5, **line_kw)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.tick_params(labelsize=10)


# FIXME absorb as much as possible into SeedIndex
def plot_seeds(ax, seeds, N0, N1, c='k'):
    idx_S, idx_T = [], []
    for seed in seeds:
        idx_S.append(seed[0])
        idx_T.append(seed[1])

    # x and y are flipped when going from matrix notation to plotting.
    ax.scatter(idx_T, idx_S, facecolor=c, edgecolor=None, lw=0, s=2, alpha=.3)
    ax.grid(True)

def plot_similar_segment(fig, ax, segment, color):
    (d_min, d_max), (a_min, a_max) = segment
    seg_ds = [d_min, d_min, d_max, d_max]
    seg_as = [a_min, a_max, a_max, a_min]

    seg_xs = [0, 0, 0, 0]
    seg_ys = [0, 0, 0, 0]
    for i in range(4):
        d, a = seg_ds[i], seg_as[i]
        x, y = (a + max(d, 0), a - min(d, 0))
        # x and y is flipped between matrix notation and plotting
        seg_xs[i], seg_ys[i] = y, x

    im = ax.fill(seg_xs, seg_ys, c=color, lw=5, alpha=.3)


# =======================
# Performance of Statistics in separating +/- cases
# Two scores are considered: z score wrt H0 and wrt H1
# =======================
# K is segment length of interest
def plot_stats_performance(fig, K, ns, n_samples, real=False, wordlen=None,
                           gap=None, subst=None, num_sds=1, banded=False,
                           replot=False, dumpfile=None):
    if replot:
        results = pickle_load(dumpfile)
    else:
        M = MutationProcess(ALPHABET, subst_probs=subst, ge_prob=gap, go_prob=gap)
        shape = (len(ns), n_samples)
        HF_kw = {'gap_prob': gap, 'subst_prob': subst,
                 'sensitivity': .9, 'alphabet': ALPHABET, 'wordlen': wordlen,
                 'path': ':memory:', 'log_level': logging.WARNING}
        results = {
            'H0': {'pos': np.zeros(shape), 'neg': np.zeros(shape)},
            'H1': {'pos': np.zeros(shape), 'neg': np.zeros(shape)},
            'banded': banded,
            'HF_kw': HF_kw,
        }
        for n_idx, n in enumerate(ns):
            print 'n = %d (%d/%d)' % (n, n_idx+1, len(ns))
            for i in range(n_samples):
                S, T = rand_seq_pair(K, real=real, related=False, mutation_process=M)
                S, T = pad_seq_pair(S, T, n-K)
                HF = HomologyFinder(S, T, **HF_kw)
                if banded:
                    radius = HF.band_radius(K)
                    A = 2 * n * radius
                    num_seeds = HF.seed_count(d_center=0, d_radius=radius)
                else:
                    A = n ** 2
                    num_seeds = HF.seed_count()
                s0, s1 = HF.score_num_seeds(num_seeds=num_seeds, area=A, seglen=K)
                results['H0']['neg'][n_idx][i] = s0
                results['H1']['neg'][n_idx][i] = s1

            for i in range(n_samples):
                S, T = rand_seq_pair(K, real=real, related=True, mutation_process=M)
                S, T = pad_seq_pair(S, T, n-K)
                HF = HomologyFinder(S, T, **HF_kw)
                if banded:
                    radius = HF.band_radius(K)
                    A = n * radius
                    num_seeds = HF.seed_count(d_center=0, d_radius=radius)
                else:
                    A = n ** 2
                    num_seeds = HF.seed_count()
                s0, s1 = HF.score_num_seeds(num_seeds=num_seeds, area=A, seglen=K)
                results['H0']['pos'][n_idx][i] = s0
                results['H1']['pos'][n_idx][i] = s1

        if dumpfile:
            pickle_dump(dumpfile, results, comment='stats performance')

    ax_H0 = fig.add_subplot(1, 2, 1)
    ax_H1 = fig.add_subplot(1, 2, 2)

    ax_H1.set_xlabel('Sequence lengths (kb)')
    ax_H0.set_xlabel('Sequence lengths (kb)')

    ax_H1.set_ylabel('Score')
    ax_H0.set_ylabel('Score')

    kw = {'marker': 'o', 'markersize': 5, 'lw': 2, 'alpha': .8}
    _ns = [n/1000. for n in ns]

    means = results['H1']['pos'].mean(axis=1)
    sds = np.sqrt(results['H1']['pos'].var(axis=1))
    ax_H1.plot(_ns, means, c='g', **kw)
    ax_H1.errorbar(_ns, means, yerr=num_sds * sds, c='g', markersize=5, marker='o')

    means = results['H1']['neg'].mean(axis=1)
    sds = np.sqrt(results['H1']['neg'].var(axis=1))
    ax_H1.plot(_ns, means, c='r', **kw)
    ax_H1.errorbar(_ns, means, yerr=num_sds * sds, c='r', markersize=5, marker='o')


    means = results['H0']['pos'].mean(axis=1)
    sds = np.sqrt(results['H0']['pos'].var(axis=1))
    ax_H0.plot(_ns, means, c='g', **kw)
    ax_H0.errorbar(_ns, means, yerr=num_sds * sds, c='g', markersize=5, marker='o')

    means = results['H0']['neg'].mean(axis=1)
    sds = np.sqrt(results['H0']['neg'].var(axis=1))
    ax_H0.plot(_ns, means, c='r', **kw)
    ax_H0.errorbar(_ns, means, yerr=num_sds * sds, c='r', markersize=5, marker='o')

    for ax in [ax_H1, ax_H0]:
        ax.grid(True)

    ax_H1.set_title('H1 %sscore' % ('banded ' if banded else ''))
    ax_H0.set_title('H0 %sscore' % ('banded ' if banded else ''))

def exp_stats_performance():
    fig = plt.figure(figsize=(14, 6))
    K = 500 # similar segment length
    ns = [i*K for i in range(1, 11)] # sequence lengths
    n_samples = 100 # number samples for each n

    wordlen = 8
    gap = .1
    subst = .1
    dumpfile = 'stats_performance_banded K=%d.txt' % K
    plot_stats_performance(fig, K, ns, n_samples, real=False, banded=True, subst=subst, gap=gap, wordlen=wordlen, dumpfile=dumpfile, replot=False)
    fig.suptitle('word len. = %d, segment len. = %d, # samples = %d' % (wordlen, K, n_samples), fontweight='bold')
    fig.tight_layout()
    fig.savefig('stats performance - banded K=%d.png' % K, dpi=300)

def exp_local_alignment():
    gap = .15
    subst = .15
    K = 500
    wordlen = 8
    # NOTE I can drive sensitivity to 0 and get decent results
    HF_kw = {'gap_prob': gap, 'subst_prob': subst,
             'sensitivity': .9, 'alphabet': ALPHABET, 'wordlen': wordlen,
             'path': 'example.db', 'log_level': logging.INFO}

    M = MutationProcess(ALPHABET, subst_probs=subst, ge_prob=gap, go_prob=gap)
    M_hgt = MutationProcess(ALPHABET, subst_probs=subst/3, ge_prob=gap/3, go_prob=gap/3)

    homs = [rand_seq(ALPHABET, i) for i in [200, 500, 1500, 2500]]
    junk = lambda: rand_seq(ALPHABET, np.random.randint(1000, 2000))

    S = junk() + homs[0] + junk() + homs[1] + junk() + homs[3] + \
        junk() + homs[2] + junk() + homs[2] + junk() + homs[0] + junk()
    homs = [M.mutate(homs[i])[0] if i != 2 else M_hgt.mutate(homs[i])[0]
            for i in range(len(homs))]
    T = junk() + homs[3] + junk() + homs[2] + junk() + homs[0] + \
        junk() + homs[3] + junk() + homs[1] + junk() + homs[2] + junk()


    import matplotlib.gridspec as gridspec
    fig = plt.figure()
    gs = gridspec.GridSpec(1, 2, width_ratios=[30,1])
    ax = plt.subplot(gs[0])
    ax_colorbar = plt.subplot(gs[1])

    import matplotlib
    cmap = plt.cm.get_cmap('jet')

    max_score = 1

    print 'finding homologies'
    HF = HomologyFinder(S, T, **HF_kw)
    scores = []
    for segment, score in HF.similar_segments(K):
        (d_min, d_max), (a_min, a_max) = segment
        print segment, score
        color = cmap(score/max_score)[:3]
        plot_similar_segment(fig, ax, segment, color=color)
        scores.append(score)

    norm = matplotlib.colors.Normalize(vmin=0, vmax=max_score)
    plot_seeds(ax, HF.seeds(), len(S), len(T))
    adjust_pw_plot(ax, len(S), len(T))

    cb1 = matplotlib.colorbar.ColorbarBase(ax_colorbar, cmap=cmap, norm=norm, orientation='vertical')

    fig.tight_layout()
    fig.savefig('local-alignment.png', dpi=300) # NOTE LOCAL ALIGNMENT WORKS !!! (Wed, 17 Jan 2018)

def exp_conserved_sequences():
    gap = .15
    subst = .15
    K = 500
    wordlen = 8
    HF_kw = {'gap_prob': gap, 'subst_prob': subst,
             'sensitivity': .9, 'alphabet': ALPHABET, 'wordlen': wordlen,
             'path': 'example.db', 'log_level': logging.INFO}

    M_std = MutationProcess(ALPHABET, subst_probs=.2, ge_prob=.15, go_prob=.15)
    M_con = MutationProcess(ALPHABET, subst_probs=.05, ge_prob=.05, go_prob=.05)

    pieces = [rand_seq(ALPHABET, 10000) for _ in range(3)]
    conserved = rand_seq(ALPHABET, 1000)

    S = pieces[0] + conserved + pieces[1] + pieces[2]
    pieces_hom = [M_std.mutate(p)[0] for p in pieces]
    conserved_hom = M_con.mutate(conserved)[0]
    T = pieces_hom[0] + pieces_hom[1] + conserved_hom + pieces_hom[2]

    import matplotlib.gridspec as gridspec
    fig = plt.figure()
    gs = gridspec.GridSpec(1, 2, width_ratios=[30,1])
    ax = plt.subplot(gs[0])
    ax_colorbar = plt.subplot(gs[1])

    import matplotlib
    cmap = plt.cm.get_cmap('jet')

    print 'finding homologies'
    HF = HomologyFinder(S, T, **HF_kw)
    scores = []
    for segment, score, match_prob in HF.similar_segments(K):
        (d_min, d_max), (a_min, a_max) = segment
        print segment, score, match_prob
        color = cmap(match_prob)[:3]
        plot_similar_segment(fig, ax, segment, color=color)
        scores.append(score)

    norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
    plot_seeds(ax, HF.seeds(), len(S), len(T))
    adjust_pw_plot(ax, len(S), len(T))

    cb1 = matplotlib.colorbar.ColorbarBase(ax_colorbar, cmap=cmap, norm=norm, orientation='vertical')

    fig.tight_layout()
    fig.savefig('conserved-sequences.png', dpi=300) # NOTE LOCAL ALIGNMENT WORKS !!! (Wed, 17 Jan 2018)


if __name__ == '__main__':
    #exp_stats_performance()
    #exp_local_alignment()
    exp_conserved_sequences()
