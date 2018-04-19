#!/usr/bin/env python
import numpy as np
import sys
import logging
from itertools import product
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from Bio import AlignIO
from time import time

from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess
from biseqt.blot import WordBlotMultiple, band_radius
from util import log, get_seqs_from_mse, with_dumpfile
from util import color_code, savefig, plot_roc, plot_with_sd


# seeds are assumed to be in standard coordinates
def plot_scored_seeds_3d(fig, ax, scored_seeds, threshold=.5):
    idx_S, idx_T1, idx_T2, cs, ss = [], [], [], [], []
    cmap = plt.cm.get_cmap('Greys')
    scores = [score for _, score in scored_seeds]
    max_score = max(scores)
    for (i, j, k), score in scored_seeds:
        idx_S.append(i)
        idx_T1.append(j)
        idx_T2.append(k)
        cs.append(cmap(score/max_score)[:3])
        ss.append(10 if score > threshold else 1)

    ax.scatter(idx_S, idx_T1, idx_T2, facecolor=cs, lw=0, s=ss, alpha=.3)
    ax.set_aspect('equal')
    ax.elev = 10
    ax.azim = 240
    norm = matplotlib.colors.Normalize(vmin=0, vmax=max_score)
    m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    m.set_array(scores)
    fig.colorbar(m, shrink=.7)


# segment is a list of 3 tuples (min, max of each coordinate)
# NOTE util.plot_similar_segment for pw takes segments in diagonal coordinates!
def plot_similar_segment_3d(ax, segment, color='k', **kw):
    assert len(segment) == 3
    assert all(len(range) == 2 for range in segment)
    ax.plot(segment[0], segment[1], segment[2], c='m', alpha=.4, lw=2)


# ========================================================
# Mostly for debugging purposes, 3d plot for simulated MSE
# ========================================================
def exp_three_syntehtic_sequences():
    subst, gap = .05, .09
    wordlen = 6
    K = 2000
    A = Alphabet('ACGT')
    M = MutationProcess(A, subst_probs=subst, ge_prob=gap, go_prob=gap)
    S = rand_seq(A, K)
    T1, _ = M.mutate(S)
    T2, _ = M.mutate(S)

    def junk(): return rand_seq(A, np.random.randint(K))

    S = junk() + S + junk()
    T1 = junk() + T1 + junk()
    T2 = junk() + T2 + junk()

    WB_kw = {'g_max': .2, 'sensitivity': .9, 'alphabet': A, 'wordlen': wordlen,
             'path': ':memory:'}
    WB = WordBlotMultiple(S, T1, T2, **WB_kw)

    p_min = (1-gap) * (1-subst)
    p_min = .9 * p_min
    scored_seeds = [(WB.to_ij_coordinates(*rec['seed']), rec['p'])
                    for rec in WB.score_seeds(K)]

    fig = plt.figure()
    ax = fig.gca(projection=Axes3D.name)
    plot_scored_seeds_3d(fig, ax, scored_seeds, threshold=p_min)
    for res in WB.similar_segments(K_min=K, p_min=p_min, score=True):
        std_ranges = WB.to_ij_coordinates_seg(res['segment'])
        plot_similar_segment_3d(ax, std_ranges)
        print K, std_ranges[0][1] - std_ranges[0][0]

    for axis in 'xyz':
        ax.tick_params(axis=axis, labelsize=5)
    ax.set_xlabel('Sequence 1')
    ax.set_ylabel('Sequence 2')
    ax.set_zlabel('Sequence 3')

    ax.set_title('estimated similarity at exactly matching %d-mers' % wordlen)

    fig.tight_layout()
    savefig(fig, 'multiple-sequence.png', dpi=300)


# ========================================================
# Biological data: compare aligned genomes
# ========================================================
def seeds_from_maf(maf_path, wordlen, ids_of_interest=[]):
    alignments = list(AlignIO.parse(maf_path, 'maf'))
    # NOTE trust the first alignment to have all the ids
    ids = set()
    for alignment in alignments:
        ids = ids.union(set(rec.id for rec in alignment))
    assert all(id_ in ids for id_ in ids_of_interest)
    ids = ids_of_interest if ids_of_interest else ids
    seqs = {id_: '' for id_ in ids}
    for alignment in alignments:
        updated = {id_: False for id_ in ids}
        line_len = len(alignment[0])
        for rec in alignment:
            if rec.id not in ids:
                continue
            seqs[rec.id] += ''.join(rec.upper())
            updated[rec.id] = True
        for id_ in seqs:
            if not updated[id_]:
                seqs[id_] += '-' * line_len

    seq_lens = set(len(seq) for seq in seqs.values())
    assert len(seq_lens) == 1, \
        'all aligned sequences must have the same length'
    seqs = np.array([list(seqs[id_]) for id_ in ids])
    num_seqs, seq_len = seqs.shape
    pos = np.zeros(num_seqs)
    for idx in range(seq_len - wordlen):
        for i in range(num_seqs):
            if seqs[i, idx] != '-':
                pos[i] += 1
        if all(len(set(seqs[i, idx + j] for i in range(num_seqs))) == 1
               for j in range(wordlen)):
            yield tuple(pos)


def exp_biological_multiple_sequences():
    maf_path = 'data/actb/actb-7vet.maf'
    wordlen = 6
    A = Alphabet('ACGT')
    WB_kw = {'g_max': .4, 'sensitivity': .9, 'alphabet': A,
             'wordlen': wordlen, 'path': ':memory:'}

    # 3 sequences for scatter plot
    ids = ['hg38.chr7',
           'panTro4.chr7',
           'canFam3.chr6',
           ]
    seqs = [A.parse(seq.upper())
            for id_, seq in get_seqs_from_mse(maf_path, fmt='maf')
            if id_ in ids]
    WB = WordBlotMultiple(*seqs, **WB_kw)
    p_min = .8
    scored_seeds = [(WB.to_ij_coordinates(*rec['seed']), rec['p'])
                    for rec in WB.score_seeds(50)]
    log('found %d seeds for %d sequences' % (len(scored_seeds), len(ids)))

    fig = plt.figure(figsize=(10, 5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[2, 2])
    ax_scatter = plt.subplot(gs[0], projection=Axes3D.name)
    plot_scored_seeds_3d(fig, ax_scatter, scored_seeds, threshold=p_min)
    for axis in 'xyz':
        ax_scatter.tick_params(axis=axis, labelsize=4)
    ax_scatter.set_xlabel(ids[0].split('.')[0])
    ax_scatter.set_ylabel(ids[1].split('.')[0])
    ax_scatter.set_zlabel(ids[2].split('.')[0])
    ax_scatter.set_title('estimated similarity at exactly matching %d-mers' %
                         wordlen, fontsize=8)

    # ============================
    # hom/non-hom seed classifier
    # ============================
    ids = ['hg38.chr7',
           'panTro4.chr7',
           'canFam3.chr6',
           'mm10.chr5',
           'rheMac3.chr3',
           'rn5.chr12',
           ]

    real_seeds = list(seeds_from_maf(maf_path, wordlen, ids_of_interest=ids))
    real_seeds = list(tuple(int(x) for x in seed) for seed in real_seeds)
    log('found %d homologous seeds for %d sequences' %
        (len(real_seeds), len(ids)))

    seqs = [A.parse(seq.upper())
            for id_, seq in get_seqs_from_mse(maf_path, fmt='maf')
            if id_ in ids]
    WB = WordBlotMultiple(*seqs, **WB_kw)
    scored_seeds = [(WB.to_ij_coordinates(*rec['seed']), rec['p'])
                    for rec in WB.score_seeds(50)]
    log('found %d seeds for %d sequences' % (len(scored_seeds), len(ids)))
    pos, neg = [], []
    for coords, p_hat in scored_seeds:
        if coords in real_seeds:
            pos.append(p_hat)
        else:
            neg.append(p_hat)

    ax_roc = plt.subplot(gs[1])
    plot_roc(ax_roc, pos, neg, color='k')
    title = 'ROC for classifing exactly matching %d-mers:' % wordlen
    title += '%d(+) %d(-) samples\n' % (len(pos), len(neg))
    title += 'species: %s' % ', '.join(x.split('.')[0] for x in ids)
    ax_roc.set_title(title, fontsize=8)

    fig.tight_layout()
    savefig(fig, 'multiple-sequence[bio].png')


# ========================================================
# Simulations (K, p, r, ds, a)
# ========================================================
@with_dumpfile
def sim_simulated_K_p(Ks, ps, **kw):
    # ps are the probability that ALL seqs match
    n_samples, n_seqs = kw['n_samples'], kw['n_seqs']

    def _zero():
        shape = (len(Ks), len(ps), n_samples)
        return {'pos': np.zeros(shape), 'neg': np.zeros(shape)}

    def _zeron():
        shape = (n_seqs - 1, len(Ks), len(ps), n_samples)
        return {'pos': np.zeros(shape), 'neg': np.zeros(shape)}

    A = Alphabet('ACGT')
    WB_kw = {
        'g_max': kw.get('g_max', .6),
        'sensitivity': kw.get('sensitivity', .99),
        'alphabet': A,
        'path': kw.get('db_path', ':memory:'),
        'log_level': kw.get('log_level', logging.WARN),
    }
    sim_data = {
        'scores': {'H0': _zero(), 'H1': _zero()},
        'K_hat': _zero(),
        'time': _zero(),
        'p_hat': _zero(),
        'd_hat': _zeron(),
        'a_hat': _zero(),
        'r_hat': _zeron(),
        'WB_kw': WB_kw,
        'Ks': Ks,
        'ps': ps,
    }
    # NOTE in multiple sequence case segments gets broken (along a) more easily
    # (this is not fragmentation which is breaking along ds). We use K / 2
    # instead.
    # K_min = 100
    # assert K_min <= min(Ks)
    p_min = .5
    assert p_min <= min(ps)

    for (K_idx, K), (p_idx, p_match) in product(enumerate(Ks), enumerate(ps)):
        # HACK adjusted wordlens
        wordlen = int(np.ceil(np.log(K)))
        WB_kw['wordlen'] = wordlen

        log('simulating (%d samples) K = %d (w=%d), p = %.2f' %
            (n_samples, K, WB_kw['wordlen'], p_match), newline=False)
        for idx in range(n_samples):
            sys.stderr.write('.')
            # distribute p_match evenly over all sequences
            p_match_indiv = np.exp(np.log(p_match) / n_seqs)
            # distribute p_match_indiv evenly over gap and subst
            subst = gap = 1 - np.sqrt(p_match_indiv)
            assert abs((1 - gap) * (1 - subst) - p_match_indiv) < 1e-3
            M = MutationProcess(A, subst_probs=subst, ge_prob=gap,
                                go_prob=gap)
            seqs_rel = [rand_seq(A, K)]
            for _ in range(1, n_seqs):
                seqs_rel.append(M.mutate(seqs_rel[0])[0])
            for seq_idx in range(n_seqs):
                seqs_rel[seq_idx] = rand_seq(A, K / 2) + \
                                    seqs_rel[seq_idx] + \
                                    rand_seq(A, K / 2)
            seqs_unrel = [rand_seq(A, 2 * K) for _ in range(n_seqs)]

            for key, seqs in zip(['pos', 'neg'], [seqs_rel, seqs_unrel]):
                t_start = time()
                WB = WordBlotMultiple(*seqs, **WB_kw)

                # calculate H0/H1 scores with perfect information:
                band_r = WB.band_radius(K)
                ds_band = [(-band_r, band_r)] * (n_seqs - 1)
                volume = K * (2 * band_r) ** (n_seqs - 1)
                num_seeds = WB.seed_count(
                    ds_band=ds_band,
                    a_band=(n_seqs * K / 2, 3 * n_seqs * K / 2)
                )
                s0, s1 = WB.score_num_seeds(num_seeds=num_seeds,
                                            volume=volume,
                                            seglen=K, p_match=p_match)
                sim_data['scores']['H0'][key][K_idx, p_idx, idx] = s0
                sim_data['scores']['H1'][key][K_idx, p_idx, idx] = s1
                sim_data['time'][key][K_idx, p_idx, idx] = time() - t_start

                def _len(seg): return (seg[1][1] - seg[1][0]) / n_seqs

                # NOTE for multiple sequences it's too much to ask for
                # at_least_one=True because unrelated sequences can have no
                # seeds.
                results = list(WB.similar_segments(K / 2, p_min))
                if results:
                    # pick the longest detected homology
                    hom = max(results, key=lambda rec: _len(rec['segment']))
                else:
                    # NOTE a fake segment with all zero properties
                    hom = {'segment': ([(0, 0)] * (n_seqs - 1), [0, 0]),
                           'p': 0}
                ds_ranges, a_range = hom['segment']
                K_hat = _len(hom['segment'])
                # FIXME K_hat doesn't grow well with K
                sim_data['K_hat'][key][K_idx, p_idx, idx] = K_hat
                sim_data['p_hat'][key][K_idx, p_idx, idx] = hom['p']
                sim_data['a_hat'][key][K_idx, p_idx, idx] = sum(a_range) / 2
                for d_idx, d_range in enumerate(ds_ranges):
                    d_min, d_max = d_range
                    d_ctr = (d_max + d_min) / 2
                    d_rad = (d_max - d_min) / 2
                    sim_data['d_hat'][key][d_idx, K_idx, p_idx, idx] = d_ctr
                    sim_data['r_hat'][key][d_idx, K_idx, p_idx, idx] = d_rad
        sys.stderr.write('\n')

    return sim_data


# FIXME make sure the definition of p_hat and p_true is the same.
def plot_simulated_K_p(sim_data, select_Ks, select_ps, suffix=''):
    Ks, ps = sim_data['Ks'], sim_data['ps']
    n_samples = sim_data['scores']['H0']['pos'].shape[2]
    n_seqs = sim_data['d_hat']['pos'].shape[0] + 1
    scores = sim_data['scores']

    kw = {'marker': 'o', 'markersize': 3, 'alpha': .7, 'lw': 1}

    # ======================================
    # varying K for select ps
    fig_by_K = plt.figure(figsize=(11, 5))
    ax_H0 = fig_by_K.add_subplot(1, 2, 1)
    ax_H1 = fig_by_K.add_subplot(1, 2, 2)
    for mode, ax in zip(['H0', 'H1'], [ax_H0, ax_H1]):
        colors = color_code(select_ps)
        for p, color in zip(select_ps, colors):
            p_idx = ps.index(p)
            for case, ls in zip(['pos', 'neg'], ['-', '--']):
                label = '' if case == 'neg' else 'p = %.2f' % p
                plot_with_sd(ax, Ks, scores[mode][case][:, p_idx, :], axis=1,
                             color=color, ls=ls, label=label, **kw)
        ax.set_ylabel('%s score' % mode, fontsize=10)
        ax.set_xlabel('similarity length', fontsize=10)
        ax.set_xticks(Ks)
        ax.set_xticklabels(Ks, rotation=90, fontsize=6)
    ax_H0.legend(loc='best', fontsize=10)
    ax_H1.legend(loc='best', fontsize=10)
    fig_by_K.suptitle('no. seqs = %d, no. samples = %d' % (n_seqs, n_samples),
                      fontsize=10)
    fig_by_K.tight_layout(rect=[0, 0.03, 1, 0.95])
    savefig(fig_by_K, 'simulations_multiple[score-by-K]%s.png' % suffix)

    # ======================================
    # varying p for select Ks
    fig_by_p = plt.figure(figsize=(11, 5))
    ax_H0 = fig_by_p.add_subplot(1, 2, 1)
    ax_H1 = fig_by_p.add_subplot(1, 2, 2)
    for mode, ax in zip(['H0', 'H1'], [ax_H0, ax_H1]):
        colors = color_code(select_Ks)
        for K, color in zip(select_Ks, colors):
            K_idx = Ks.index(K)
            for case, ls in zip(['pos', 'neg'], ['-', '--']):
                label = '' if case == 'neg' else 'K = %d' % K
                plot_with_sd(ax, ps, scores[mode][case][K_idx, :, :], axis=1,
                             color=color, ls=ls, label=label, **kw)
        ax.set_ylabel('%s score' % mode, fontsize=10)
        ax.set_xlabel('similarity match probability', fontsize=10)
        ax.set_xticks(ps)
        ax.set_xticklabels(ps, rotation=90, fontsize=6)
    ax_H0.legend(loc='best', fontsize=10)
    ax_H1.legend(loc='best', fontsize=10)
    fig_by_p.suptitle('no. seqs = %d, no. samples = %d' % (n_seqs, n_samples),
                      fontsize=10)
    fig_by_p.tight_layout(rect=[0, 0.03, 1, 0.95])
    savefig(fig_by_p, 'simulations_multiple[score-by-p]%s.png' % suffix)

    # ======================================
    K_hats = sim_data['K_hat']
    p_hats = sim_data['p_hat']
    # estimated Ks for select ps
    fig_hat = plt.figure(figsize=(11, 5))
    ax_K = fig_hat.add_subplot(1, 2, 1)
    ax_p = fig_hat.add_subplot(1, 2, 2)
    colors = color_code(select_ps)
    for p, color in zip(select_ps, colors):
        p_idx = ps.index(p)
        plot_with_sd(ax_K, Ks, K_hats['pos'][:, p_idx, :], axis=1,
                     color=color, label='p = %.2f' % p, **kw)
        plot_with_sd(ax_K, Ks, K_hats['neg'][:, p_idx, :], axis=1,
                     color=color, ls='--', **kw)
        ax_K.set_xticks(Ks)
        ax_K.set_xticklabels(Ks, rotation=90, fontsize=6)
    ax_K.set_ylabel('estimated similarity length', fontsize=10)
    ax_K.set_xlabel('true similarity length', fontsize=10)
    ax_K.plot(Ks, Ks, ls='--', c='k', lw=3, alpha=.4)
    ax_K.legend(loc='upper left', fontsize=10)

    # estimated ps for select Ks
    colors = color_code(select_Ks)
    for K, color in zip(select_Ks, colors):
        K_idx = Ks.index(K)
        plot_with_sd(ax_p, ps, p_hats['pos'][K_idx, :, :], axis=1,
                     color=color, label='K = %d' % K, **kw)
        plot_with_sd(ax_p, ps, p_hats['neg'][K_idx, :, :], axis=1,
                     color=color, ls='--', **kw)
        ax_p.set_xticks(ps)
        ax_p.set_xticklabels(ps, rotation=90, fontsize=6)
    ax_p.set_ylabel('estimated match probability', fontsize=10)
    ax_p.set_xlabel('true match probability', fontsize=10)
    ax_p.plot(ps, ps, ls='--', c='k', lw=3, alpha=.4)
    ax_p.set_ylim(0, 1)
    ax_p.legend(loc='upper left', fontsize=10)

    fig_hat.suptitle('no. seqs = %d, no. samples = %d' % (n_seqs, n_samples),
                     fontsize=10)
    fig_hat.tight_layout(rect=[0, 0.03, 1, 0.95])
    savefig(fig_hat, 'simulations_multiple[K,p-hat]%s.png' % suffix)

    # ======================================
    # estimated diagonal and antidiagonal position and band radius for select
    # match probabilities (select_ps), as a function of K
    d_hats = sim_data['d_hat']['pos']
    a_hats = sim_data['a_hat']['pos']
    r_hats = sim_data['r_hat']['pos']
    fig_hat = plt.figure(figsize=(13, 5))
    ax_r = fig_hat.add_subplot(1, 3, 1)
    ax_d = fig_hat.add_subplot(1, 3, 2)
    ax_a = fig_hat.add_subplot(1, 3, 3)

    colors = color_code(select_ps)
    for p, color in zip(select_ps, colors):
        kw = {'color': color, 'alpha': .3, 'lw': 1, 'label': 'p = %.2f' % p,
              'marker': 'o', 'markersize': 2}
        truth_kw = {'ls': '--', 'alpha': .1, 'color': 'k', 'lw': 3}
        p_idx = ps.index(p)

        plot_with_sd(ax_a, Ks, a_hats[:, p_idx, :], axis=1, **kw)

        for idx in range(n_seqs - 1):
            plot_with_sd(ax_r, Ks, r_hats[idx, :, p_idx, :], axis=1, **kw)
            plot_with_sd(ax_d, Ks, d_hats[idx, :, p_idx, :], axis=1, **kw)
            if 'label' in kw:
                kw.pop('label')
        for ax in [ax_r, ax_d, ax_a]:
            ax.set_xticks(Ks)
            ax.set_xticklabels(Ks, rotation=90, fontsize=6)

    g_max = sim_data['WB_kw']['g_max']
    sensitivity = sim_data['WB_kw']['sensitivity']

    ax_r.set_xlabel('similarity length')
    ax_r.set_ylabel('estimated diagonal band widths of similarity')
    ax_r.plot(Ks, [band_radius(K, g_max, sensitivity) for K in Ks], **truth_kw)
    ax_r.legend(loc='best', fontsize=8)
    ax_r.set_ylim(0, min(Ks))

    ax_d.set_xlabel('similarity length')
    ax_d.set_ylabel('estimatedd diagonal positions of similarity')
    ax_d.plot(Ks, [0] * len(Ks), **truth_kw)
    ax_d.legend(loc='best', fontsize=8)
    ax_d.set_ylim(-min(Ks) / 2, min(Ks) / 2)

    ax_a.plot(Ks, [n_seqs * K for K in Ks], **truth_kw)
    ax_a.set_xlabel('similarity length')
    ax_a.set_ylabel('estimated antidiagonal position of similarity')
    ax_a.legend(loc='best', fontsize=8)

    savefig(fig_hat, 'simulations_multiple[r,d,a-hat][by-K]%s.png' % suffix)

    # ======================
    # time to score seeds as a function of K (for select ps) and p (for select
    # Ks)
    fig_t = plt.figure(figsize=(8, 5))
    ax_t_K = fig_t.add_subplot(1, 2, 1)
    ax_t_p = fig_t.add_subplot(1, 2, 2)

    kw = {'alpha': .5, 'lw': 1, 'marker': 'o', 'markersize': 2}

    # times as a function of p for select Ks
    colors = color_code(select_Ks)
    for K, color in zip(select_Ks, colors):
        K_idx = Ks.index(K)
        label = 'K = %d (wordlen = %d)' % (K, int(np.ceil(np.log(K))))
        plot_with_sd(ax_t_p, ps, 1000 * sim_data['time']['pos'][K_idx, :, :],
                     axis=1, color=color, ls='-', label=label, **kw)
    ax_t_p.set_xticks(ps)
    ax_t_p.set_xticklabels(ps, rotation=90, fontsize=6)
    # NOTE make room for legend manually
    ax_t_p.set_ylim(None, 400)
    ax_t_p.set_ylabel('time (ms) to score all seeds', fontsize=10)
    ax_t_p.set_xlabel('true match probability', fontsize=10)
    ax_t_p.legend(loc='upper left', fontsize=8)

    # times as a function of K for select ps
    colors = color_code(select_ps)
    for p_match, color in zip(select_ps, colors):
        p_idx = ps.index(p_match)
        label = 'p = %.2f' % p_match
        plot_with_sd(ax_t_K, Ks, 1000 * sim_data['time']['pos'][:, p_idx, :],
                     axis=1, color=color, ls='-', label=label, **kw)
    ax_t_K.set_xticks(Ks)
    ax_t_K.set_xticklabels(Ks, rotation=90, fontsize=6)
    # NOTE make room for legend manually
    ax_t_K.set_ylim(None, 400)
    ax_t_K.set_ylabel('time (ms) to score all seeds', fontsize=10)
    ax_t_K.set_xlabel('similarity lengths (half sequence length)', fontsize=10)
    ax_t_K.legend(loc='upper left', fontsize=8)

    fig_t.suptitle('no. seqs = %d, no. samples = %d' % (n_seqs, n_samples),
                   fontsize=10)
    fig_t.tight_layout(rect=[0, 0.03, 1, 0.95])
    savefig(fig_t, 'simulations_multiple[times]%s.png' % suffix)


def exp_simulated_K_p():
    Ks = [200 * i for i in range(1, 9)]
    select_Ks = Ks[1], Ks[3], Ks[5]

    ps = [round(1 - .06 * i, 2) for i in range(1, 8)]
    select_ps = ps[0], ps[3], ps[5]

    n_samples = 50
    n_seqs = 6

    suffix = ''
    dumpfile = 'simulations_multiple%s.txt' % suffix
    sim_data = sim_simulated_K_p(
        Ks, ps, n_samples=n_samples, n_seqs=n_seqs,
        dumpfile=dumpfile, ignore_existing=False)
    plot_simulated_K_p(sim_data, select_Ks, select_ps, suffix=suffix)


# ========================================================
# Simulated Rearrangements
# ========================================================
@with_dumpfile
def sim_rearranged_sequences(**kw):
    A = Alphabet('ACGT')
    region_len, subst, gap = kw['region_len'], kw['subst'], kw['gap']
    n_regions, n_indivs = kw['n_regions'], kw['n_indivs']
    wordlen = kw['wordlen']
    M = MutationProcess(A, subst_probs=subst, ge_prob=gap, go_prob=gap)

    # we represnet a sequence by a list of region copies. In each generation
    # random mutations are applied to each region copy and a random pair of
    # regions swap place.
    def _evolve(region_copies):
        new_gen = [None] * len(region_copies)
        for idx, region in enumerate(region_copies):
            new_gen[idx], _ = M.mutate(region)
            new_gen[idx], _ = M.mutate(region)
        i, j = np.random.choice(range(n_regions), 2, replace=False)
        new_gen[i], new_gen[j] = region_copies[j], region_copies[i]
        return new_gen

    generations = [[rand_seq(A, region_len) for _ in range(n_regions)]]
    for _ in range(n_indivs - 1):
        generations.append(_evolve(generations[-1]))
    # each individual is a sequence of region copies
    seqs = [sum(individual, A.parse('')) for individual in generations]

    # NOTE if g is too small we get massive fragmentation that grows
    # exponentially with the number of sequences. The same presumably happens
    # in pairwaise but less noticable.
    # if g is too large we avoid fragmentation but we start not seeing
    # similarities because volume increases (and then what?)
    WB_kw = {'g_max': .6, 'sensitivity': .9, 'alphabet': A,
             'wordlen': wordlen, 'path': ':memory:'}

    K_min = .9 * region_len
    p_min = .6
    sim_data = {
        'n_indivs': n_indivs,
        'n_regions': n_regions,
        'region_len': region_len,
        'generations': generations,
        'seqs': seqs,
        'WB_kw': WB_kw,
        'p_min': p_min,
        'K_min': K_min,
        'segments': [],
    }

    WB = WordBlotMultiple(*seqs, **WB_kw)
    for res in WB.similar_segments(K_min=K_min, p_min=p_min, score=False):
        std_ranges = WB.to_ij_coordinates_seg(res['segment'])
        p_hat = round(res['p'], 3)
        sim_data['segments'].append({'p': p_hat, 'segment': std_ranges})

    return sim_data


def plot_rearranged_sequences(sim_data):
    seqs = sim_data['seqs']
    n_indivs = len(seqs)
    width = .5 * n_indivs
    height = sim_data['n_regions'] * sim_data['region_len'] / 400.
    fig = plt.figure(figsize=(width, height))
    ax = fig.add_subplot(1, 1, 1)

    colors = color_code(range(len(sim_data['segments'])), cmap='tab20')

    for color, rec in zip(colors, sim_data['segments']):
        seg = rec['segment']
        ctrs = [sum(r) / 2 for r in seg]
        xs = range(n_indivs)
        for idx, x in enumerate(xs):
            # ys = [seg[idx][0], sum(seg[idx]) / 2, seg[idx][1]]
            ax.plot([x, x], seg[idx], c=color, lw=2, alpha=.4)
        ax.plot(xs, ctrs, lw=1, marker='o', markersize=5, c=color, alpha=.6)

    for idx, indiv in enumerate(sim_data['generations']):
        positions = np.cumsum([len(seq) for seq in indiv])
        for pos in positions:
            ax.plot([idx - .02, idx + .02], [pos, pos], c='k', alpha=.7)

    ax.set_xticks(range(n_indivs))
    ax.set_xticklabels(['seq %d' % (i + 1) for i in range(n_indivs)],
                       fontsize=8)
    ax.set_xlim(-1, n_indivs)

    fig.tight_layout()
    savefig(fig, 'multiple_rearrangement.png')


def exp_rearrangement():
    suffix = ''
    dumpfile = 'multiple_rearangements%s.txt' % suffix
    subst = .03
    gap = .02

    # NOTE adjust wordlen according to sequence length; good pairs: region_len
    # (100), wordlen (10). region_len(1000), wordlen(12)
    region_len = 200
    n_regions = 10
    wordlen = 8
    n_indivs = 10
    sim_data = sim_rearranged_sequences(
        region_len=region_len, n_regions=n_regions, n_indivs=n_indivs,
        subst=subst, gap=gap, wordlen=wordlen,
        dumpfile=dumpfile, ignore_existing=True
    )
    plot_rearranged_sequences(sim_data)


if __name__ == '__main__':
    exp_three_syntehtic_sequences()
    exp_biological_multiple_sequences()
    exp_rearrangement()
    exp_simulated_K_p()
