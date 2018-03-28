import numpy as np
from matplotlib import pyplot as plt
from itertools import product
from scipy.ndimage.filters import gaussian_filter1d

import logging

from biseqt.blot import WordBlot, find_peaks, band_radii, band_radius
from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess

from util import plot_with_sd, color_code, plot_classifier
from util import with_dumpfile, log, savefig, load_fasta
from util import seq_pair

from util import plot_global_alignment, adjust_pw_plot
from util import estimate_match_probs_in_opseq, fill_in_unknown
from util import plot_scored_seeds, seeds_from_opseq
from util import opseq_path, plot_local_alignment
from biseqt.pw import Aligner, STD_MODE, LOCAL


@with_dumpfile
def sim_stats_varying_K_p(Ks, ps, n_samples, **kw):
    def _zero():
        shape = (len(Ks), len(ps), n_samples)
        return {'pos': np.zeros(shape), 'neg': np.zeros(shape)}

    A = Alphabet('ACGT')
    wordlen = kw['wordlen']
    WB_kw = {
        'g_max': kw.get('g_max', .6),
        'sensitivity': kw.get('sensitivity', .99),
        'wordlen': wordlen,
        'alphabet': A,
        'path': kw.get('db_path', ':memory:'),
        'log_level': kw.get('log_level', logging.WARNING),
    }
    sim_data = {
        'scores': {'H0': _zero(), 'H1': _zero()},
        'K_hat': _zero(),
        'p_hat': _zero(),
        'd_hat': _zero(),
        'a_hat': _zero(),
        'r_hat': _zero(),
        'WB_kw': WB_kw,
        'Ks': Ks,
        'ps': ps,
    }
    K_min = 100
    p_min = .5
    assert K_min <= min(Ks)
    assert p_min <= min(ps)

    for (K_idx, K), (p_idx, p_match) in product(enumerate(Ks), enumerate(ps)):
        log('evaluating scores for K = %d, p = %.2f' % (K, p_match))
        for idx in range(n_samples):
            # distribute p_match evenly over gap and subst
            subst = gap = 1 - np.sqrt(p_match)
            assert abs((1 - gap) * (1 - subst) - p_match) < 1e-3
            M = MutationProcess(A, subst_probs=subst, ge_prob=gap,
                                go_prob=gap)
            S_rel, T_rel = seq_pair(K, A, mutation_process=M)
            S_rel = rand_seq(A, K / 2) + S_rel + rand_seq(A, K / 2)
            T_rel = rand_seq(A, K / 2) + T_rel + rand_seq(A, K / 2)
            S_urel, T_urel = rand_seq(A, K), rand_seq(A, K)

            for key, (S, T) in zip(['pos', 'neg'],
                                   [(S_rel, T_rel), (S_urel, T_urel)]):
                WB = WordBlot(S, T, **WB_kw)

                def _len(seg): return seg[1][1] - seg[1][0]

                results = list(WB.similar_segments(K_min, p_min,
                                                   at_least_one=True))
                # pick the longest detected homology
                sim_data['K_hat'][key][K_idx][p_idx][idx] = max(
                    _len(rec['segment']) for rec in results
                )
                hom = max(results, key=lambda rec: _len(rec['segment']))
                (d_min, d_max), (a_min, a_max) = hom['segment']
                s0, s1 = hom['scores']
                sim_data['p_hat'][key][K_idx][p_idx][idx] = hom['p']
                sim_data['d_hat'][key][K_idx][p_idx][idx] = (d_min + d_max) / 2
                sim_data['a_hat'][key][K_idx][p_idx][idx] = (a_min + a_max) / 2
                sim_data['r_hat'][key][K_idx][p_idx][idx] = (d_max - d_min) / 2
                sim_data['scores']['H0'][key][K_idx][p_idx][idx] = s0
                sim_data['scores']['H1'][key][K_idx][p_idx][idx] = s1

    return sim_data


@with_dumpfile
def sim_stats_real_homologies(seqs, pws, **kw):
    A = Alphabet('ACGT')
    wordlen, K_min, p_min = kw['wordlen'], kw['K_min'], kw['p_min']
    WB_kw = {
        'g_max': kw.get('g_max', .6),
        'sensitivity': kw.get('sensitivity', .99),
        'wordlen': wordlen,
        'alphabet': A,
        'path': kw.get('db_path', ':memory:'),
        'log_level': kw.get('log_level', logging.WARNING),
        # HACK mask CG-rich and homopolymeric regions
        # 'mask': [set(x) for x in [[0], [1], [2], [3], [1, 2], [0, 3]]],
    }
    qM = 1
    qS = qG = -1
    aligner_kw = {
        'match_score': qM,
        'mismatch_score': qS,
        'ge_score': qG,
        'go_score': 0,
        'alnmode': STD_MODE,
        'alntype': LOCAL,
        # 'alnmode': BANDED_MODE,
        # 'alntype': B_GLOBAL,
    }

    similar_segments_kw = {'K_min': K_min, 'p_min': p_min}
    sim_data = {
        'pws': pws,
        'seqlens': {name: len(seqs[name]) for name in seqs},
        'seed_pos': [],
        'seed_neg': [],
        'seeds': {key: {} for key in pws},
        'opseq_seeds': {key: [] for key in pws},
        'homologous_tpr': {key: 0 for key in pws},
        'segments': {key: {} for key in pws},
        'WB_kw': WB_kw,
        'similar_segments_kw': similar_segments_kw,
        'aligner_kw': aligner_kw,
    }
    for idx, ((id1, id2), opseq) in enumerate(pws.items()):
        log('finding local homologies between %s (%d) and %s (%d)' %
            (id1, sim_data['seqlens'][id1], id2, sim_data['seqlens'][id2]))
        S = seqs[id1]
        T = seqs[id2]

        WB = WordBlot(S, T, **WB_kw)
        sim_data['seeds'][(id1, id2)] = {
            WB.to_ij_coordinates(*rec['seed']): rec['p']
            for rec in WB.score_seeds(K_min)
        }
        log('-> found %d exactly matching %d-mers' %
            (len(sim_data['seeds'][(id1, id2)]), wordlen))

        # separate +/- seeds with their estimated probabilities
        band_r = band_radius(K_min, sim_data['WB_kw']['g_max'],
                             sim_data['WB_kw']['sensitivity'])
        in_band = np.zeros((len(S), len(T)))
        opseq_seeds = list(seeds_from_opseq(opseq, wordlen))
        for _, seed in opseq_seeds:
            d, a = WB.to_diagonal_coordinates(*seed)
            for d_ in range(d - band_r, d + band_r):
                # FIXME what should the range be?
                for a_ in range(a - band_r, a + band_r):
                    i_, j_ = WB.to_ij_coordinates(d_, a_)
                    if 0 <= i_ < len(S) and 0 <= j_ < len(T):
                        in_band[i_, j_] = 1
        sim_data['opseq_seeds'][(id1, id2)] = opseq_seeds

        for seed, p in sim_data['seeds'][(id1, id2)].items():
            if in_band[seed[0], seed[1]]:
                sim_data['seed_pos'].append(p)
            else:
                sim_data['seed_neg'].append(p)
        log('separated +/- seeds')

        sim_data['segments'][(id1, id2)] = list(
            WB.similar_segments(**similar_segments_kw)
        )
        log('found %d similar segments' %
            len(sim_data['segments'][(id1, id2)]))

        # separate +/- coords
        probs = estimate_match_probs_in_opseq(opseq, wordlen)
        xs, ys = opseq_path(opseq, x0=0, y0=0)
        segments_on_opseq = find_peaks(probs, 5, p_min)
        hom_coords = np.zeros((len(S) + 1, len(T) + 1))
        for start, end in segments_on_opseq:
            for aln_pos in range(start, end + 1):
                i, j = xs[aln_pos], ys[aln_pos]
                hom_coords[i, j] = 1

        in_band_segs = np.zeros((len(S) + 1, len(T) + 1))
        for seg_info in sim_data['segments'][(id1, id2)]:
            (d_min, d_max), (a_min, a_max) = seg_info['segment']
            for d, a in product(range(d_min, d_max), range(a_min, a_max)):
                i, j = WB.to_ij_coordinates(d, a)
                if 0 <= i < len(S) and 0 <= j < len(T):
                    in_band_segs[i, j] = 1
        true_positive = np.sum(in_band_segs * hom_coords)
        positive = np.sum(hom_coords)
        sim_data['homologous_tpr'][(id1, id2)] = true_positive / positive
        log('tpr = %.2f' % sim_data['homologous_tpr'][(id1, id2)])

        # actual alignments
        for idx, seg_info in enumerate(sim_data['segments'][(id1, id2)]):
            (d_min, d_max), (a_min, a_max) = seg_info['segment']
            corners = [WB.to_ij_coordinates(d, a)
                       for d, a in product(*seg_info['segment'])]
            i_start = min(i for i, _ in corners)
            j_start = min(j for _, j in corners)
            i_end = max(i for i, j in corners)
            j_end = max(j for _, j in corners)
            S_ = S[i_start: i_end]
            T_ = T[j_start: j_end]
            aligner = Aligner(S_, T_, **aligner_kw)
            with aligner:
                aligner.solve()
                alignment = aligner.traceback()
                sim_data['segments'][(id1, id2)][idx]['alignment'] = alignment
                sim_data['segments'][(id1, id2)][idx]['frame'] = \
                    (i_start, i_end), (j_start, j_end)
                # if alignment is not None:
                #     tx = alignment.transcript
                #     print seg_info['p'], 1. * tx.count('M') / len(tx)
    return sim_data


def plot_stats_varying_K_p(sim_data, select_Ks, select_ps, suffix=''):
    Ks, ps = sim_data['Ks'], sim_data['ps']
    wordlen = sim_data['WB_kw']['wordlen']
    n_samples = sim_data['scores']['H0']['pos'].shape[2]
    scores = sim_data['scores']

    assert all(K in Ks for K in select_Ks)
    assert all(p in ps for p in select_ps)

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
    ax_H0.legend(loc='best', fontsize=10)
    ax_H1.legend(loc='best', fontsize=10)
    fig_by_K.suptitle('w = %d, no. samples = %d' % (wordlen, n_samples),
                      fontsize=10)
    fig_by_K.tight_layout(rect=[0, 0.03, 1, 0.95])
    savefig(fig_by_K, 'stats[score-by-K]%s.png' % suffix)

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
    ax_H0.legend(loc='best', fontsize=10)
    ax_H1.legend(loc='best', fontsize=10)
    fig_by_p.suptitle('w = %d, no. samples = %d' % (wordlen, n_samples),
                      fontsize=10)
    fig_by_p.tight_layout(rect=[0, 0.03, 1, 0.95])
    savefig(fig_by_p, 'stats[score-by-p]%s.png' % suffix)

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
    ax_p.set_ylabel('estimated match probability', fontsize=10)
    ax_p.set_xlabel('true match probability', fontsize=10)
    ax_p.plot(ps, ps, ls='--', c='k', lw=3, alpha=.4)
    ax_p.set_ylim(0, 1)
    ax_p.legend(loc='upper left', fontsize=10)

    fig_hat.suptitle('w = %d, no. samples = %d' % (wordlen, n_samples),
                     fontsize=10)
    fig_hat.tight_layout(rect=[0, 0.03, 1, 0.95])
    savefig(fig_hat, 'stats[K,p-hat]%s.png' % suffix)

    # ======================================
    # estimated diagonal and antidiagonal position and band radius for select
    # match probabilities (select_ps), as a function of K
    d_hats = sim_data['d_hat']['pos']
    a_hats = sim_data['a_hat']['pos']
    r_hats = sim_data['r_hat']['pos']
    fig_hat = plt.figure(figsize=(12, 4))
    ax_r = fig_hat.add_subplot(1, 3, 1)
    ax_d = fig_hat.add_subplot(1, 3, 2)
    ax_a = fig_hat.add_subplot(1, 3, 3)

    colors = color_code(select_ps)
    for p, color in zip(select_ps, colors):
        kw = {'color': color, 'alpha': .9, 'lw': 1, 'label': 'p = %.2f' % p}
        truth_kw = {'ls': '--', 'alpha': .1, 'color': 'k', 'lw': 3}
        p_idx = ps.index(p)

        plot_with_sd(ax_r, Ks, r_hats[:, p_idx, :], axis=1, **kw)
        ax_r.set_ylim(-min(Ks), min(Ks))
        ax_r.set_xlabel('similarity length K')
        ax_r.set_ylabel('estimated diagonal band width of similarity')
        g_max = sim_data['WB_kw']['g_max']
        sensitivity = sim_data['WB_kw']['sensitivity']
        ax_r.plot(Ks, band_radii(Ks, g_max, sensitivity), **truth_kw)
        ax_r.legend(loc='best', fontsize=8)

        plot_with_sd(ax_d, Ks, d_hats[:, p_idx, :], axis=1, **kw)
        ax_d.set_ylim(-min(Ks), min(Ks))
        ax_d.set_xlabel('similarity length K')
        ax_d.set_ylabel('estimated diagonal position of similarity')
        ax_d.plot(Ks, [0] * len(Ks), **truth_kw)
        ax_d.legend(loc='best', fontsize=8)

        plot_with_sd(ax_a, Ks, a_hats[:, p_idx, :], axis=1, **kw)
        ax_a.plot(Ks, Ks, **truth_kw)
        ax_a.set_xlabel('similarity length K')
        ax_a.set_ylabel('estimated antidiagonal position of similarity')
        ax_a.legend(loc='best', fontsize=8)
    savefig(fig_hat, 'stats[r,d,a-hat]%s.png' % suffix)


def plot_stats_real_homologies(sim_data, suffix='', naming_style=None):
    seqlens = sim_data['seqlens']
    pws = sim_data['pws']
    K_min = sim_data['similar_segments_kw']['K_min']
    fig_num = int(np.ceil(np.sqrt(len(pws))))

    fig_seeds = plt.figure(figsize=(6 * fig_num, 5 * fig_num))
    fig_profiles = plt.figure(figsize=(8 * fig_num, 4 * fig_num))
    fig_p = plt.figure(figsize=(6, 4))
    fig_p_aln = plt.figure(figsize=(6, 4))
    fig_K_aln = plt.figure(figsize=(6, 4))
    fig_coord_classifier = plt.figure(figsize=(6, 4))

    ax_coord_classifier = fig_coord_classifier.add_subplot(1, 1, 1)

    ax_p = fig_p.add_subplot(1, 1, 1)

    ax_p_aln = fig_p_aln.add_subplot(1, 1, 1)

    ax_K_aln = fig_K_aln.add_subplot(1, 1, 1)

    p_min = sim_data['similar_segments_kw']['p_min']

    labels = []
    for idx, ((id1, id2), opseq) in enumerate(pws.items()):
        key = (id1, id2)

        if naming_style == 'ensembl':
            id1 = id1.split('/')[0].replace('_', ' ')
            id2 = id2.split('/')[0].replace('_', ' ')
        elif naming_style == 'ucsc':
            id1 = id1.split('.')[0]
            id2 = id2.split('.')[0]

        labels.append((id1, id2))

        # ============================
        # match probability estimation
        # ============================
        radius = K_min / 2

        ax_profiles = fig_profiles.add_subplot(fig_num, fig_num, idx + 1)
        ps_true = estimate_match_probs_in_opseq(opseq, radius)
        # smooth for ease of visual inspection
        ps_true_smooth = gaussian_filter1d(ps_true, 10)
        ax_profiles.plot(range(len(ps_true)), ps_true_smooth, c='g', alpha=.8,
                         lw=1)
        ax_profiles.set_title('%s vs %s' % (id1, id2))
        ax_profiles.set_xlabel('position along global alignment', fontsize=4)
        ax_profiles.set_ylabel('match probability', fontsize=4)

        ps_hat_pos, ps_hat = [], []
        for pos, seed in sim_data['opseq_seeds'][key]:
            if seed not in sim_data['seeds'][key]:
                # some seeds are masked
                continue
            if ps_true[pos] == 0:
                # flanking zero regions in estimated gap probabilities
                continue
            ps_hat.append(sim_data['seeds'][key][seed])
            ps_hat_pos.append(pos)
        ax_profiles.scatter(ps_hat_pos, ps_hat, lw=0, c='k', s=4, alpha=.4)

        ax_p.scatter([ps_hat[i] for i in range(len(ps_hat))],
                     [ps_true[ps_hat_pos[i]] for i in range(len(ps_hat))],
                     lw=0, color='g', s=4, alpha=.2)

        # =============
        # Dot Plots
        # =============
        ax_seeds = fig_seeds.add_subplot(fig_num, fig_num, idx + 1)
        plot_scored_seeds(ax_seeds,
                          sim_data['seeds'][key].items(),
                          threshold=p_min, alpha=.7, zorder=9)
        plot_global_alignment(ax_seeds, opseq, lw=7, alpha=.4, color='g')
        adjust_pw_plot(ax_seeds, seqlens[key[0]], seqlens[key[1]])
        ax_seeds.set_ylabel(id1, fontsize=10)
        ax_seeds.set_xlabel(id2, fontsize=10)

        # ================================
        # match probability estimation vs local alignments
        # ================================
        for seg_info in sim_data['segments'][key]:
            p_hat = seg_info['p']
            alignment = seg_info['alignment']
            K_hat = seg_info['segment'][1][1] - seg_info['segment'][1][0]
            if alignment is None:
                p_aln = 0
                K_aln = 0
            else:
                transcript = alignment.transcript
                p_aln = 1. * transcript.count('M') / len(transcript)
                plot_local_alignment(ax_seeds, transcript,
                                     seg_info['frame'][0][0],
                                     seg_info['frame'][1][0],
                                     lw=7, alpha=.2, color='b')
                K_aln = len(transcript)
            ax_p_aln.scatter([p_hat], [p_aln], lw=0, color='b', s=10, alpha=.5)
            ax_K_aln.scatter([K_hat], [K_aln], lw=0, color='b', s=10, alpha=.5)

    tprs = [sim_data['homologous_tpr'][key_] for key_ in pws]
    ax_coord_classifier.bar(range(len(labels)), tprs, color='g', alpha=.9,
                            align='center', width=.4)
    ax_coord_classifier.set_xticks(range(len(pws)))
    ax_coord_classifier.set_xticklabels([' vs. '.join(x) for x in labels],
                                        rotation=30, fontsize=4, ha='right')
    ax_coord_classifier.set_ylabel('true positive rate')
    ax_coord_classifier.set_title('homologous coordinates in similar segments')

    ax_p.plot([0, 1], [0, 1], lw=1, ls='--', alpha=.8, c='k')
    ax_p.set_xlabel('estimated match probability')
    ax_p.set_ylabel('global alignment match probability')
    ax_p.set_xlim(-.1, 1.1)
    ax_p.set_ylim(-.1, 1.1)
    ax_p.set_aspect('equal')
    ax_p.set_title('Match probability at homologous seeds')

    ax_p_aln.plot([0, 1], [0, 1], lw=1, ls='--', alpha=.8, c='k')
    ax_p_aln.set_xlabel('estimated match probability', fontsize=8)
    ax_p_aln.set_ylabel('local alignment match probability', fontsize=8)
    ax_p_aln.set_xlim(-.1, 1.1)
    ax_p_aln.set_ylim(-.1, 1.1)
    ax_p_aln.set_aspect('equal')
    ax_p_aln.set_title('Match probability at identified segments')

    min_K, max_K = ax_K_aln.get_xlim()
    ax_K_aln.plot([min_K, max_K], [min_K, max_K], lw=1, ls='--', alpha=.8,
                  c='k')
    ax_K_aln.set_xlabel('estimated similarity length', fontsize=8)
    ax_K_aln.set_ylabel('local alignment length', fontsize=8)
    ax_K_aln.set_aspect('equal')
    ax_K_aln.set_title('Similarity length at identified segments')

    plot_classifier('real_homologies[seed-classifier]%s.png' % suffix,
                    sim_data['seed_pos'], sim_data['seed_neg'],
                    labels=['homologous', 'non-homologous'])

    for fig in [fig_profiles, fig_p, fig_p_aln, fig_K_aln, fig_seeds,
                fig_coord_classifier]:
        fig.tight_layout()
    savefig(fig_profiles, 'real_homologies[profiles]%s.png' % suffix)
    savefig(fig_p, 'real_homologies[p-hat]%s.png' % suffix)
    savefig(fig_p_aln, 'real_homologies[p-hat-aln]%s.png' % suffix)
    savefig(fig_K_aln, 'real_homologies[K-hat-aln]%s.png' % suffix)
    savefig(fig_seeds, 'real_homologies[seeds]%s.png' % suffix)
    savefig(fig_coord_classifier,
            'real_homologies[coords-classifier]%s.png' % suffix)


def exp_stats_performance_varying_K_p():
    Ks = [200 * i for i in range(1, 9)]
    select_Ks = Ks[1], Ks[3], Ks[5]

    ps = [1 - .06 * i for i in range(1, 9)]
    select_ps = ps[1], ps[3], ps[5]

    n_samples = 50  # HACK
    wordlen = 5

    suffix = '[varying-K-p]'
    dumpfile = 'stats%s.txt' % suffix
    sim_data = sim_stats_varying_K_p(
        Ks, ps, n_samples, wordlen=wordlen,
        dumpfile=dumpfile, ignore_existing=False)
    plot_stats_varying_K_p(sim_data, select_Ks, select_ps, suffix=suffix)


def exp_stats_real_homologies():
    p_min = .6
    A = Alphabet('ACGT')

    # wordlen = 5
    # K_min = 100
    # style = 'ucsc'
    # suffix = '[actb][p=%.2f]' % p_min
    # dumpfile = 'real_homologies%s.txt' % suffix
    # seqs_path = 'data/actb/actb-7vet.fa'
    # pws_path = 'data/actb/actb-7vet-pws.fa'

    wordlen = 8
    K_min = 100
    style = 'ensembl'
    suffix = '[irx1][p=%.2f]' % p_min
    dumpfile = 'real_homologies%s.txt' % suffix
    seqs_path = 'data/irx1/irx1-vert-amniota-indiv.fa'
    pws_path = 'data/irx1/irx1-vert-amniota-pws.fa'

    with open(seqs_path) as f:
        seqs = {name: A.parse(fill_in_unknown(seq.upper(), A))
                # for seq, name, _ in load_fasta(f) if name in names}
                for seq, name, _ in load_fasta(f)}
    with open(pws_path) as f:
        pws = {tuple(name.split(':')): seq for seq, name, _ in load_fasta(f)}
        pws = {key: value for key, value in pws.items()
               if key[0] in seqs and key[1] in seqs}
    sim_data = sim_stats_real_homologies(
        seqs, pws, wordlen=wordlen, p_min=p_min, K_min=K_min,
        dumpfile=dumpfile, ignore_existing=False)
    plot_stats_real_homologies(sim_data, suffix=suffix, naming_style=style)


if __name__ == '__main__':
    exp_stats_performance_varying_K_p()
    exp_stats_real_homologies()
