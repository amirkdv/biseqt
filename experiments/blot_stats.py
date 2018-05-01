import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
from itertools import product
from scipy.ndimage.filters import gaussian_filter1d

import sys
import logging

from biseqt.blot import WordBlot, find_peaks, band_radius
from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess

from util import plot_with_sd, color_code, plot_classifier
from util import with_dumpfile, log, savefig, load_fasta
from util import seq_pair

from util import plot_global_alignment, adjust_pw_plot
from util import estimate_match_probs_in_opseq, fill_in_unknown
from util import plot_scored_seeds, seeds_from_opseq
from util import opseq_path, plot_local_alignment
from util import plot_cdf
from util import estimate_gap_probs_in_opseq
from biseqt.pw import Aligner, BANDED_MODE, B_GLOBAL


@with_dumpfile
def sim_simulated_K_p(Ks, ps, n_samples, **kw):
    shape = (len(Ks), len(ps), n_samples)

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
        'scores': {'H0': {'pos': np.zeros(shape), 'neg': np.zeros(shape)},
                   'H1': {'pos': np.zeros(shape), 'neg': np.zeros(shape)}},
        'K_hat': np.zeros(shape),
        'p_hat': np.zeros(shape),
        'a_hat': np.zeros(shape),
        'd_min': np.zeros(shape),
        'd_max': np.zeros(shape),
        'WB_kw': WB_kw,
        'Ks': Ks,
        'ps': ps,
    }
    K_min = 100
    assert K_min <= min(Ks)

    for (K_idx, K), (p_idx, p_match) in product(enumerate(Ks), enumerate(ps)):
        p_min = p_match

        log('simulating (%d samples) K = %d, p = %.2f' %
            (n_samples, K, p_match), newline=False)
        for idx in range(n_samples):
            sys.stderr.write('.')
            # distribute p_match evenly over gap and subst
            subst = gap = 1 - np.sqrt(p_match)
            assert abs((1 - gap) * (1 - subst) - p_match) < 1e-3
            M = MutationProcess(A, subst_probs=subst, ge_prob=gap,
                                go_prob=gap)
            S_rel, T_rel = seq_pair(K, A, mutation_process=M)
            S_rel = rand_seq(A, K / 2) + S_rel + rand_seq(A, K / 2)
            T_rel = rand_seq(A, K / 2) + T_rel + rand_seq(A, K / 2)
            S_urel, T_urel = rand_seq(A, 2 * K), rand_seq(A, 2 * K)

            for key, (S, T) in zip(['pos', 'neg'],
                                   [(S_rel, T_rel), (S_urel, T_urel)]):
                WB = WordBlot(S, T, **WB_kw)

                # calculate H0/H1 scores with perfect information:
                band_r = WB.band_radius(K)
                num_seeds = WB.seed_count(d_band=(-band_r, band_r),
                                          a_band=(K, 3 * K))
                s0, s1 = WB.score_num_seeds(num_seeds=num_seeds,
                                            area=2 * band_r * K,
                                            seglen=K, p_match=p_match)
                sim_data['scores']['H0'][key][K_idx, p_idx, idx] = s0
                sim_data['scores']['H1'][key][K_idx, p_idx, idx] = s1

                results = list(WB.similar_segments(K_min, p_min,
                                                   at_least_one=True))

                if key == 'neg':
                    continue

                # sum of K_hat, average of a_hat, d_min, d_max
                sim_data['K_hat'][K_idx, p_idx, idx] = sum(
                    r['segment'][1][1] - r['segment'][1][0]
                    for r in results) / 2
                sim_data['a_hat'][K_idx, p_idx, idx] = sum(
                    (r['segment'][1][1] + r['segment'][1][0]) / 2
                    for r in results) / len(results)
                sim_data['d_min'][K_idx, p_idx, idx] = sum(
                    r['segment'][0][0] for r in results) / len(results)
                sim_data['d_max'][K_idx, p_idx, idx] = sum(
                    r['segment'][0][1] for r in results) / len(results)

                # pick the longest detected homology for p_hat
                hom = max(
                    results,
                    key=lambda r: r['segment'][1][1] - r['segment'][1][0]
                )
                sim_data['p_hat'][key][K_idx, p_idx, idx] = hom['p']

        sys.stderr.write('\n')
    return sim_data


def plot_simulated_K_p(sim_data, select_Ks, select_ps, suffix=''):
    Ks, ps = sim_data['Ks'], sim_data['ps']
    scores = sim_data['scores']

    kw = {'marker': 'o', 'markersize': 3, 'alpha': .7, 'lw': 1}
    truth_kw = {'ls': '--', 'alpha': .6, 'color': 'k', 'lw': 1}

    # ======================================
    # score by varying K for select ps
    # score varying p for select Ks
    fig_scores = plt.figure(figsize=(8, 7))
    ax_H0_K = fig_scores.add_subplot(2, 2, 1)
    ax_H0_p = fig_scores.add_subplot(2, 2, 2)
    ax_H1_K = fig_scores.add_subplot(2, 2, 3)
    ax_H1_p = fig_scores.add_subplot(2, 2, 4)

    for mode, ax in zip(['H0', 'H1'], [ax_H0_K, ax_H1_K]):
        colors = color_code(select_ps)
        for p, color in zip(select_ps, colors):
            p_idx = ps.index(p)
            for case, ls in zip(['pos', 'neg'], ['-', '--']):
                label = '' if case == 'neg' else 'p = %.2f' % p
                plot_with_sd(ax, Ks, scores[mode][case][:, p_idx, :], axis=1,
                             color=color, ls=ls, label=label, **kw)
        ax.set_ylabel('%s score' % mode)
        ax.set_xlabel('similarity length')
        ax.set_xscale('log')
        ax.set_xticks(Ks)
        ax.set_xticklabels(Ks)
        ax.legend(loc='best')

    for mode, ax in zip(['H0', 'H1'], [ax_H0_p, ax_H1_p]):
        colors = color_code(select_Ks)
        for K, color in zip(select_Ks, colors):
            K_idx = Ks.index(K)
            for case, ls in zip(['pos', 'neg'], ['-', '--']):
                label = '' if case == 'neg' else 'K = %d' % K
                plot_with_sd(ax, ps, scores[mode][case][K_idx, :, :], axis=1,
                             color=color, ls=ls, label=label, **kw)
        ax.set_ylabel('%s score' % mode)
        ax.set_xlabel('similarity match probability')
        ax.set_xticks(ps)
        ax.set_xticklabels(ps)
        ax.legend(loc='best')
    savefig(fig_scores, 'simulations[scores]%s.png' % suffix)

    # =====================================
    fig_hats = plt.figure(figsize=(8, 8))
    ax_K = fig_hats.add_subplot(2, 2, 1)
    ax_p = fig_hats.add_subplot(2, 2, 3)
    ax_d = fig_hats.add_subplot(2, 2, 2)
    ax_a = fig_hats.add_subplot(2, 2, 4)

    K_hats = sim_data['K_hat']
    p_hats = sim_data['p_hat']
    # wordblot p_hats are with respect to projected alignment lengths whereas
    # p_true which is a property of MutationProcess is with respect to
    # alignment length; correction via p_hats *= (K/L) = (1 - g/2)
    p_hats *= (1 + p_hats) / 2

    # estimated Ks for select ps
    colors = color_code(select_ps)
    for p, color in zip(select_ps, colors):
        p_idx = ps.index(p)
        plot_with_sd(ax_K, Ks, K_hats[:, p_idx, :], axis=1,
                     color=color, label='p = %.2f' % p, **kw)
        ax_K.set_xscale('log')
        ax_K.set_yscale('log')
        ax_K.set_xticks(Ks)
        ax_K.set_xticklabels(Ks, rotation=90)
    ax_K.set_ylabel('estimated similarity length')
    ax_K.set_xlabel('true similarity length')
    ax_K.plot(Ks, Ks, **truth_kw)
    ax_K.legend(loc='upper left')

    # estimated ps for select Ks
    colors = color_code(select_Ks)
    for K, color in zip(select_Ks, colors):
        K_idx = Ks.index(K)
        plot_with_sd(ax_p, ps, p_hats[K_idx, :, :], axis=1,
                     color=color, label='K = %d' % K, **kw)
        ax_p.set_xticks(ps)
        ax_p.set_xticklabels(ps, rotation=90)
    ax_p.set_ylabel('estimated match probability')
    ax_p.set_xlabel('true match probability')
    ax_p.plot(ps, ps, **truth_kw)
    ax_p.set_ylim(0.3, 1)
    ax_p.legend(loc='lower left')

    # ======================================
    # estimated diagonal and antidiagonal position and band radius for select
    # match probabilities (select_ps), as a function of K
    d_mins = sim_data['d_min']
    d_maxs = sim_data['d_max']
    a_hats = sim_data['a_hat']

    colors = color_code(select_ps)
    for p, color in zip(select_ps, colors):
        label = 'p = %.2f' % p
        kw = {'color': color, 'alpha': .6, 'lw': 1,
              'marker': 'o', 'markersize': 3}
        p_idx = ps.index(p)
        d_ctrs = (d_mins[:, p_idx, :] + d_maxs[:, p_idx, :]) / 2
        plot_with_sd(ax_d, Ks, d_ctrs, axis=1, label=label, **kw)
        plot_with_sd(ax_a, Ks, a_hats[:, p_idx, :], axis=1, label=label, **kw)
        for ax in [ax_d, ax_a]:
            ax.set_xscale('log')
            ax.set_xticks(Ks)
            ax.set_xticklabels(Ks, rotation=90, fontsize=6)
            ax.set_xlabel('similarity length')
        ax_a.set_yscale('log')

    ax_d.set_ylabel('estimated diagonal position')
    ax_d.legend(loc='best', fontsize=8)

    ax_a.plot(Ks, [2 * K for K in Ks], **truth_kw)
    ax_a.set_ylabel('estimated antidiagonal position')
    ax_a.legend(loc='best', fontsize=8)

    savefig(fig_hats, 'simulations[estimates]%s.png' % suffix)


def exp_simulated_K_p():
    """Performance of Word-Blot and associated statistical scores for pairwise
    local similarity search in *simulations*:

    * z-scores assigned to diagonal strips with and without
      similarities under the corresponding limiting Normal distributions
      (:math:`H0, H1`), for varying similarity lengths and match probabilities.
    * Estimated coordinates, lengths and match probabilities of local
      similarities for varying similarity lengths and match probabilities. In
      each trial with similarity length :math:`K` and match probability
      :math:`p`, two pairs of input sequences are provided to Word-Blot:
      a pair of sequences of length :math:`2K` whose substrings at positions
      :math:`[\\frac{K}{2}, \\frac{3K}{2}]` are homologies of length :math:`K`
      and match probability :math:`p`., and two unrelated sequences of length
      :math:`2K`. are considered.

    **Supported Claims**

    * z-scores calculated by :func:`biseqt.blot.WordBlot.score_num_seeds`
      against the limiting normal distributions :math:`H_0, H_1` (unrelated and
      related) are stable and comparable across different similarity lengths
      and match probabilities; thus both scores are reliable statistics.
    * Local similarities found by Word-Blot as per
      :func:`biseqt.blot.WordBlot.similar_segments` are accurate in length,
      estimated match probability, and coordinates.

    .. figure::
        https://www.dropbox.com/s/s38hi2oo9pp78ul/
        simulations%5Bscores%5D.png?raw=1
       :target:
        https://www.dropbox.com/s/s38hi2oo9pp78ul/
        simulations%5Bscores%5D.png?raw=1
       :alt: lightbox

       Z-scores of the number of seeds in diagonal strip for related (solid
       lines) and unrelated (dashed lines) segments of varying lengths,
       calculated against the limiting distribution for unrelated pairs of
       sequences (*left*) and related pairs (*right*) as a function of
       similarity length (*top*) and match probability (*bottom*), n=50
       samples, shaded regions indicate one standard deviation. Note that, as
       desired, H0 score is stable for unrelated sequences of any length and H1
       score is stable for related sequences of any length or match
       probability.

    .. figure::
        https://www.dropbox.com/s/lqm4s4evmbkk4as/
        simulations%5Bestimates%5D.png?raw=1
       :target:
        https://www.dropbox.com/s/lqm4s4evmbkk4as/
        simulations%5Bestimates%5D.png?raw=1
       :alt: lightbox

       Estimated length (*top left*), match probability (*bottom left*),
       diagonal position (*top right*) and antidiagonal position (*bottom
       right*) of local similarities, n=50 samples, shaded regions indicate one
       standard deviation. Dashed black lines are ground truth for comparison.
       Length estimates are the sum of the lengths of all reported similar
       segments. Diagonal and antidiagonal coordinates are averaged over all
       reported segments. Match probability is taken from the longest similar
       segment.
    """
    Ks = [100 * 2 ** i for i in range(1, 8)]
    select_Ks = Ks[1], Ks[3], Ks[5]

    ps = [round(1 - .06 * i, 2) for i in range(1, 8)]
    select_ps = ps[0], ps[3], ps[5]

    n_samples = 50
    wordlen = 6

    suffix = ''
    dumpfile = 'simulations%s.txt' % suffix
    sim_data = sim_simulated_K_p(
        Ks, ps, n_samples, wordlen=wordlen,
        dumpfile=dumpfile, ignore_existing=False)
    plot_simulated_K_p(sim_data, select_Ks, select_ps, suffix=suffix)


@with_dumpfile
def sim_comp_aligned_genes(seqs, pws, **kw):
    A = Alphabet('ACGT')
    wordlen, K_min, p_min = kw['wordlen'], kw['K_min'], kw['p_min']
    WB_kw = {
        'g_max': kw.get('g_max', .3),
        'sensitivity': kw.get('sensitivity', .99),
        'wordlen': wordlen,
        'alphabet': A,
        'path': kw.get('db_path', ':memory:'),
        'log_level': kw.get('log_level', logging.WARNING),
    }
    qM = 1 / p_min - 1
    qS = qG = -1

    aligner_kw = {
        'match_score': qM,
        'mismatch_score': qS,
        'ge_score': qG,
        'go_score': 0,
        'alnmode': BANDED_MODE,
        'alntype': B_GLOBAL,
    }

    similar_segments_kw = {'K_min': K_min, 'p_min': p_min}
    sim_data = {
        'pws': pws,
        'seqlens': {name: len(seqs[name]) for name in seqs},
        'seed_pos': [],
        'seed_neg': [],
        'seeds': {key: {} for key in pws},
        'opseq_seeds': {key: [] for key in pws},
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
        sim_data['opseq_seeds'][(id1, id2)] = opseq_seeds
        for _, seed in opseq_seeds:
            d, a = WB.to_diagonal_coordinates(*seed)
            for d_ in range(d - band_r, d + band_r):
                for a_ in range(a - band_r, a + band_r):
                    i_, j_ = WB.to_ij_coordinates(d_, a_)
                    if 0 <= i_ < len(S) and 0 <= j_ < len(T):
                        in_band[i_, j_] = 1

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

        # run DP alignment on found segments
        for idx, seg_info in enumerate(sim_data['segments'][(id1, id2)]):
            seg = seg_info['segment']
            (i_start, i_end), (j_start, j_end) = WB.to_ij_coordinates_seg(seg)
            # to_ij_coordinates_seg might overflow; fix it
            i_end = min(sim_data['seqlens'][id1] - 1, i_end)
            j_end = min(sim_data['seqlens'][id2] - 1, j_end)
            # allow longer alignments to be found
            start_shift = min(i_start, j_start, 2 * wordlen)
            end_shift = min(len(seqs[id1]) - i_end,
                            len(seqs[id2]) - j_end, 2 * wordlen)
            i_start, j_start = i_start - start_shift, j_start - start_shift
            i_end, j_end = i_end + end_shift, j_end + end_shift
            S_ = S[i_start: i_end]
            T_ = T[j_start: j_end]
            rad = (seg[0][1] - seg[0][0]) / 2
            rad = min(len(S_), len(T_), max(rad, abs(len(S_) - len(T_)) + 2))
            aligner_kw['diag_range'] = (-rad, rad)
            aligner = Aligner(S_, T_, **aligner_kw)
            with aligner:
                aligner.solve()
                alignment = aligner.traceback()
                len_orig = alignment.projected_len(alignment.transcript)
                alignment = alignment.truncate_to_match()
                sim_data['segments'][(id1, id2)][idx]['alignment'] = alignment
                sim_data['segments'][(id1, id2)][idx]['frame'] = \
                    (i_start, i_end), (j_start, j_end)
                if alignment is not None:
                    tx = alignment.transcript
                    proj_len = alignment.projected_len(tx)
                    print round(seg_info['p'], 3), \
                        round(1. * tx.count('M') / proj_len, 3), \
                        proj_len, len_orig
    return sim_data


def plot_comp_aligned_genes(sim_data, suffix='', example_pair=None,
                            naming_style=None):
    seqlens = sim_data['seqlens']
    pws = sim_data['pws']
    K_min = sim_data['similar_segments_kw']['K_min']
    p_min = sim_data['similar_segments_kw']['p_min']

    # dims of plot containing subplots for all pairs of sequences; these plots
    # are there just for posterity
    fig_num = int(np.ceil(np.sqrt(len(pws))))
    fig_seeds = plt.figure(figsize=(6 * fig_num, 5 * fig_num))
    fig_profiles = plt.figure(figsize=(8 * fig_num, 4 * fig_num))

    fig_main = plt.figure(figsize=(13, 7))
    grids = gridspec.GridSpec(2, 6, height_ratios=[3, 1.3])

    ax_seeds_ex = fig_main.add_subplot(grids[0, :2])
    ax_p = fig_main.add_subplot(grids[0, 2:4])
    ax_p_aln = fig_main.add_subplot(grids[0, 4:])
    ax_profile_ex = fig_main.add_subplot(grids[1, :3])
    ax_p_cdf = fig_main.add_subplot(grids[1, 3:])

    labels = []
    p_alns = []
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

        ps_true = estimate_match_probs_in_opseq(opseq, radius)
        ps_true_smooth = gaussian_filter1d(ps_true, 10)
        ps_hat_pos, ps_hat = [], []
        for pos, seed in sim_data['opseq_seeds'][key]:
            if ps_true[pos] == 0:
                # flanking zero regions in estimated gap probabilities
                continue
            ps_hat.append(sim_data['seeds'][key][seed])
            ps_hat_pos.append(pos)

        ax_profile = fig_profiles.add_subplot(fig_num, fig_num, idx + 1)
        axes_profile = [ax_profile]
        if key == example_pair:
            axes_profile.append(ax_profile_ex)

        for ax in axes_profile:
            ax.plot(range(len(ps_true)), ps_true_smooth, c='g', alpha=.8, lw=1)
            ax.set_title('similarity profile (%s vs. %s)' % (id1, id2))
            ax.set_xlabel('position along global alignment')
            ax.set_ylabel('match probability')
            ax.scatter(ps_hat_pos, ps_hat, lw=0, c='k', s=7, alpha=.4)

        # ps_true is calculated from the alignment and is against alignment
        # length; whereas our p_hat is estimated againsted the projected
        # length. Correct ps_true to be against projected length too.
        gs_true = estimate_gap_probs_in_opseq(opseq, radius)
        ps_true_proj = [
            (1 / (1 - gs_true[ps_hat_pos][i])) * ps_true[ps_hat_pos[i]]
            for i in range(len(ps_hat))]
        ax_p.scatter(ps_hat, ps_true_proj, lw=0, color='g', s=3, alpha=.5)

        # =============
        # Dot Plots
        # =============
        ax_seeds = fig_seeds.add_subplot(fig_num, fig_num, idx + 1)
        axes_seeds = [ax_seeds]
        if key == example_pair:
            axes_seeds.append(ax_seeds_ex)
        for ax in axes_seeds:
            plot_scored_seeds(ax, sim_data['seeds'][key].items(),
                              threshold=p_min, alpha=.7, zorder=9)
            plot_global_alignment(ax, opseq, lw=5, alpha=.2, color='g')
            adjust_pw_plot(ax, seqlens[key[0]], seqlens[key[1]])
            ax.set_ylabel(id1)
            ax.set_xlabel(id2)

        # ================================
        # match probability estimation vs local alignments
        # ================================
        for seg_info in sim_data['segments'][key]:
            p_hat = seg_info['p']
            alignment = seg_info['alignment']
            if alignment is None:
                p_aln = 0
            else:
                transcript = alignment.transcript
                K_aln = alignment.projected_len(transcript)
                p_aln = 1. * transcript.count('M') / K_aln
                for ax in axes_seeds:
                    plot_local_alignment(ax, transcript,
                                         seg_info['frame'][0][0],
                                         seg_info['frame'][1][0],
                                         lw=1, alpha=.9, color='b', zorder=10)
            p_alns.append(p_aln)
            ax_p_aln.scatter([p_hat], [p_aln], lw=0, color='b', s=10, alpha=.3)

    for ax in [ax_p, ax_p_aln]:
        ax.plot([0, 1], [0, 1], lw=2, ls='--', alpha=.6, c='k')
        ax.set_xlabel('estimated match probability')
        ax.set_ylabel('true match probability')
        ax.set_xlim(-.1, 1.1)
        ax.set_ylim(-.1, 1.1)
        ax.set_aspect('equal')
    ax_p.set_title('homologous seeds')
    ax_p_aln.set_title('similar segments')
    ax_p_aln.plot([0, 1], [p_min, p_min], lw=2, ls='--', alpha=.6, c='k')

    plot_cdf(ax_p_cdf, p_alns, c='b', lw=1, alpha=.9, smooth_radius=.5)
    ax_p_cdf.axvline(x=p_min, ls='--', lw=2, c='k', alpha=.6)
    ax_p_cdf.set_xlabel('local alignment match probability')
    ax_p_cdf.set_ylabel('cumulative distribution')
    ax_p_cdf.set_xlim(0.4, 1)

    for ax in [ax_p_cdf, ax_profile_ex]:
        ticks = [i * .1 for i in range(11)]
        ax.set_yticks(ticks)
        ax.set_yticklabels(ticks, fontsize=6)
    ax_profile_ex.axhline(y=p_min, ls='--', lw=2, c='k', alpha=.6)

    for fig in [fig_main, fig_seeds, fig_profiles]:
        fig.tight_layout()
    savefig(fig_main, 'comp_aligned_genes[estimates]%s.png' % suffix)
    savefig(fig_seeds, 'comp_aligned_genes[seeds]%s.png' % suffix)
    savefig(fig_profiles, 'comp_aligned_genes[profiles]%s.png' % suffix)

    fig, _ = plot_classifier(sim_data['seed_pos'], sim_data['seed_neg'],
                             labels=['homologous', 'non-homologous'],
                             mark_threshold=.8)
    savefig(fig, 'comp_aligned_genes[classifier]%s.png' % suffix)


def exp_comp_aligned_genes():
    """Performance of Word-Blot and associated statistical scores for pairwise
    local similarity search on *biological data*. Given aligned copies of a
    gene in multiple species we consider the following for each pair:

    * estimated match probability :math:`\hat{p}` at *homologous seeds* (i.e.
      those exactly matching kmers that lie on the known global pairwise
      alignment) compared to local match probability of said global alignment.
    * discrimination power of :math:`\hat{p}` to separate homologous seeds from
      non-homologous seeds.
    * all identified similar segments between the two sequences verified by
      comparison to banded dynamic programming alignment in the identified
      diagonal strip.

    **Supported Claims**

    * Word-Blot accurately estimates match probabilities, length, and
      coordinates of local similarities in biological data.
    * estimated match probabilities :math:`\hat{p}` are effective statistics to
      separate homologous and non-homologous seeds.
    * Dynamic programming alignment confirms that Word-Blot correctly
      identifies local similarities, even those that are inevitably missed by
      global alignments.

    .. figure::
        https://www.dropbox.com/s/8ajed2uo0qobyjw/
        comp_aligned_genes%5Bestimates%5D%5Birx1%5D%5Bp%3D0.80%5D.png?raw=1
       :target:
        https://www.dropbox.com/s/8ajed2uo0qobyjw/
        comp_aligned_genes%5Bestimates%5D%5Birx1%5D%5Bp%3D0.80%5D.png?raw=1
       :alt: lightbox

       Word-Blot performance in identifying local similarities among *Iroquois
       homeobox protein 1 (IRX1)* genes for 10 amniota obtain from Ensembl with
       word length 8, minimum match probability 0.8, and minimum similarity
       length 200nt. An example pairwise comparison (*top left*) shows seeds
       (grey) color coded by intensity according to their estimated match
       probability with those exceeding the threshold shown larger for clarity.
       Known global alignment is shown in green and aligned local similarities
       identified by Word-Blot are shown in blue. Estimated similarity profile
       for the same example pair (*bottom left*) is shown at positions
       exceeding the match probability minimum (black dots) and compared to the
       profile obtained from the known global alignment (green). For all
       pairwise comparisons estimated and true match probability at homologous
       seeds are shown (*top middle*, green), the diagonal shows ground truth
       for comparison (dashed black).  Similarly for all identified similar
       segments, estimated match probability and true probability obtained from
       banded global alignment are shown (*top right*, blue), the diagonal
       showing ground truth (dashed black) and the horizontal line shows the
       cutoff used by Word-Blot (dashed black).

    .. figure::
        https://www.dropbox.com/s/2mc1rlnazodeu6l/
        comp_aligned_genes%5Bclassifier%5D%5Birx1%5D%5Bp%3D0.80%5D.png?raw=1
       :target:
        https://www.dropbox.com/s/2mc1rlnazodeu6l/
        comp_aligned_genes%5Bclassifier%5D%5Birx1%5D%5Bp%3D0.80%5D.png?raw=1
       :alt: lightbox

       Performance of Word-Blot estimated match probability at discriminating
       between homologous and non-homologous seeds in the same context as the
       figure above. The cumulative distribution of estimated match probability
       in each case (*bottom right*), the resulting ROC curve (*left*) and the
       predictive ROC curve (*top right*). The threshold 0.8 is indicated in
       blue in all three plots.


    .. figure::
        https://www.dropbox.com/s/38xgx8gbvziq126/
        comp_aligned_genes%5Bestimates%5D%5Bactb%5D%5Bp%3D0.50%5D.png?raw=1
       :target:
        https://www.dropbox.com/s/38xgx8gbvziq126/
        comp_aligned_genes%5Bestimates%5D%5Bactb%5D%5Bp%3D0.50%5D.png?raw=1
       :alt: lightbox

       Word-Blot performance in identifying local similarities among *Beta
       Actin (ACTB)* genes for 7 vertebrates obtain from UCSC genome browser
       with word length 6, minimum match probability 0.5, and minimum
       similarity length 200nt. All subplots are similar to those of *IRX1*
       experiment shown above.

    .. figure::
        https://www.dropbox.com/s/t5pff0nk54m8swr/
        comp_aligned_genes%5Bclassifier%5D%5Bactb%5D%5Bp%3D0.50%5D.png?raw=1
       :target:
        https://www.dropbox.com/s/t5pff0nk54m8swr/
        comp_aligned_genes%5Bclassifier%5D%5Bactb%5D%5Bp%3D0.50%5D.png?raw=1
       :alt: lightbox

       Performance of Word-Blot estimated match probability at discriminating
       between homologous and non-homologous seeds in the same context as the
       figure above.
    """
    A = Alphabet('ACGT')

    # p_min = .5
    # K_min = 200
    # wordlen = 6
    # style = 'ucsc'
    # suffix = '[actb][p=%.2f]' % p_min
    # dumpfile = 'comp_aligned_genes%s.txt' % suffix
    # seqs_path = 'data/actb/actb-7vet.fa'
    # pws_path = 'data/actb/actb-7vet-pws.fa'
    # example_pair = ('hg38.chr7', 'mm10.chr5')

    p_min = .8
    wordlen = 8
    K_min = 200
    style = 'ensembl'
    suffix = '[irx1][p=%.2f]' % p_min
    dumpfile = 'comp_aligned_genes%s.txt' % suffix
    seqs_path = 'data/irx1/irx1-vert-amniota-indiv.fa'
    pws_path = 'data/irx1/irx1-vert-amniota-pws.fa'
    example_pair = ('homo_sapiens/1-10720', 'bos_taurus/1-10720')

    with open(seqs_path) as f:
        seqs = {name: A.parse(fill_in_unknown(seq.upper(), A))
                for seq, name, _ in load_fasta(f)}
    with open(pws_path) as f:
        pws = {tuple(name.split(':')): seq for seq, name, _ in load_fasta(f)}
        pws = {key: value for key, value in pws.items()}
    sim_data = sim_comp_aligned_genes(
        seqs, pws, wordlen=wordlen, p_min=p_min, K_min=K_min,
        dumpfile=dumpfile, ignore_existing=False)
    plot_comp_aligned_genes(sim_data, suffix=suffix, example_pair=example_pair,
                            naming_style=style)


if __name__ == '__main__':
    exp_simulated_K_p()
    exp_comp_aligned_genes()
