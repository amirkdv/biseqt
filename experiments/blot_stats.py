import numpy as np
from matplotlib import pyplot as plt
from itertools import product
from scipy.ndimage.filters import gaussian_filter1d

import sys
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
from util import plot_cdf
from util import estimate_gap_probs_in_opseq
from biseqt.pw import Aligner, BANDED_MODE, B_GLOBAL


@with_dumpfile
def sim_simulated_K_p(Ks, ps, n_samples, **kw):
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
        'a_hat': _zero(),
        'd_min': _zero(),
        'd_max': _zero(),
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

                # sum of K_hat, average of a_hat, d_min, d_max
                sim_data['K_hat'][key][K_idx, p_idx, idx] = sum(
                    r['segment'][1][1] - r['segment'][1][0]
                    for r in results) / 2
                sim_data['a_hat'][key][K_idx, p_idx, idx] = sum(
                    (r['segment'][1][1] + r['segment'][1][0]) / 2
                    for r in results) / len(results)
                sim_data['d_min'][key][K_idx, p_idx, idx] = sum(
                    r['segment'][0][0] for r in results) / len(results)
                sim_data['d_max'][key][K_idx, p_idx, idx] = sum(
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
    wordlen = sim_data['WB_kw']['wordlen']
    # n_samples = sim_data['scores']['H0']['pos'].shape[2]
    scores = sim_data['scores']

    kw = {'marker': 'o', 'markersize': 3, 'alpha': .7, 'lw': 1}
    truth_kw = {'ls': '--', 'alpha': .6, 'color': 'k', 'lw': 1}

    # ======================================
    # score by varying K for select ps
    # score varying p for select Ks
    fig_scores = plt.figure(figsize=(8, 8))
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
        ax.set_ylabel('%s score' % mode, fontsize=10)
        ax.set_xlabel('similarity length', fontsize=10)
        ax.set_xscale('log')
        ax.set_xticks(Ks)
        ax.set_xticklabels(Ks, fontsize=10)
        ax.legend(loc='best', fontsize=10)

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
        ax.legend(loc='best', fontsize=10)
    savefig(fig_scores, 'simulations[scores]%s.png' % suffix)

    # =====================================
    fig_hats = plt.figure(figsize=(8, 8))
    ax_K = fig_hats.add_subplot(2, 2, 1)
    ax_p = fig_hats.add_subplot(2, 2, 3)
    ax_d = fig_hats.add_subplot(2, 2, 2)
    ax_a = fig_hats.add_subplot(2, 2, 4)

    K_hats = sim_data['K_hat']
    p_hats = sim_data['p_hat']

    # estimated Ks for select ps
    colors = color_code(select_ps)
    for p, color in zip(select_ps, colors):
        p_idx = ps.index(p)
        plot_with_sd(ax_K, Ks, K_hats['pos'][:, p_idx, :], axis=1,
                     color=color, label='p = %.2f' % p, **kw)
        # plot_with_sd(ax_K, Ks,
                     # np.maximum(K_hats['neg'][:, p_idx, :], 2 * wordlen),
                     # axis=1, color=color, ls='--', **kw)
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
        plot_with_sd(ax_p, ps, p_hats['pos'][K_idx, :, :], axis=1,
                     color=color, label='K = %d' % K, **kw)
        # plot_with_sd(ax_p, ps, p_hats['neg'][K_idx, :, :], axis=1,
                     # color=color, ls='--', **kw)
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
    d_mins = sim_data['d_min']['pos']
    d_maxs = sim_data['d_max']['pos']
    a_hats = sim_data['a_hat']['pos']

    g_max = sim_data['WB_kw']['g_max']
    sensitivity = sim_data['WB_kw']['sensitivity']

    colors = color_code(select_ps)
    for p, color in zip(select_ps, colors):
        label = 'p = %.2f' % p
        kw = {'color': color, 'alpha': .6, 'lw': 1,
              'marker': 'o', 'markersize': 3}
        p_idx = ps.index(p)

        plot_with_sd(ax_d, Ks, d_mins[:, p_idx, :], axis=1, label=label, **kw)
        plot_with_sd(ax_d, Ks, d_maxs[:, p_idx, :], axis=1, **kw)
        plot_with_sd(ax_a, Ks, a_hats[:, p_idx, :], axis=1, label=label, **kw)
        ax_d.plot(Ks, band_radii(Ks, g_max, sensitivity), **truth_kw)
        ax_d.plot(Ks, -band_radii(Ks, g_max, sensitivity), **truth_kw)
        for ax in [ax_d, ax_a]:
            ax.set_xscale('log')
            ax.set_xticks(Ks)
            ax.set_xticklabels(Ks, rotation=90, fontsize=6)
            ax.set_xlabel('similarity length')
        ax_a.set_yscale('log')

    ax_d.set_ylabel('estimated diagonal position of similarity')
    ax_d.legend(loc='best', fontsize=8)

    ax_a.plot(Ks, [2 * K for K in Ks], **truth_kw)
    ax_a.set_ylabel('estimated antidiagonal position of similarity')
    ax_a.legend(loc='best', fontsize=8)

    savefig(fig_hats, 'simulations[estimates]%s.png' % suffix)


def exp_simulated_K_p():
    """Performance of Word-Blot and associated statistical scores for pairwise
    local similarity search:

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

    # qM, qS, qG = 1, 0, -1
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
        sim_data['opseq_seeds'][(id1, id2)] = opseq_seeds
        for _, seed in opseq_seeds:
            d, a = WB.to_diagonal_coordinates(*seed)
            for d_ in range(d - band_r, d + band_r):
                # FIXME what should the range be?
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
        # FIXME report a global (over all pairs) TPR too
        log('tpr = %.2f' % sim_data['homologous_tpr'][(id1, id2)])

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


def plot_comp_aligned_genes(sim_data, suffix='', naming_style=None):
    seqlens = sim_data['seqlens']
    pws = sim_data['pws']
    K_min = sim_data['similar_segments_kw']['K_min']
    p_min = sim_data['similar_segments_kw']['p_min']
    fig_num = int(np.ceil(np.sqrt(len(pws))))

    fig_seeds = plt.figure(figsize=(6 * fig_num, 5 * fig_num))
    fig_profiles = plt.figure(figsize=(8 * fig_num, 4 * fig_num))
    fig_p = plt.figure(figsize=(6, 4))
    fig_p_aln = plt.figure(figsize=(10, 4))
    fig_K_aln = plt.figure(figsize=(6, 4))
    fig_coord_classifier = plt.figure(figsize=(6, 4))

    ax_coord_classifier = fig_coord_classifier.add_subplot(1, 1, 1)
    ax_p = fig_p.add_subplot(1, 1, 1)
    ax_p_aln = fig_p_aln.add_subplot(1, 2, 1)
    ax_p_aln_cdf = fig_p_aln.add_subplot(1, 2, 2)
    ax_K_aln = fig_K_aln.add_subplot(1, 1, 1)

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
            if ps_true[pos] == 0:
                # flanking zero regions in estimated gap probabilities
                continue
            ps_hat.append(sim_data['seeds'][key][seed])
            ps_hat_pos.append(pos)
        ax_profiles.scatter(ps_hat_pos, ps_hat, lw=0, c='k', s=4, alpha=.4)
        # ps_true is calculated from the alignment and is against alignment
        # length; whereas our p_hat is estimated againsted the projected
        # length. Correct ps_true to be against projected length too.
        gs_true = estimate_gap_probs_in_opseq(opseq, radius)
        ps_true_proj = [
            (1 / (1 - gs_true[ps_hat_pos][i])) * ps_true[ps_hat_pos[i]]
            for i in range(len(ps_hat))]
        ax_p.scatter(ps_hat, ps_true_proj, lw=0, color='g', s=4, alpha=.2)

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
            K_hat = (seg_info['segment'][1][1] - seg_info['segment'][1][0]) / 2
            if alignment is None:
                p_aln = 0
                K_aln = 0
            else:
                transcript = alignment.transcript
                K_aln = alignment.projected_len(transcript)
                p_aln = 1. * transcript.count('M') / K_aln
                plot_local_alignment(ax_seeds, transcript,
                                     seg_info['frame'][0][0],
                                     seg_info['frame'][1][0],
                                     lw=7, alpha=.2, color='b')
            p_alns.append(p_aln)
            ax_p_aln.scatter([p_hat], [p_aln], lw=0, color='k', s=4, alpha=.6)
            ax_K_aln.scatter([K_hat], [K_aln], lw=0, color='k', s=10, alpha=.6)

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

    plot_cdf(ax_p_aln_cdf, p_alns, c='k', lw=1, alpha=.7, smooth_radius=.5)
    ax_p_aln_cdf.axvline(x=p_min, c='k', alpha=.3, lw=3)
    ax_p_aln_cdf.set_xlabel('local alignment match probability')
    ax_p_aln_cdf.set_xlim(0.4, 1.1)
    ax_p_aln.plot([0, 1], [0, 1], lw=1, ls='--', alpha=.6, c='k')
    ax_p_aln.plot([0, 1], [p_min, p_min], lw=1, ls='--', alpha=.6, c='k')
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

    fig, _ = plot_classifier(sim_data['seed_pos'], sim_data['seed_neg'],
                             labels=['homologous', 'non-homologous'],
                             mark_threshold=.8)
    savefig(fig, 'comp_aligned_genes[seed-classifier]%s.png' % suffix)

    for fig in [fig_profiles, fig_p, fig_p_aln, fig_K_aln, fig_seeds,
                fig_coord_classifier]:
        fig.tight_layout()
    savefig(fig_profiles, 'comp_aligned_genes[profiles]%s.png' % suffix)
    savefig(fig_p, 'comp_aligned_genes[p-hat]%s.png' % suffix)
    savefig(fig_p_aln, 'comp_aligned_genes[p-hat-aln]%s.png' % suffix)
    savefig(fig_K_aln, 'comp_aligned_genes[K-hat-aln]%s.png' % suffix)
    savefig(fig_seeds, 'comp_aligned_genes[seeds]%s.png' % suffix)
    savefig(fig_coord_classifier,
            'comp_aligned_genes[coords-classifier]%s.png' % suffix)


def exp_comp_aligned_genes():
    A = Alphabet('ACGT')

    # FIXME seed plots seem weird (seeds above threshold without segments).
    # FIXME the homologous classifier ROC needs work to justify (whasn't this
    # the point of including irx1?). Also, I think we can get rid of TPR plot.
    # p_min = .5
    # K_min = 200
    # wordlen = 6
    # style = 'ucsc'
    # suffix = '[actb][p=%.2f]' % p_min
    # dumpfile = 'comp_aligned_genes%s.txt' % suffix
    # seqs_path = 'data/actb/actb-7vet.fa'
    # pws_path = 'data/actb/actb-7vet-pws.fa'

    p_min = .8
    wordlen = 8
    K_min = 100
    style = 'ensembl'
    suffix = '[irx1][p=%.2f]' % p_min
    dumpfile = 'comp_aligned_genes%s.txt' % suffix
    seqs_path = 'data/irx1/irx1-vert-amniota-indiv.fa'
    pws_path = 'data/irx1/irx1-vert-amniota-pws.fa'

    with open(seqs_path) as f:
        seqs = {name: A.parse(fill_in_unknown(seq.upper(), A))
                for seq, name, _ in load_fasta(f)}
    with open(pws_path) as f:
        pws = {tuple(name.split(':')): seq for seq, name, _ in load_fasta(f)}
        pws = {key: value for key, value in pws.items()}
    sim_data = sim_comp_aligned_genes(
        seqs, pws, wordlen=wordlen, p_min=p_min, K_min=K_min,
        dumpfile=dumpfile, ignore_existing=False)
    plot_comp_aligned_genes(sim_data, suffix=suffix, naming_style=style)


if __name__ == '__main__':
    exp_simulated_K_p()
    exp_comp_aligned_genes()
