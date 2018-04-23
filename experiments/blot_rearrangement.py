import numpy as np
from matplotlib import gridspec
import matplotlib
import logging
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from biseqt.pw import Aligner, STD_MODE, LOCAL
from biseqt.blot import WordBlot
from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess
from util import plot_scored_seeds, plot_seeds
from util import plot_similar_segment, adjust_pw_plot
from util import log, savefig


def exp_rearrangement():
    """Example demonstrating of Word-Blot for pairwise local similarity search on
    two randomly generated sequencees with motif sequences violating
    collinearity :math:`S=M_1M_2M_3, T=M'_1M'_1M'_3M'_2` where motif pairs
    :math:`(M_i, M'_i)_{i=1,2,3}` have lengths 200, 400, 600 and are related by
    match probabilities 0.95, 0.85, and 0.75, respectively.

    .. figure::
        https://www.dropbox.com/s/nsvsf5gaui6t9ww/rearrangement.png?raw=1
       :target:
        https://www.dropbox.com/s/nsvsf5gaui6t9ww/rearrangement.png?raw=1
       :alt: lightbox

       Dynamic programming scores of the forward pass of Smith Waterman are
       shown in color code (*left*) with seeds (word length 6) grey intensity
       coded according to the local match probability assigned by Word-Blot
       (minimum similarity length 200). Similar segments reported by Word-Blot
       are shown as grey diagonal strips (*left*) and schematically (*right*)
       color coded by their Word-Blot estimated match probabilities (note
       agreement with true match probabilities).
    """
    # NOTE we are running whole table DP later here; be careful with size
    K = 200
    wordlen = 6
    A = Alphabet('ACGT')

    WB_kw = {'g_max': .2, 'sensitivity': .9, 'alphabet': A, 'wordlen': wordlen,
             'path': ':memory:', 'log_level': logging.INFO}

    # homologies
    Hs = [rand_seq(A, i) for i in [i * K for i in range(1, 4)]]
    ps = [.95, .85, .75]
    Ms = []
    for p_match in ps:
        subst = gap = 1 - np.sqrt(p_match)
        print subst, gap
        Ms.append(
            MutationProcess(A, subst_probs=subst, ge_prob=gap, go_prob=gap)
        )

    # connector junk
    def J(): return rand_seq(A, 2 * K)

    S = J() + Hs[0] + J() + Hs[1] + J() + Hs[2] + J()
    Hs = [M.mutate(hom)[0] for hom, M in zip(Hs, Ms)]
    T = J() + Hs[0] + J() + Hs[0] + Hs[2] + J() + Hs[1] + J()

    fig = plt.figure(figsize=(9, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
    ax_seeds = plt.subplot(gs[0])
    ax_mapping = plt.subplot(gs[1])

    WB = WordBlot(S, T, **WB_kw)

    p_min = .95 * min(ps)
    scored_seeds = WB.score_seeds(K)
    scored_seeds = [(WB.to_ij_coordinates(*rec['seed']), rec['p'])
                    for rec in scored_seeds]
    plot_seeds(ax_seeds, [x[0] for x in scored_seeds])

    cmap = plt.cm.get_cmap('plasma')
    sim_segments = list(WB.similar_segments(K_min=K, p_min=p_min))
    min_p_obs = min(rec['p'] for rec in sim_segments)
    max_p_obs = max(rec['p'] for rec in sim_segments)
    for rec in sim_segments:
        print rec
        seg = rec['segment']
        (i_start, i_end), (j_start, j_end) = WB.to_ij_coordinates_seg(seg)
        i_ctr, j_ctr = (i_start + i_end) / 2, (j_start + j_end) / 2
        color = cmap((rec['p'] - min_p_obs) / (max_p_obs - min_p_obs))[:3]
        plot_similar_segment(ax_seeds, seg, lw=5, alpha=.1, c='k')
        ax_mapping.plot([1, 1], [i_start, i_end], lw=3, c=color, alpha=.7)
        ax_mapping.plot([2, 2], [j_start, j_end], lw=3, c=color, alpha=.7)
        ax_mapping.plot([1, 2], [i_ctr, j_ctr], marker='o', markersize=7, lw=2,
                        c=color, alpha=.4)

    ax_mapping.set_xticks([1, 2])
    ax_mapping.set_xticklabels(['S', 'T'])
    ax_mapping.set_xlim(0, 3)
    ax_mapping.set_ylim(0, None)
    ax_c = make_axes_locatable(ax_mapping).append_axes('right', size='4%',
                                                       pad=0.05)
    norm = matplotlib.colors.Normalize(vmin=min_p_obs, vmax=max_p_obs)
    matplotlib.colorbar.ColorbarBase(ax_c, cmap=cmap, norm=norm,
                                     orientation='vertical')

    aligner_kw = {
        'match_score': 1 / p_min - 1,
        'mismatch_score': -1,
        'ge_score': -1,
        'go_score': 0,
        'alnmode': STD_MODE,
        'alntype': LOCAL,
    }

    print len(S), len(T)
    with Aligner(S, T, **aligner_kw) as aligner:
        aligner.solve()
        scores = np.array(aligner.table_scores())
        min_score = min(scores.flatten())
        max_score = max(scores.flatten())
        ax_seeds.imshow(scores, cmap='plasma', alpha=.3)
        ax_c = make_axes_locatable(ax_seeds).append_axes('right', size='4%',
                                                         pad=0.05)
        norm = matplotlib.colors.Normalize(vmin=min_score, vmax=max_score)
        matplotlib.colorbar.ColorbarBase(ax_c, cmap='plasma', norm=norm,
                                         orientation='vertical')

    adjust_pw_plot(ax_seeds, len(S), len(T))
    ax_seeds.set_xlabel('T')
    ax_seeds.set_ylabel('S')

    fig.tight_layout()
    savefig(fig, 'rearrangement.png')


def exp_repeat_regions():
    gap = .2
    subst = .1
    K = 500
    wordlen = 8
    A = Alphabet('ACGT')

    # NOTE I can drive sensitivity to 0 and get decent results
    WB_kw = {'g_max': .2, 'sensitivity': .9, 'alphabet': A, 'wordlen': wordlen,
             'path': ':memory:', 'log_level': logging.INFO}

    M = MutationProcess(A, subst_probs=subst, ge_prob=gap, go_prob=gap)

    homs = [rand_seq(A, i) for i in [K/2, K, 2 * K, 4 * K]]

    def junk(): return rand_seq(A, np.random.randint(2 * K, 4 * K))

    junks = [junk() for _ in range(3 * len(homs))]
    S = sum([junks[3 * i] + R + junks[3 * i + 1] + R + junks[3 * i + 2] + R
             for i, R in enumerate(homs)], A.parse('')) + junk()
    homs = [M.mutate(homs[i])[0] for i in range(len(homs))]
    T = S

    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(1, 1, 1)

    log('finding repeat regions')
    WB = WordBlot(S, T, **WB_kw)
    match = (1 - gap) * (1 - subst)

    scored_seeds = WB.score_seeds(K)
    # convert to ij coordinates and exclude half of the table (mirror image)
    scored_seeds = [(WB.to_ij_coordinates(*rec['seed']), rec['p'])
                    for rec in scored_seeds if rec['seed'][0] <= 0]

    plot_scored_seeds(ax, scored_seeds)
    for rec in WB.similar_segments(K_min=K, p_min=match):
        segment, scores, p_hat = rec['segment'], rec['scores'], rec['p']
        (d_min, d_max), (a_min, a_max) = rec['segment']
        if d_min > 0:
            # self comparison, exclude half of the table
            continue
        log('repeat region %s, scores = (%.2f, %.2f), p = %.2f' %
            (str(segment), scores[0], scores[1], p_hat))
        plot_similar_segment(ax, segment, c='b', lw=5, alpha=.2)
    adjust_pw_plot(ax, len(S), len(T))
    ax.set_title('Repeat regions', y=1.05, fontsize=10)

    fig.suptitle('word len. = %d, min. hom. len = %d, min.match = %.2f' %
                 (wordlen, K, match), fontsize=8)
    fig.tight_layout()
    savefig(fig, 'repeat_regions.png')


if __name__ == '__main__':
    exp_rearrangement()
    exp_repeat_regions()
