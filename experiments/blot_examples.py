import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib
import logging
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from biseqt.blot import WordBlot
from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess
from util import plot_seeds, plot_scored_seeds
from util import plot_similar_segment, adjust_pw_plot
from util import log, savefig


def exp_recombination():
    K = 500
    wordlen = 8
    A = Alphabet('ACGT')

    # NOTE I can drive sensitivity to 0 and get decent results
    WB_kw = {'g_max': .2, 'sensitivity': .9, 'alphabet': A, 'wordlen': wordlen,
             'path': ':memory:', 'log_level': logging.INFO}

    homs = [rand_seq(A, i) for i in [i * K for i in range(1, 5)]]
    ps = [.01, .06, .12, .18]
    Ms = [MutationProcess(A, subst_probs=p, ge_prob=p, go_prob=p) for p in ps]

    def junk(): return rand_seq(A, np.random.randint(K / 2, K))

    S = junk() + homs[0] + junk() + homs[1] + junk() + homs[3] + \
        junk() + homs[2] + junk() + homs[0] + junk()
    homs = [M.mutate(hom)[0] for hom, M in zip(homs, Ms)]
    T = junk() + homs[3] + junk() + homs[2] + junk() + homs[0] + \
        junk() + homs[1] + junk() + homs[2] + junk()

    fig = plt.figure(figsize=(9, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
    ax_seeds = plt.subplot(gs[0])
    ax_mapping = plt.subplot(gs[1])

    ax_mapping.plot([1, 1], [0, len(S)], lw=2, c='k', alpha=.8)
    ax_mapping.plot([2, 2], [0, len(T)], lw=2, c='k', alpha=.8)

    WB = WordBlot(S, T, **WB_kw)

    p_min = (1 - max(ps)) ** 2
    scored_seeds = WB.score_seeds(K)
    # convert to ij coordinates and leave only the H1 score
    scored_seeds = [(WB.to_ij_coordinates(*rec['seed']), rec['p'])
                    for rec in scored_seeds]
    plot_scored_seeds(ax_seeds, scored_seeds, extent=[0, 1], threshold=p_min)

    cmap = plt.cm.get_cmap('jet')
    for rec in WB.similar_segments(K_min=K, p_min=p_min):
        seg = rec['segment']
        (i_start, i_end), (j_start, j_end) = WB.to_ij_coordinates_seg(seg)
        i_ctr, j_ctr = (i_start + i_end) / 2, (j_start + j_end) / 2
        color = cmap((rec['p'] - p_min) / (1 - p_min))[:3]
        plot_similar_segment(ax_seeds, seg, lw=5, alpha=.4, c=color)
        ax_mapping.plot([1, 1], [i_start, i_end], lw=10, c=color, alpha=.3)
        ax_mapping.plot([2, 2], [j_start, j_end], lw=10, c=color, alpha=.3)
        ax_mapping.plot([1, 2], [i_ctr, j_ctr], lw=1, c=color, alpha=.7)

    ax_mapping.set_xticks([1, 2])
    ax_mapping.set_xticklabels(['sequence 1', 'sequence 2'], fontsize=8)
    ax_mapping.set_xlim(0, 3)
    ax_c = make_axes_locatable(ax_mapping).append_axes('right', size='4%',
                                                       pad=0.05)
    norm = matplotlib.colors.Normalize(vmin=p_min, vmax=1)
    matplotlib.colorbar.ColorbarBase(ax_c, cmap=cmap, norm=norm,
                                     orientation='vertical')

    adjust_pw_plot(ax_seeds, len(S), len(T))

    fig.tight_layout()
    savefig(fig, 'rearrangement_duplication.png')


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


def exp_conserved_sequences():
    gap = .15
    subst = .15
    K = 500
    wordlen = 8
    A = Alphabet('ACGT')
    WB_kw = {'g_max': gap, 'sensitivity': .9, 'alphabet': A,
             'wordlen': wordlen, 'path': ':memory:', 'log_level': logging.INFO}

    M_std = MutationProcess(A, subst_probs=.2, ge_prob=.15, go_prob=.15)
    M_con = MutationProcess(A, subst_probs=.05, ge_prob=.05, go_prob=.05)

    pieces = [rand_seq(A, 10000) for _ in range(3)]
    conserved = rand_seq(A, 1000)

    S = pieces[0] + conserved + pieces[1] + pieces[2]
    pieces_hom = [M_std.mutate(p)[0] for p in pieces]
    conserved_hom = M_con.mutate(conserved)[0]
    T = pieces_hom[0] + pieces_hom[1] + conserved_hom + pieces_hom[2]

    fig = plt.figure()
    gs = gridspec.GridSpec(1, 2, width_ratios=[30, 1])
    ax = plt.subplot(gs[0])
    ax_colorbar = plt.subplot(gs[1])

    cmap = plt.cm.get_cmap('jet')

    WB = WordBlot(S, T, **WB_kw)
    subst = (1 - gap) * (1 - subst)
    for rec in WB.similar_segments(K_min=K, p_min=subst):
        segment, scores, p_hat = rec['segment'], rec['scores'], rec['p']
        (d_min, d_max), (a_min, a_max) = segment
        log('homologous segment: %s, scores = %.2f, %.2f, p = %.2f' %
            (str(segment), scores[0], scores[1], p_hat))
        color = cmap(p_hat)[:3]
        plot_similar_segment(ax, segment, color=color, alpha=.3, lw=10)

    norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
    plot_seeds(ax, WB.seeds())
    adjust_pw_plot(ax, len(S), len(T))

    matplotlib.colorbar.ColorbarBase(ax_colorbar, cmap=cmap, norm=norm,
                                     orientation='vertical')

    fig.suptitle('Conserved sequences with variable sequence similarity',
                 fontsize=10)
    fig.tight_layout()
    savefig(fig, 'conserved_sequences.png')


if __name__ == '__main__':
    exp_recombination()
    exp_repeat_regions()
    exp_conserved_sequences()
