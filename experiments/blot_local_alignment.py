import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib
import logging
from matplotlib import pyplot as plt
from biseqt.blot import HomologyFinder
from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess
from util import plot_seeds, plot_scored_seeds
from util import plot_similar_segment, adjust_pw_plot
from util import log, savefig


def exp_local_alignment():
    gap = .2
    subst = .1
    K = 500
    wordlen = 8
    A = Alphabet('ACGT')

    # NOTE I can drive sensitivity to 0 and get decent results
    HF_kw = {'g_max': .2, 'sensitivity': .9, 'alphabet': A, 'wordlen': wordlen,
             'path': ':memory:', 'log_level': logging.INFO}

    M = MutationProcess(A, subst_probs=subst, ge_prob=gap, go_prob=gap)

    homs = [rand_seq(A, i) for i in [K/2, K, 2 * K, 4 * K]]

    def junk(): return rand_seq(A, np.random.randint(2 * K, 4 * K))

    S = junk() + homs[0] + junk() + homs[1] + junk() + homs[3] + \
        junk() + homs[2] + junk() + homs[2] + junk() + homs[0] + junk()
    homs = [M.mutate(homs[i])[0] for i in range(len(homs))]
    T = junk() + homs[3] + junk() + homs[2] + junk() + homs[0] + \
        junk() + homs[3] + junk() + homs[1] + junk() + homs[2] + junk()

    fig = plt.figure()
    gs = gridspec.GridSpec(1, 2, width_ratios=[30, 1])
    ax = plt.subplot(gs[0])
    ax_colorbar = plt.subplot(gs[1])

    log('finding homologies')
    HF = HomologyFinder(S, T, **HF_kw)
    match = (1 - gap) * (1 - subst)

    scored_seeds = HF.score_seeds(K_min=K, p_min=match)
    # convert to ij coordinates and leave only the H0 score
    scored_seeds = [(HF.to_ij_coordinates(d, a), s0)
                    for (d, a), (s0, s1) in scored_seeds]
    plot_scored_seeds(ax, ax_colorbar, scored_seeds)

    for segment, score, match_p in HF.similar_segments(K_min=K, p_min=match):
        if segment is None:
            continue
        (d_min, d_max), (a_min, a_max) = segment
        log('homologous segment %s, score = %.2f' % (str(segment), score))
        plot_similar_segment(ax, segment, lw=5, alpha=.2)

    adjust_pw_plot(ax, len(S), len(T))

    fig.tight_layout()
    savefig(fig, 'local_alignment.png')


def exp_conserved_sequences():
    gap = .15
    subst = .15
    K = 500
    wordlen = 8
    A = Alphabet('ACGT')
    HF_kw = {'g_max': gap, 'sensitivity': .9, 'alphabet': A,
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

    log('finding homologies')
    HF = HomologyFinder(S, T, **HF_kw)
    scores = []
    subst = (1 - gap) * (1 - subst)
    for segment, score, match_p in HF.similar_segments(K_min=K, p_min=subst):
        (d_min, d_max), (a_min, a_max) = segment
        log('homologous segment: %s, score = %.2f, estim. match = %.2f' %
            (str(segment), score, match_p))
        color = cmap(match_p)[:3]
        plot_similar_segment(ax, segment, color=color, alpha=.3)
        scores.append(score)

    norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
    plot_seeds(ax, HF.seeds())
    adjust_pw_plot(ax, len(S), len(T))

    matplotlib.colorbar.ColorbarBase(ax_colorbar, cmap=cmap, norm=norm,
                                     orientation='vertical')

    fig.tight_layout()
    savefig(fig, 'conserved_sequences.png')


if __name__ == '__main__':
    exp_local_alignment()
    exp_conserved_sequences()
