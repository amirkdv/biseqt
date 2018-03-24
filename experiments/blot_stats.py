import numpy as np
from matplotlib import pyplot as plt
from itertools import product
from scipy.ndimage.filters import gaussian_filter1d

import logging

from biseqt.blot import HomologyFinder, find_peaks
from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess

from util import plot_with_sd, color_code
from util import with_dumpfile, log, savefig, load_fasta
from util import seq_pair, sample_bio_opseqs, sample_bio_seqs, apply_opseq

from util import plot_global_alignment
from util import adjust_pw_plot
from util import estimate_match_probs_in_opseq, fill_in_unknown
from util import plot_scored_seeds


@with_dumpfile
def sim_stats_fixed_K(K, ns, n_samples, **kw):
    def _zero():
        shape = (len(ns), n_samples)
        return {'pos': np.zeros(shape), 'neg': np.zeros(shape)}

    A = Alphabet('ACGT')
    wordlen, bio, p_match = kw['wordlen'], kw['bio'], kw['p_match']
    HF_kw = {
        'g_max': kw.get('g_max', .6),
        'sensitivity': kw.get('sensitivity', .99),
        'wordlen': wordlen,
        'alphabet': A,
        'path': kw.get('db_path', ':memory:'),
        'log_level': kw.get('log_level', logging.WARNING),
    }
    sim_data = {
        'scores': {
            'band': {'H0': _zero(), 'H1': _zero()},
            'segment': {'H0': _zero(), 'H1': _zero()},
        },
        'HF_kw': HF_kw,
        'K': K,
        'ns': ns,
        'bio': bio,
        'p_match': p_match,
    }
    if bio:
        bio_source_seq = kw['bio_source_seq']
        bio_source_opseq = kw['bio_source_opseq']
        opseqs = list(sample_bio_opseqs(bio_source_opseq, K, n_samples,
                                        gap=None, match=p_match))
    else:
        # NOTE distribute the burden of p_match evenly on gap and subst:
        subst = gap = 1 - np.sqrt(p_match)
        assert abs((1 - gap) * (1 - subst) - p_match) < 1e-3
        M = MutationProcess(A, subst_probs=subst, ge_prob=gap, go_prob=gap)

    for n_idx, n in enumerate(ns):
        log('evaluating scores for K = %d, n = %d (%d/%d)' %
            (K, n, n_idx + 1, len(ns)))
        for idx in range(n_samples):
            if bio:
                S_urel, T_urel = sample_bio_seqs([n, n], bio_source_seq)
                S_rel = sample_bio_seqs([K], bio_source_seq)[0]
                T_rel = apply_opseq(S_rel, opseqs[idx])
                junk = sample_bio_seqs([n - K, n - K], bio_source_seq)
                S_rel += junk[0]
                T_rel += junk[1]
            else:
                S_rel, T_rel = seq_pair(K, A, mutation_process=M)
                S_rel += rand_seq(A, n - K)
                T_rel += rand_seq(A, n - K)
                S_urel, T_urel = rand_seq(A, n), rand_seq(A, n)

            for key, (S, T) in zip(['pos', 'neg'],
                                   [(S_rel, T_rel), (S_urel, T_urel)]):
                HF = HomologyFinder(S, T, **HF_kw)

                # band score
                # NOTE assume exact knowledge of K and p_match
                scores_by_d_ = HF.score_diagonal_bands(K, p_match)
                s0, s1 = scores_by_d_[HF.d0]
                sim_data['scores']['band']['H0'][key][n_idx][idx] = s0
                sim_data['scores']['band']['H1'][key][n_idx][idx] = s1

                # segment score
                d_radius = int(np.ceil(HF.band_radius(K)))
                d_band = (-d_radius, d_radius)
                _, scores_by_a = HF.score_segments(K, p_match, d_band=d_band)
                s0, s1 = scores_by_a[int(K/2)]
                sim_data['scores']['segment']['H0'][key][n_idx][idx] = s0
                sim_data['scores']['segment']['H1'][key][n_idx][idx] = s1
    return sim_data


@with_dumpfile
def sim_stats_varying_K_p(Ks, ps, n_samples, **kw):
    def _zero():
        shape = (len(Ks), len(ps), n_samples)
        return {'pos': np.zeros(shape), 'neg': np.zeros(shape)}

    A = Alphabet('ACGT')
    wordlen, bio = kw['wordlen'], kw['bio']
    HF_kw = {
        'g_max': kw.get('g_max', .6),
        'sensitivity': kw.get('sensitivity', .99),
        'wordlen': wordlen,
        'alphabet': A,
        'path': kw.get('db_path', ':memory:'),
        'log_level': kw.get('log_level', logging.WARNING),
    }
    sim_data = {
        'scores': {'H0': _zero(), 'H1': _zero()},
        'K_hat': {'H0': _zero(), 'H1': _zero()},
        'p_hat': {'H0': _zero(), 'H1': _zero()},
        'HF_kw': HF_kw,
        'Ks': Ks,
        'ps': ps,
        'bio': bio,
    }
    if bio:
        bio_source_seq = kw['bio_source_seq']
        bio_source_opseq = kw['bio_source_opseq']

    for (K_idx, K), (p_idx, p_match) in product(enumerate(Ks), enumerate(ps)):
        if bio:
            opseqs = list(sample_bio_opseqs(bio_source_opseq, K, n_samples,
                                            gap=None, match=p_match))
        log('evaluating scores for K = %d, p = %.2f' % (K, p_match))
        for idx in range(n_samples):
            if bio:
                S_urel, T_urel = sample_bio_seqs([K, K], bio_source_seq)
                S_rel = sample_bio_seqs([K], bio_source_seq)[0]
                T_rel = apply_opseq(S_rel, opseqs[idx])
                junk = sample_bio_seqs([K / 2] * 4, bio_source_seq)
                S_rel = junk[0] + S_rel + junk[1]
                T_rel = junk[2] + T_rel + junk[3]
            else:
                # NOTE distribute the burden of p_match evenly on gap and subst
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
                HF = HomologyFinder(S, T, **HF_kw)
                scored_seeds = HF.score_seeds(K, p_match)

                if scored_seeds:
                    s0 = max(s0 for _, _, (s0, s1) in scored_seeds)
                    s1 = max(s1 for _, _, (s0, s1) in scored_seeds)
                else:
                    s0 = s1 = float('-inf')
                sim_data['scores']['H0'][key][K_idx][p_idx][idx] = s0
                sim_data['scores']['H1'][key][K_idx][p_idx][idx] = s1

                if key == 'neg':
                    continue
                # estimate segment probabilities with unknown K, p
                p_min = .4
                K_min = 100
                for mode in ['H0', 'H1']:
                    results = list(HF.similar_segments(K_min, p_min,
                                                       mode=mode))
                    scores = [score for _, score, _ in results]
                    if not scores:
                        continue
                    segment, _, p_hat = results[np.argmax(scores)]
                    _, (a_min, a_max) = segment
                    K_hat = a_max - a_min
                    sim_data['K_hat'][mode][key][K_idx][p_idx][idx] = K_hat
                    sim_data['p_hat'][mode][key][K_idx][p_idx][idx] = p_hat
    return sim_data


@with_dumpfile
def sim_stats_real_homologies(seqs, pws, **kw):
    A = Alphabet('ACGT')
    wordlen, K_min, p_min = kw['wordlen'], kw['K_min'], kw['p_min']
    HF_kw = {
        'g_max': kw.get('g_max', .6),
        'sensitivity': kw.get('sensitivity', .99),
        'wordlen': wordlen,
        'alphabet': A,
        'path': kw.get('db_path', ':memory:'),
        'log_level': kw.get('log_level', logging.WARNING),
        # HACK mask CG-rich and homopolymeric regions
        'mask': [set(x) for x in [[0], [1], [2], [3], [1, 2], [0, 3]]],
    }
    similar_segments_kw = {'K_min': K_min, 'p_min': p_min}
    sim_data = {
        'pws': pws,
        'seqlens': {name: len(seqs[name]) for name in seqs},
        'similarity': {key: [np.zeros(len(seqs[key[0]])),
                             np.zeros(len(seqs[key[1]]))]
                       for key in pws},
        'seeds': {key: [] for key in pws},
        'HF_kw': HF_kw,
        'similar_segments_kw': similar_segments_kw
    }
    for idx, (id1, id2) in enumerate(pws):
        log('finding local homologies between %s (%d) and %s (%d)' %
            (id1, sim_data['seqlens'][id1], id2, sim_data['seqlens'][id2]))
        S = seqs[id1]
        T = seqs[id2]

        HF = HomologyFinder(S, T, **HF_kw)
        sim_data['seeds'][(id1, id2)] = {
            HF.to_ij_coordinates(*rec['seed']): rec['p']
            for rec in HF.score_seeds(K_min)
        }
        log('-> found %d exactly matching %d-mers' %
            (len(sim_data['seeds'][(id1, id2)]), wordlen))
        for (pos1, pos2), p_hat in sim_data['seeds'][(id1, id2)].items():
            res = sim_data['similarity'][(id1, id2)][0]
            for i in range(pos1, pos1 + wordlen):
                res[i] = max(res[i], p_hat)
            res = sim_data['similarity'][(id1, id2)][1]
            for i in range(pos2, pos2 + wordlen):
                res[i] = max(res[i], p_hat)
    return sim_data


def plot_stats_fixed_K(sim_data, suffix=''):
    K, ns, p_match = sim_data['K'], sim_data['ns'], sim_data['p_match']
    wordlen = sim_data['HF_kw']['wordlen']
    n_samples = sim_data['scores']['band']['H0']['pos'].shape[1]
    scores = sim_data['scores']

    for key in ['band', 'segment']:
        fig = plt.figure(figsize=(11, 5))
        ax_H0 = fig.add_subplot(1, 2, 1)
        ax_H1 = fig.add_subplot(1, 2, 2)
        scores = sim_data['scores'][key]

        kw = {'marker': 'o', 'markersize': 3, 'alpha': .8}
        for mode, ax in zip(['H0', 'H1'], [ax_H0, ax_H1]):
            for case, color in zip(['pos', 'neg'], ['g', 'r']):
                label = ('' if case == 'pos' else 'non-') + 'homologous'
                plot_with_sd(ax, ns, scores[mode][case], axis=1, color=color,
                             label=label, **kw)
            ax.set_ylabel('%s %s score' % (mode, key), fontsize=10)
            ax.set_xlabel('sequence length', fontsize=10)
        ax_H0.legend(loc='upper right', fontsize=10)

        fig.suptitle('K = %d, w = %d, no. samples = %d, match prob= %.2f' %
                     (K, wordlen, n_samples, p_match), fontsize=10)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        savefig(fig, 'stats[%s]%s.png' % (key, suffix))


def plot_stats_varying_K_p(sim_data, select_Ks, select_ps, suffix=''):
    Ks, ps = sim_data['Ks'], sim_data['ps']
    wordlen = sim_data['HF_kw']['wordlen']
    n_samples = sim_data['scores']['H0']['pos'].shape[2]
    scores = sim_data['scores']
    K_hats = sim_data['K_hat']
    p_hats = sim_data['p_hat']

    assert all(K in Ks for K in select_Ks)
    assert all(p in ps for p in select_ps)

    kw = {'marker': 'o', 'markersize': 3, 'alpha': .8}

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
        ax.set_xlabel('homology length', fontsize=10)
    ax_H0.legend(loc='upper left', fontsize=10)
    ax_H1.legend(loc='lower left', fontsize=10)
    fig_by_K.suptitle('w = %d, no. samples = %d' % (wordlen, n_samples),
                      fontsize=10)
    fig_by_K.tight_layout(rect=[0, 0.03, 1, 0.95])
    savefig(fig_by_K, 'stats[score-by-K]%s.png' % suffix)

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
        ax.set_xlabel('homology match probability', fontsize=10)
    ax_H0.legend(loc='upper left', fontsize=10)
    ax_H1.legend(loc='lower left', fontsize=10)
    fig_by_p.suptitle('w = %d, no. samples = %d' % (wordlen, n_samples),
                      fontsize=10)
    fig_by_p.tight_layout(rect=[0, 0.03, 1, 0.95])
    savefig(fig_by_p, 'stats[score-by-p]%s.png' % suffix)

    # estimated Ks for select ps
    fig_K_hat = plt.figure(figsize=(11, 5))
    ax_H0 = fig_K_hat.add_subplot(1, 2, 1)
    ax_H1 = fig_K_hat.add_subplot(1, 2, 2)
    for mode, ax in zip(['H0', 'H1'], [ax_H0, ax_H1]):
        colors = color_code(select_ps)
        for p, color in zip(select_ps, colors):
            p_idx = ps.index(p)
            plot_with_sd(ax, Ks, K_hats[mode]['pos'][:, p_idx, :], axis=1,
                         color=color, label='p = %.2f' % p, **kw)
        ax.set_ylabel('%s estimated homology length' % mode, fontsize=10)
        ax.set_xlabel('true homology length', fontsize=10)
        ax.plot(Ks, Ks, ls='--', c='k', lw=5, alpha=.2)
        ax.legend(loc='upper left', fontsize=10)
    fig_K_hat.suptitle('w = %d, no. samples = %d' % (wordlen, n_samples),
                       fontsize=10)
    fig_K_hat.tight_layout(rect=[0, 0.03, 1, 0.95])
    savefig(fig_K_hat, 'stats[K-hat]%s.png' % suffix)

    # estimated ps for select Ks
    fig_p_hat = plt.figure(figsize=(11, 5))
    ax_H0 = fig_p_hat.add_subplot(1, 2, 1)
    ax_H1 = fig_p_hat.add_subplot(1, 2, 2)
    for mode, ax in zip(['H0', 'H1'], [ax_H0, ax_H1]):
        colors = color_code(select_Ks)
        for K, color in zip(select_Ks, colors):
            K_idx = Ks.index(K)
            plot_with_sd(ax, ps, p_hats[mode]['pos'][K_idx, :, :], axis=1,
                         color=color, label='K = %d' % K, **kw)
        ax.set_ylabel('%s estimated match probability' % mode,
                      fontsize=10)
        ax.set_xlabel('true match probability', fontsize=10)
        ax.plot(ps, ps, ls='--', c='k', lw=5, alpha=.2)
        ax.set_ylim(0, 1)
        ax.legend(loc='upper left', fontsize=10)
    fig_p_hat.suptitle('w = %d, no. samples = %d' % (wordlen, n_samples),
                       fontsize=10)
    fig_p_hat.tight_layout(rect=[0, 0.03, 1, 0.95])
    savefig(fig_p_hat, 'stats[p-hat]%s.png' % suffix)


# NOTE ROC for actb is, but for acta2 it's TERRIBLE because there
# are a lot more similarities than there are in the alignment. We can slightly
# offset this by increasing K_min dramatically (e.g. upto 30,000) but instead
# we just show the distribution of scores of "real seeds" across different
# pairs of sequences to make a case for the score being "stable" (i.e. one
# cutoff works for all situations).
def plot_stats_real_homologies(sim_data, suffix=''):
    seqlens = sim_data['seqlens']
    pws = sim_data['pws']
    K_min = sim_data['similar_segments_kw']['K_min']
    fig_num = int(np.ceil(np.sqrt(len(pws))))

    fig_seeds = plt.figure(figsize=(6 * fig_num, 5 * fig_num))
    fig_sim = plt.figure(figsize=(8 * fig_num, 4 * fig_num))
    fig_cons = plt.figure(figsize=(6, 4))
    ax_cons = fig_cons.add_subplot(1, 1, 1)

    colors_cons = color_code(range(len(pws)))
    p_cons = .7

    def _conserved(ps):
        return sum([range(*r) for r in find_peaks(ps, 10, p_cons)], [])

    for idx, ((id1, id2), opseq) in enumerate(pws.items()):
        ax_seeds = fig_seeds.add_subplot(fig_num, fig_num, idx + 1)

        ax_sim = fig_sim.add_subplot(fig_num, fig_num, idx + 1)
        radius = K_min / 2

        ps_hat = sim_data['similarity'][(id1, id2)][0]
        ps_true = estimate_match_probs_in_opseq(opseq, radius, projection=1)

        log(id1 + '/' + id2)
        cons_hat = _conserved(ps_hat)
        cons_true = _conserved(ps_true)
        tp = len(set(cons_hat).intersection(set(cons_true)))
        tpr = 1. * tp / len(cons_true)
        ppv = 1. * tp / len(cons_hat)
        print tpr, ppv
        color = colors_cons[idx]
        ax_cons.scatter([tpr], [ppv], s=20, c=color, lw=0, alpha=.4,
                        label='%s/%s' % (id1.split('.')[0], id2.split('.')[0]))

        # smooth for ease of visual inspection
        ps_hat = gaussian_filter1d(ps_hat, radius)
        ps_true = gaussian_filter1d(ps_true, radius)
        ax_sim.plot(range(len(ps_hat)), ps_hat, c='k', alpha=.8, lw=1)
        ax_sim.plot(range(len(ps_true)), ps_true, c='g', alpha=.3, lw=2)
        ax_sim.set_title('%s/%s' % (id1, id2))
        # ax_sim.fill_between(range(len(ps_true)), ps_true,
        #                     where=ps_true >= p_cons, color='g', alpha=.1)
        ax_sim.set_xlabel('position in %s' % id1)
        ax_sim.set_ylabel('match probability (window: %d)' % K_min)

        plot_scored_seeds(ax_seeds,
                          sim_data['seeds'][(id1, id2)].items(),
                          threshold=.7)
        plot_global_alignment(ax_seeds, opseq, c='k', lw=5, alpha=.2)
        adjust_pw_plot(ax_seeds, seqlens[id1], seqlens[id2])
        ax_seeds.set_ylabel(id1, fontsize=10)
        ax_seeds.set_xlabel(id2, fontsize=10)

    ax_cons.set_xlim(-.1, 1.1)
    ax_cons.set_ylim(-.1, 1.1)
    ax_cons.legend(loc='best', fontsize=4)
    ax_cons.set_xlabel('True positive rate')
    ax_cons.set_ylabel('Positive predictive value')
    ax_cons.set_title('Highly conserved ($\hat{p} > %.2f$) regions' % p_cons)
    savefig(fig_seeds, 'real_homologies[seeds]%s.png' % suffix)
    savefig(fig_sim, 'real_homologies[similarity]%s.png' % suffix)
    savefig(fig_cons, 'real_homologies[conservation]%s.png' % suffix)


# the point of this experiment is: band score is unreliable and segment score
# is reliable. After this experiment we exclusively look at segment score.
def exp_stats_performance_fixed_K():
    K = 200  # similar segment length
    ns = [K * 2 ** i for i in range(8)]  # sequence lengths
    n_samples = 50  # number samples for each n

    wordlen = 8
    p_match = .8
    suffix = '[K=%d]' % K
    dumpfile = 'stats%s.txt' % suffix
    sim_data = sim_stats_fixed_K(
        K, ns, n_samples, p_match=p_match, wordlen=wordlen, bio=False,
        dumpfile=dumpfile, ignore_existing=False)
    plot_stats_fixed_K(sim_data, suffix=suffix)

    # ===============
    # Biological Data
    # ===============
    suffix = '[K=%d][bio]' % K
    dumpfile = 'stats%s.txt' % suffix
    A = Alphabet('ACGT')
    with open('data/acta2/acta2-7vet.fa') as f:
        bio_source_seq = ''.join(s for s, _, _ in load_fasta(f))
        bio_source_seq = bio_source_seq.upper()
        bio_source_seq = ''.join(x for x in bio_source_seq if x in 'ACGT')
        bio_source_seq = A.parse(bio_source_seq)
    with open('data/acta2/acta2_opseqs.fa') as f:
        bio_source_opseq = ''.join(s for s, _, _ in load_fasta(f))
    sim_data = sim_stats_fixed_K(
        K, ns, n_samples, p_match=p_match, wordlen=wordlen, bio=True,
        bio_source_opseq=bio_source_opseq, bio_source_seq=bio_source_seq,
        dumpfile=dumpfile, ignore_existing=False)
    plot_stats_fixed_K(sim_data, suffix=suffix)


def exp_stats_performance_varying_K_p():
    Ks = [200 * i for i in range(1, 9)]
    select_Ks = Ks[0], Ks[2], Ks[4]

    ps = [1 - .06 * i for i in range(1, 9)]
    select_ps = ps[0], ps[2], ps[4]

    n_samples = 50  # number samples for each n

    wordlen = 8
    suffix = '[varying-K-p]'
    dumpfile = 'stats%s.txt' % suffix
    sim_data = sim_stats_varying_K_p(
        Ks, ps, n_samples, wordlen=wordlen, bio=False,
        dumpfile=dumpfile, ignore_existing=False)
    plot_stats_varying_K_p(sim_data, select_Ks, select_ps, suffix=suffix)

    # ===============
    # Biological Data
    # ===============
    suffix = '[varying-K-p][bio]'
    dumpfile = 'stats%s.txt' % suffix
    A = Alphabet('ACGT')
    with open('data/acta2/acta2-7vet.fa') as f:
        bio_source_seq = ''.join(s for s, _, _ in load_fasta(f))
        bio_source_seq = bio_source_seq.upper()
        bio_source_seq = ''.join(x for x in bio_source_seq if x in 'ACGT')
        bio_source_seq = A.parse(bio_source_seq)
    with open('data/acta2/acta2_opseqs.fa') as f:
        bio_source_opseq = ''.join(s for s, _, _ in load_fasta(f))
    sim_data = sim_stats_varying_K_p(
        Ks, ps, n_samples, wordlen=wordlen, bio=True,
        bio_source_opseq=bio_source_opseq, bio_source_seq=bio_source_seq,
        dumpfile=dumpfile, ignore_existing=False)
    plot_stats_varying_K_p(sim_data, select_Ks, select_ps, suffix=suffix)


def exp_stats_real_homologies():
    p_min = .5
    A = Alphabet('ACGT')

    wordlen = 5
    K_min = 100
    suffix = '[actb]'
    dumpfile = 'real_homologies%s.txt' % suffix
    seqs_path = 'data/actb/actb-7vet.fa'
    pws_path = 'data/actb/actb-7vet-pws.fa'

    # wordlen = 12
    # K_min = 200
    # suffix = '[acta2]'
    # dumpfile = 'real_homologies%s.txt' % suffix
    # seqs_path = 'data/acta2/acta2-7vet.fa'
    # pws_path = 'data/acta2/acta2-7vet-pws.fa'

    # wordlen = 8
    # K_min = 50
    # suffix = '[acta1]'
    # dumpfile = 'real_homologies%s.txt' % suffix
    # seqs_path = 'data/acta1/acta1-7vet.fa'
    # pws_path = 'data/acta1/acta1-7vet-pws.fa'

    # wordlen = 8
    # K_min = 50
    # suffix = '[ngf]'
    # dumpfile = 'real_homologies%s.txt' % suffix
    # seqs_path = 'data/ngf/ngf-7vet.fa'
    # pws_path = 'data/ngf/ngf-7vet-pws.fa'

    # NOTE anything "acta2" is actually "actn2"
    # FIXME the acta2 thing we're using is not actually a gene (what is it?)
    # also, all our examples (aside from acta2) are short. I'm not sure
    # if we're looking at full genes or just exons because the
    # extents on genome browser suggest longer sequences than we actually get
    # from maf (assuming our conversion code is correct)
    # names = ['hg38.chr1', 'rn5.chr17', 'canFam3.chr4', 'mm10.chr13']
    # names = ['hg38.chr1', 'panTro4.chr1', 'canFam3.chr4']

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
    plot_stats_real_homologies(sim_data, suffix=suffix)


if __name__ == '__main__':
    exp_stats_performance_fixed_K()
    exp_stats_performance_varying_K_p()
    exp_stats_real_homologies()
