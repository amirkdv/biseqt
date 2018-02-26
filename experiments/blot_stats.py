import numpy as np
from matplotlib import pyplot as plt
from itertools import product

import logging

from biseqt.blot import HomologyFinder
from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess

from util import plot_with_sd, color_code
from util import with_dumpfile, log, savefig, load_fasta
from util import seq_pair, sample_bio_opseqs, sample_bio_seqs, apply_opseq

from util import plot_global_alignment, plot_similar_segment
from util import adjust_pw_plot, plot_seeds, opseq_path


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
                    s0 = max(s0 for _, (s0, s1) in scored_seeds)
                    s1 = max(s1 for _, (s0, s1) in scored_seeds)
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
    mode, threshold = kw['mode'], kw['threshold']
    wordlen, K_min, p_min = kw['wordlen'], kw['K_min'], kw['p_min']
    HF_kw = {
        'g_max': kw.get('g_max', .3),
        'sensitivity': kw.get('sensitivity', .99),
        'wordlen': wordlen,
        'alphabet': A,
        'path': kw.get('db_path', ':memory:'),
        'log_level': kw.get('log_level', logging.WARNING),
    }
    similar_segments_kw = {
        'K_min': K_min, 'p_min': p_min, 'mode': mode, 'threshold': threshold,
    }
    sim_data = {
        'pws': pws,
        'seqlens': {name: len(seqs[name]) for name in seqs},
        'segments': {key: [] for key in pws},
        'seeds': {key: [] for key in pws},
        'HF_kw': HF_kw,
        'similar_segments_kw': similar_segments_kw
    }
    for idx, (id1, id2) in enumerate(pws):
        log('finding all local homologies between %s and %s' %
            (id1, id2))
        S = seqs[id1]
        T = seqs[id2]
        HF = HomologyFinder(S, T, **HF_kw)
        sim_data['seeds'][(id1, id2)] = list(HF.seeds())
        for segment, _, _ in HF.similar_segments(**similar_segments_kw):
            sim_data['segments'][(id1, id2)].append(segment)
            # print segment
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


def plot_stats_real_homologies(sim_data, suffix=''):
    seqlens = sim_data['seqlens']
    pws = sim_data['pws']
    fig_num = int(np.ceil(np.sqrt(len(pws))))
    fig = plt.figure(figsize=(6*fig_num, 5*fig_num))
    tpr, ppv = {}, {}
    for idx, (id1, id2) in enumerate(pws):
        ax = fig.add_subplot(fig_num, fig_num, idx+1)
        opseq = pws[(id1, id2)]
        segments = sim_data['segments'][(id1, id2)]
        ms_len = opseq.count('M') + opseq.count('S')
        hom_len = sum(a_max - a_min for _, (a_min, a_max) in segments)
        path = zip(*opseq_path(opseq))
        ms_in_seg = 0
        for op, (i, j) in zip(opseq, path):
            if op not in 'MS':
                continue
            d, a = i - j, min(i, j)
            for (d_min, d_max), (a_min, a_max) in segments:
                if d_min <= d <= d_max and a_min <= a <= a_max:
                    ms_in_seg += 1
        tpr[(id1, id2)] = ms_in_seg * 1. / ms_len
        ppv[(id1, id2)] = ms_in_seg * 1. / hom_len
        log('comparison between %s and %s, tpr = %.2f, ppv = %.2f' %
            (id1, id2, tpr[(id1, id2)], ppv[(id1, id2)]))

        for segment in segments:
            plot_similar_segment(ax, segment, lw=3, alpha=.2)
        plot_seeds(ax, sim_data['seeds'][(id1, id2)])
        plot_global_alignment(ax, opseq, c='g', lw=5, alpha=.3)
        adjust_pw_plot(ax, seqlens[id1], seqlens[id2])
        ax.set_xlabel(id1, fontsize=10)
        ax.set_ylabel(id2, fontsize=10)
    savefig(fig, 'real_homologies[seeds]%s.png' % suffix)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    pairs = pws.keys()
    colors = color_code(range(len(pairs)), cmap='brg')
    for color, pair in zip(colors, pairs):
        ax.scatter([tpr[pair]], [ppv[pair]], s=30, alpha=.8, c=color, lw=0)
    for color, pair in zip(colors, pairs):
        wobble = np.random.randint(-10, 10) / 500.
        rotation = [0, 90][np.random.randint(0, 2)]
        label = '/'.join(p.split('.')[0] for p in pair)
        # dx = wobble if rotation == 0 else -.04
        # dy = wobble if rotation == 90 else .01
        rotation = 0
        dx = .01
        dy = wobble
        ax.text(tpr[pair] + dx, ppv[pair] + dy, label,
                color=color, fontsize=6, alpha=.5, rotation=rotation)
    ax.set_xlabel('true positive rate')
    ax.set_ylabel('positive predictive value')
    ax.set_xlim(-.1, 1.2)
    ax.set_ylim(-.1, 1.2)
    savefig(fig, 'real_homologies[performance]%s.png' % suffix, dpi=500)


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
    wordlen = 8
    p_min = .6
    K_min = 50
    threshold = 1.5
    mode = 'H1'
    A = Alphabet('ACGT')

    suffix = '[actb]'
    dumpfile = 'real_homologies%s.txt' % suffix
    seqs = 'data/actb/actb-7vet.fa'
    pws = 'data/actb/actb-7vet-pws.fa'

    # suffix = '[acta2]'
    # dumpfile = 'real_homologies%s.txt' % suffix
    # seqs = 'data/acta2/acta2-7vet.fa'
    # pws = 'data/acta2/acta2-7vet-pws.fa'

    with open(seqs) as f:
        seqs = {name: A.parse(seq.upper())
                for seq, name, _ in load_fasta(f)}
    with open(pws) as f:
        pws = {tuple(name.split(':')): seq for seq, name, _ in load_fasta(f)}
    sim_data = sim_stats_real_homologies(
        seqs, pws, wordlen=wordlen, p_min=p_min, K_min=K_min, mode=mode,
        threshold=threshold, dumpfile=dumpfile, ignore_existing=False)
    plot_stats_real_homologies(sim_data, suffix=suffix)


if __name__ == '__main__':
    exp_stats_performance_fixed_K()
    exp_stats_performance_varying_K_p()
    exp_stats_real_homologies()
