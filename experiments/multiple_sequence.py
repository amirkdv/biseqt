#!/usr/bin/env python
import numpy as np
import sys
import logging
from itertools import product
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D
from Bio import AlignIO
from time import time

from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess
from biseqt.blot import WordBlotMultiple
from biseqt.blot import WordBlotMultipleFast
from util import log, get_seqs_from_mse, with_dumpfile, fill_in_unknown
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
def plot_similar_segment_3d(ax, segment, **kw):
    assert len(segment) == 3
    assert all(len(range) == 2 for range in segment)
    ax.plot(segment[0], segment[1], segment[2], **kw)


# ========================================================
# Biological data: compare aligned genomes
# ========================================================
def seeds_from_maf(mse_path, wordlen, fmt='maf', ids=None):
    alignments = list(AlignIO.parse(mse_path, fmt))
    # NOTE trust the first alignment to have all the ids
    all_ids = set()
    for alignment in alignments:
        all_ids = all_ids.union(set(rec.id for rec in alignment))
    if ids is None:
        ids = all_ids
    else:
        assert set(ids).issubset(all_ids)
    assert len(ids)
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
    """Multiple sequence similarities with Word-Blot on *biological data*.

    **Supported Claimes**

    * Word-Blot estimated match probabilities are good classifiers for
      homologous and non-homologous seeds in multiple sequences.

    .. figure::
        https://www.dropbox.com/s/qhgl9vtomhfzomx/
        multiple-sequence%5Bactb%5D.png?raw=1
       :target:
        https://www.dropbox.com/s/qhgl9vtomhfzomx/
        multiple-sequence%5Bactb%5D.png?raw=1
       :alt: lightbox

       Multiple sequence comparison on *Beta Actin (ACTB)* gene on the 7
       vertebrate USCS dataset (excluding *monodelphis domestica* in which the
       gene is split between chromosomes 1 and 6) with word length 6, minimum
       similarity length 100nt, and minimum match probability 0.75. For an
       arbitrary group of 3 sequences a dot plot similar to those of pairwise
       comparisons is shown (*left*). Seeds are color coded by intensity to
       represent their estimated match probability, blue ligns show local
       similarities identified by Word-Blot. The performance of estimated match
       probabilities in the full data set to discriminate between homologous
       and non-homologous seeds is shown by an ROC curve (*right*).

    .. figure::
        https://www.dropbox.com/s/96xf8p42gbazjyw/
        multiple-sequence%5Birx1%5D.png?raw=1
       :target:
        https://www.dropbox.com/s/96xf8p42gbazjyw/
        multiple-sequence%5Birx1%5D.png?raw=1
       :alt: lightbox

       Same as figure above but with the *IRX1* data set from Ensembl
       (excluding *rattus norvegicus* whose gene copy is much shorter than the
       rest) with word length 10.
    """
    # mse_path = 'data/actb/actb-7vet.maf'
    # gene = 'ACTB'
    # mse_fmt = 'maf'
    # mse_names = 'ucsc'
    # suffix = '[%s]' % gene.lower()
    # wordlen = 6
    # ids_to_exclude = ['monDom5.chr1', 'monDom5.chr6']  # split gene
    # scatter_ids = ['hg38.chr7', 'panTro4.chr7', 'canFam3.chr6']

    mse_path = 'data/irx1/irx1-vert-amniota.aln'
    gene = 'IRX1'
    mse_fmt = 'clustal'
    mse_names = 'ensembl'
    wordlen = 10
    ids_to_exclude = ['rattus_norvegicus/1-10720']  # much shorter than rest
    scatter_ids = ['homo_sapiens/1-10720', 'canis_familiaris/1-10720',
                   'bos_taurus/1-10720']

    p_min = .75
    K_min = 100

    A = Alphabet('ACGT')
    WB_kw = {'g_max': .4, 'sensitivity': .9, 'alphabet': A,
             'wordlen': wordlen, 'path': ':memory:'}

    # 3 sequences for scatter plot
    assert len(scatter_ids) == 3
    seqs = {id_: A.parse(fill_in_unknown(seq.upper(), A))
            for id_, seq in get_seqs_from_mse(mse_path, fmt=mse_fmt)
            if id_ not in ids_to_exclude}
    ids = seqs.keys()
    WB = WordBlotMultipleFast(*[seq for id_, seq in seqs.items()
                                if id_ in scatter_ids], **WB_kw)
    scored_seeds = [(WB.to_ij_coordinates(*rec['seed']), rec['p'])
                    for rec in WB.score_seeds(K_min)]
    log('found %d seeds for 3 sequences' % len(scored_seeds))

    fig = plt.figure(figsize=(8, 5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1.3])
    ax_scatter = plt.subplot(gs[0], projection=Axes3D.name)
    plot_scored_seeds_3d(fig, ax_scatter, scored_seeds, threshold=p_min)
    for axis in 'xyz':
        ax_scatter.tick_params(axis=axis, labelsize=4)

    if mse_names == 'ucsc':
        ax_scatter.set_xlabel(scatter_ids[0].split('.')[0])
        ax_scatter.set_ylabel(scatter_ids[1].split('.')[0])
        ax_scatter.set_zlabel(scatter_ids[2].split('.')[0])
    if mse_names == 'ensembl':
        ax_scatter.set_xlabel(scatter_ids[0].split('/')[0].replace('_', ' '))
        ax_scatter.set_ylabel(scatter_ids[1].split('/')[0].replace('_', ' '))
        ax_scatter.set_zlabel(scatter_ids[2].split('/')[0].replace('_', ' '))
    ax_scatter.set_title('Estimated similarity at %d ' % len(scored_seeds) +
                         'seeds in 3 copiees of %s' % gene, fontsize=8)

    for res in WB.similar_segments(K_min=K_min, p_min=p_min, score=False):
        std_ranges = WB.to_ij_coordinates_seg(res['segment'])
        plot_similar_segment_3d(ax_scatter, std_ranges, c='b', lw=1, zorder=10)

    # ============================
    # hom/non-hom seed classifier
    # ============================
    real_seeds = list(seeds_from_maf(mse_path, wordlen, fmt=mse_fmt, ids=ids))
    real_seeds = list(tuple(int(x) for x in seed) for seed in real_seeds)
    log('found %d homologous seeds for %d sequences' %
        (len(real_seeds), len(ids)))

    WB = WordBlotMultipleFast(*[seq for id_, seq in seqs.items()
                                if id_ in ids], **WB_kw)
    scored_seeds = [(WB.to_ij_coordinates(*rec['seed']), rec['p'])
                    for rec in WB.score_seeds(K_min)]
    log('found %d seeds for %d sequences' % (len(scored_seeds), len(ids)))
    pos, neg = [], []
    for coords, p_hat in scored_seeds:
        if coords in real_seeds:
            pos.append(p_hat)
        else:
            neg.append(p_hat)

    ax_roc = plt.subplot(gs[1])
    plot_roc(ax_roc, pos, neg, color='k', lw=1)
    avg_len = sum(len(seq) for seq in seqs.values()) / len(seqs)
    title = 'Classifing %d(+) and %d(-) seeds\n' % (len(pos), len(neg))
    title += 'in %d copies of %s (avg. len. %dnt)' % (len(ids), gene, avg_len)
    ax_roc.set_title(title, fontsize=8)

    fig.tight_layout()
    savefig(fig, 'multiple-sequence[%s].png' % gene.lower())


# ========================================================
# Simulations (K, p, r, ds, a)
# ========================================================
@with_dumpfile
def sim_simulated_K_p(Ks, ps, **kw):
    # ps are the probability that ALL seqs match
    n_samples, n_seqs = kw['n_samples'], kw['n_seqs']
    shape = (len(Ks), len(ps), n_samples)

    A = Alphabet('ACGT')
    WB_kw = {
        'g_max': kw.get('g_max', .6),
        'sensitivity': kw.get('sensitivity', .99),
        'alphabet': A,
        'path': kw.get('db_path', ':memory:'),
        'log_level': kw.get('log_level', logging.WARN),
    }
    sim_data = {
        'n_seqs': n_seqs,
        'scores': {'H0': {'pos': np.zeros(shape), 'neg': np.zeros(shape)},
                   'H1': {'pos': np.zeros(shape), 'neg': np.zeros(shape)}},
        'K_hat': np.zeros(shape),
        'p_hat': np.zeros(shape),
        'a_hat': np.zeros(shape),
        'd_min': np.zeros(shape),
        'd_max': np.zeros(shape),
        'time': {'pos': np.zeros(shape), 'neg': np.zeros(shape)},
        'WB_kw': WB_kw,
        'Ks': Ks,
        'ps': ps,
    }
    # NOTE in multiple sequence case segments gets broken (along a) more easily
    # (this is not fragmentation which is breaking along ds). We use K / 2
    # instead.
    # K_min = 100
    # assert K_min <= min(Ks)
    K_min = 100

    for (K_idx, K), (p_idx, p_match) in product(enumerate(Ks), enumerate(ps)):
        p_min = p_match
        # adjusted wordlens
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

                # calculate H1/H1 scores with perfect information:
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

                # we only need scores from negative cases
                if key == 'neg':
                    continue

                def _len(seg): return (seg[1][1] - seg[1][0]) / n_seqs

                # for multiple sequences it's too much to ask for
                # at_least_one=True
                results = list(WB.similar_segments(K_min, p_min))
                if results:
                    # pick the longest detected homology
                    hom = max(results, key=lambda rec: _len(rec['segment']))
                    sim_data['p_hat'][K_idx, p_idx, idx] = hom['p']
                else:
                    # NOTE a fake segment with all zero properties
                    sim_data['p_hat'][K_idx, p_idx, idx] = 0
                    results = [{'segment': ([(1, 0)] * (n_seqs - 1), [0, 0]),
                                'p': 0}]

                # sum of K_hat, average of a_hat, d_min, d_max
                sim_data['K_hat'][K_idx, p_idx, idx] = sum(
                    _len(r['segment']) for r in results)
                sim_data['a_hat'][K_idx, p_idx, idx] = sum(
                    (r['segment'][1][1] + r['segment'][1][0]) / 2
                    for r in results) / len(results)
                sim_data['d_min'][K_idx, p_idx, idx] = sum(
                    d_range[0] for r in results
                    for d_range in r['segment'][0]) \
                    / (len(results) * (n_seqs - 1))
                sim_data['d_max'][K_idx, p_idx, idx] = sum(
                    d_range[1] for r in results
                    for d_range in r['segment'][0]) \
                    / (len(results) * (n_seqs - 1))
        sys.stderr.write('\n')

    return sim_data


def plot_simulated_K_p(sim_data, select_Ks, select_ps, suffix=''):
    Ks, ps = sim_data['Ks'], sim_data['ps']
    n_seqs = sim_data['n_seqs']
    scores = sim_data['scores']

    kw = {'marker': 'o', 'markersize': 3, 'alpha': .6, 'lw': 1.5}
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
                ax.plot(Ks, scores[mode][case][:, p_idx, :].mean(axis=1),
                        color=color, ls=ls, label=label, **kw)
        ax.set_xscale('log')
        ax.set_ylabel('%s score' % mode, fontsize=10)
        ax.set_xlabel('similarity length', fontsize=10)
        ax.set_xticks(Ks)
        ax.set_xticklabels(Ks)
    ax_H0_K.legend(loc='lower right', fontsize=10)
    ax_H1_K.legend(loc='upper right', fontsize=10)

    for mode, ax in zip(['H0', 'H1'], [ax_H0_p, ax_H1_p]):
        colors = color_code(select_Ks)
        for K, color in zip(select_Ks, colors):
            K_idx = Ks.index(K)
            for case, ls in zip(['pos', 'neg'], ['-', '--']):
                label = '' if case == 'neg' else 'K = %d' % K
                ax.plot(ps, scores[mode][case][K_idx, :, :].mean(axis=1),
                        color=color, ls=ls, label=label, **kw)
        ax.set_ylabel('%s score' % mode)
        ax.set_xlabel('similarity match probability')
        ax.set_xticks(ps)
        ax.set_xticklabels(ps)
        ax.legend(loc='lower left')

    for ax in [ax_H0_K, ax_H0_p]:
        ax.set_yscale('symlog')

    savefig(fig_scores, 'simulations_multiple[scores]%s.png' % suffix)

    # ======================================
    fig_hats = plt.figure(figsize=(8, 7))
    ax_K = fig_hats.add_subplot(2, 2, 1)
    ax_p = fig_hats.add_subplot(2, 2, 3)
    ax_d = fig_hats.add_subplot(2, 2, 2)
    ax_a = fig_hats.add_subplot(2, 2, 4)

    K_hats = sim_data['K_hat']
    p_hats = sim_data['p_hat']
    p_hats *= (1 + p_hats) / 2

    # estimated Ks for select ps
    colors = color_code(select_ps)
    for p, color in zip(select_ps, colors):
        p_idx = ps.index(p)
        ax_K.plot(Ks, K_hats[:, p_idx, :].mean(axis=1), color=color,
                  label='p = %.2f' % p, **kw)
        ax_K.set_xscale('log')
        ax_K.set_yscale('log')
        ax_K.set_xticks(Ks)
        ax_K.set_xticklabels(Ks)
    ax_K.set_ylabel('estimated similarity length', fontsize=10)
    ax_K.set_xlabel('true similarity length', fontsize=10)
    ax_K.plot(Ks, Ks, **truth_kw)
    ax_K.legend(loc='lower right', fontsize=8)

    # estimated ps for select Ks
    colors = color_code(select_Ks)
    for K, color in zip(select_Ks, colors):
        K_idx = Ks.index(K)
        ax_p.plot(ps, p_hats[K_idx, :, :].mean(axis=1), color=color,
                  label='K = %d' % K, **kw)
        ax_p.set_xticks(ps)
        ax_p.set_xticklabels(ps)
    ax_p.set_ylabel('estimated match probability', fontsize=10)
    ax_p.set_xlabel('true match probability', fontsize=10)
    ax_p.plot(ps, ps, **truth_kw)
    ax_p.legend(loc='lower right', fontsize=8)

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
        ax_d.plot(Ks, d_ctrs.mean(axis=1), label=label, **kw)
        ax_a.plot(Ks, a_hats[:, p_idx, :].mean(axis=1), label=label, **kw)
        for ax in [ax_d, ax_a]:
            ax.set_xscale('log')
            ax.set_xticks(Ks)
            ax.set_xticklabels(Ks)
            ax.set_xlabel('similarity length')
        ax_a.set_yscale('log')

    # FIXME manually set range
    ax_d.set_ylim(-200, 200)
    ax_d.set_xlabel('similarity length')
    ax_d.set_ylabel('estimated diagonal positions')
    ax_d.plot(Ks, [0] * len(Ks), **truth_kw)
    ax_d.legend(loc='upper left', fontsize=8)

    ax_a.plot(Ks, [n_seqs * K for K in Ks], **truth_kw)
    ax_a.set_xlabel('similarity length')
    ax_a.set_ylabel('estimated antidiagonal position')
    ax_a.legend(loc='upper left', fontsize=8)

    savefig(fig_hats, 'simulations_multiple[estimates]%s.png' % suffix)

    # ======================
    # time to score seeds as a function of K (for select ps) and p (for select
    # Ks)
    fig_t = plt.figure(figsize=(4, 8))
    ax_t_K = fig_t.add_subplot(2, 1, 1)
    ax_t_p = fig_t.add_subplot(2, 1, 2)

    kw = {'alpha': .5, 'lw': 1, 'marker': 'o', 'markersize': 2}

    # times as a function of p for select Ks
    colors = color_code(select_Ks)
    for K, color in zip(select_Ks, colors):
        K_idx = Ks.index(K)
        label = 'K = %d (wordlen = %d)' % (K, int(np.ceil(np.log(K))))
        plot_with_sd(ax_t_p, ps, 1000 * sim_data['time']['pos'][K_idx, :, :],
                     axis=1, color=color, ls='-', label=label, **kw)
    ax_t_p.set_xticks(ps)
    ax_t_p.set_xticklabels(ps)
    # NOTE make room for legend manually
    ax_t_p.set_ylim(None, 1400)
    ax_t_p.set_ylabel('time (ms) to score all seeds')
    ax_t_p.set_xlabel('true match probability')
    ax_t_p.legend(loc='upper left', fontsize=6)

    # times as a function of K for select ps
    colors = color_code(select_ps)
    for p_match, color in zip(select_ps, colors):
        p_idx = ps.index(p_match)
        label = 'p = %.2f' % p_match
        plot_with_sd(ax_t_K, Ks, 1000 * sim_data['time']['pos'][:, p_idx, :],
                     axis=1, color=color, ls='-', label=label, **kw)
    ax_t_K.set_xticks(Ks)
    ax_t_K.set_xticklabels(Ks, rotation=90, fontsize=8)
    # NOTE make room for legend manually
    ax_t_K.set_ylim(None, 1400)
    ax_t_K.set_ylabel('time (ms) to score all seeds')
    ax_t_K.set_xlabel('similarity lengths')
    ax_t_K.legend(loc='upper left', fontsize=6)

    savefig(fig_t, 'simulations_multiple[times]%s.png' % suffix)


def exp_simulated_K_p():
    """Performance of Word-Blot and associated statistical scores for multiple
    local similarity search in *simulations*:

    * z-scores assigned to diagonal n-dimensional volumes with and without
      similarities under the corresponding limiting Normal distributions
      (:math:`H0, H1`), for varying similarity lengths and match probabilities.
    * Estimated coordinates, lengths and match probabilities of local
      similarities for varying similarity lengths and match probabilities. Same
      procedure is used as in pairwise comparison simulations (i.e. trials
      containing a similarity consists of multiple sequences of length
      :math:`2K` whose substrings at positions :math:`[\\frac{K}{2},
      \\frac{3K}{2}]` are homologies of length :math:`K` and match probability
      :math:`p` in the multiple sequence sense).
    * Word lengths are automatically adjusted to ensure effectively linear
      growth of complexity with respect to sequence lengths.

    **Supported Claims**

    * z-scores calculated by
      :func:`biseqt.blot.WordBlotMultiuple.score_num_seeds` against the
      limiting normal distributions for multiple sequences are stable and
      comparable across different similarity lengths and match probabilities;
      thus both scores are reliable statistics.
    * Local similarities found by Word-Blot as per
      :func:`biseqt.blot.WordBlotMultiple.similar_segments` are accurate in
      length, estimated match probability, and coordinates.
    * Adjusting kmer word lengths in logarithmic relationship with sequence
      length maintains effectively linear complexity and up to sequences of
      length 6.4 kbp behaves well.

    .. figure::
        https://www.dropbox.com/s/un7v14rnng0gxru/
        simulations_multiple%5Bscores%5D.png?raw=1
       :target:
        https://www.dropbox.com/s/un7v14rnng0gxru/
        simulations_multiple%5Bscores%5D.png?raw=1
       :alt: lightbox

       Z-scores of the number of seeds in diagonal volumes for related (solid
       lines) and unrelated (dashed lines) multiple sequence (N=6) segments of
       varying lengths, calculated against the limiting distribution for
       unrelated pairs of sequences (*left*) and related pairs (*right*) as a
       function of similarity length (*top*) and match probability (*bottom*),
       n=50 samples. Note that, as desired, H0 score is stable for unrelated
       sequences of any length and H1 score is stable for related sequences of
       any length or match probability.

    .. figure::
        https://www.dropbox.com/s/owyhm9ks1rt5v2c/
        simulations_multiple%5Bestimates%5D.png?raw=1
       :target:
        https://www.dropbox.com/s/owyhm9ks1rt5v2c/
        simulations_multiple%5Bestimates%5D.png?raw=1
       :alt: lightbox

       Estimated length (*top left*), match probability (*bottom left*),
       diagonal position (*top right*) and antidiagonal position (*bottom
       right*) of multiple sequence (N=6) local similarities, n=50 samples,
       diagonal positions are averaged over N-1 values. Dashed black lines are
       ground truth for comparison. Length estimates are the sum of the
       lengths of all reported similar segments. Diagonal and antidiagonal
       coordinates are averaged over all reported segments. Match probability
       is taken from the longest similar segment.

    .. figure::
        https://www.dropbox.com/s/ticgyshozt6epi0/
        simulations_multiple%5Btimes%5D.png?raw=1
       :target:
        https://www.dropbox.com/s/ticgyshozt6epi0/
        simulations_multiple%5Btimes%5D.png?raw=1
       :alt: lightbox

       Time to find and score all multiple sequence (N=6) seeds as a function
       of varying similarity length (*left*) and varying match probability
       (*right*), n=20 samples, shaded regions indicated one standard
       deviation. Word lengths are logarithmically adjusted to sequence length.
    """
    Ks = [100 * 2 ** i for i in range(1, 7)]
    select_Ks = Ks[0], Ks[2], Ks[4]

    ps = [round(1 - .06 * i, 2) for i in range(1, 8)]
    select_ps = ps[0], ps[3], ps[5]

    n_samples = 20
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
    exp_biological_multiple_sequences()
    exp_rearrangement()
    exp_simulated_K_p()
