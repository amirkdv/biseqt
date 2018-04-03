import numpy as np
import sys
import matplotlib.gridspec as gridspec
import matplotlib
import logging
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from biseqt.pw import Aligner, BANDED_MODE, B_LOCAL
from biseqt.blot import WordBlot, WordBlotLocalRef
from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess
from util import plot_scored_seeds
from util import plot_similar_segment, adjust_pw_plot
from util import log, savefig
from util import fill_in_unknown
from util import load_fasta, with_dumpfile


@with_dumpfile
def sim_ig_genotyping(reads, genes, **kw):
    wordlens, p_min = kw['wordlens'], kw['p_min']
    mismatch_scores = kw['mismatch_scores']
    gap_open_score = kw['gap_open_score']
    gap_extend_score = kw['gap_extend_score']
    minlens = kw['minlens']

    A = Alphabet('ACGT')
    WB_kw = {'g_max': .5, 'sensitivity': .9, 'alphabet': A,
             'log_level': logging.WARN}

    from time import time
    sim_data = {
        'genes': genes,
        'wordlens': wordlens,
        'minlens': minlens,
        'reads': reads,
        'mappings': {
            read_name: {'V': None, 'D': None, 'J': None,
                        'D_start': None, 'D_end': None,
                        'time': 0}
            for read_name in reads
        },
        'p_min': p_min,
        'WB_kw': WB_kw,
        'mismatch_scores': mismatch_scores,
        'gap_open_score': gap_open_score,
        'gap_extend_score': gap_extend_score,
    }

    def _matches_in_seg(similar_segment):
        p_hat, seg = similar_segment['p'], similar_segment['segment']
        seglen = seg[1][1] - seg[1][0]
        return p_hat * seglen

    def _map_gene_type(read_, gene_type):
        t_start = time()
        assert gene_type in 'VJD'
        WB_kw['wordlen'] = wordlens[gene_type]
        WB = WordBlotLocalRef(read_, **WB_kw)
        candidates = {}
        for idx, (gene, gene_rec) in enumerate(genes[gene_type].items()):
            gene_len = len(gene_rec['seq'])

            K_min = gene_len / 2
            K_min = minlens[gene_type]
            similarities = list(
                WB.similar_segments(gene_rec['seq'], K_min, p_min)
            )
            if not similarities:
                continue
            res = max(similarities, key=lambda rec: _matches_in_seg(rec))
            candidates[gene] = res
        if not candidates:
            sys.stderr.write(' no Ig %s gene found!\n' % gene_type)
            return None, time() - t_start

        chosen_genes = dict(
            sorted(candidates.items(),
                   key=lambda rec: -_matches_in_seg(rec[1]))[:3]
        )
        return chosen_genes, time() - t_start

    def _add_aln(read_, gene_seq, gene_type, rec):
        d_band = rec['segment'][0]
        aligner_kw = {
            'match_score': 1,
            'mismatch_score': mismatch_scores[gene_type],
            'ge_score': gap_extend_score,
            'go_score': gap_open_score,
            'alnmode': BANDED_MODE,
            'alntype': B_LOCAL,
            'diag_range': (int(d_band[0]), int(d_band[1])),
        }

        with Aligner(read_, gene_seq, **aligner_kw) as aligner:
            aligner.solve()
            alignment = aligner.traceback()
            tx = alignment.transcript
            len_on_gene = sum(tx.count(op) for op in 'MSI')
            num_matches = tx.count('M')
            p_aln = 1. * num_matches / len_on_gene
            rec['p_aln'] = round(p_aln, 2)
            rec['len_aln'] = len_on_gene
            rec['alignment'] = alignment

        return rec

    for read_idx, (read_name, read_rec) in enumerate(reads.items()):
        read_seq = read_rec['seq']
        sys.stderr.write('%d/%d %s\n' % (read_idx + 1, len(reads), read_name))
        start_pos = 0
        for gene_type in 'VDJ':
            mapped_genes, t_elapsed = _map_gene_type(read_seq[start_pos:],
                                                     gene_type)
            sim_data['mappings'][read_name]['time'] += t_elapsed
            igblast = read_rec['igblast'][gene_type]
            if not igblast:
                # if igblast didn't map a read to any genes don't bother;
                # we don't have a ground truth.
                continue
            min_igblast_p = round(min(rec['p'] for rec in igblast.values()), 2)
            min_igblast_L = min(rec['length'] for rec in igblast.values())
            min_igblast_m = min(rec['num_matches'] for rec in igblast.values())
            sim_data['mappings'][read_name][gene_type] = {
                'min_igblast_p': min_igblast_p,
                'min_igblast_L': min_igblast_L,
                'min_igblast_m': min_igblast_m,
            }

            if not mapped_genes:
                continue
            for gene in mapped_genes:
                # run overlap NW for all chosen genes
                gene_seq = genes[gene_type][gene]['seq']
                # print gene
                mapped_genes[gene] = _add_aln(read_seq[start_pos:], gene_seq,
                                              gene_type, mapped_genes[gene])
                aln = mapped_genes[gene]['alignment']
                mapped_genes[gene]['start_pos'] = start_pos
                mapped_genes[gene]['end_pos'] = start_pos \
                    + aln.origin_start \
                    + aln.projected_len(aln.transcript, on='origin')

                true_pos = gene in reads[read_name]['igblast'][gene_type]

                ours_m = aln.transcript.count('M')
                p_aln = mapped_genes[gene]['p_aln']
                higher_p = p_aln >= min_igblast_p and \
                    ours_m >= min_igblast_m - 2
                higher_m = ours_m >= min_igblast_m and \
                    p_aln >= min_igblast_p - .01
                true_pos_forgiving = higher_p or higher_m
                sys.stderr.write('      %s: %s(%s) ' %
                                 (gene_type,
                                  '+' if true_pos else '-',
                                  '+' if true_pos_forgiving else '-'))
                sys.stderr.write(
                    '%s match=%d (%d),p=%.2f(%.2f)\n' %
                    (gene, ours_m, min_igblast_m, p_aln, min_igblast_p)
                )

            start_pos = min(rec['end_pos']
                            for gene, rec in mapped_genes.items())
            sim_data['mappings'][read_name][gene_type].update(mapped_genes)

        sys.stderr.write('      * %.2f s\n' %
                         (sim_data['mappings'][read_name]['time']))
    return sim_data


def plot_ig_genotyping(sim_data, suffix=''):
    reads = sim_data['reads']

    comparison = {
        'p': {'V': [], 'D': [], 'J': []},
        'K': {'V': [], 'D': [], 'J': []},
        'num_match': {'V': [], 'D': [], 'J': []},
    }
    accuracy = {
        'strict': {'V': 0, 'D': 0, 'J': 0},
        'forgiving': {'V': 0, 'D': 0, 'J': 0},
    }

    elapsed_times = []
    for read, mappings in sim_data['mappings'].items():
        elapsed_times.append(mappings['time'])

    for gene_type in 'VDJ':
        num_agreements_strict = 0
        num_agreements_forgiving = 0
        total_predictions = 0
        for read, mappings in sim_data['mappings'].items():
            if mappings[gene_type] is None:
                continue
            # pop these metrics so we can iterate over genes
            min_igblast_p = mappings[gene_type].pop('min_igblast_p')
            min_igblast_L = mappings[gene_type].pop('min_igblast_L')
            min_igblast_m = mappings[gene_type].pop('min_igblast_m')

            # HACK to see what happens if we pick the best performing match;
            # result on first 1000 reads: % 86.86 (% 99.80)
            # mappings[gene_type] = dict(sorted(
            #     mappings[gene_type].items(), key=lambda x: x[1]['p']
            # )[-1:])

            for gene, rec in mappings[gene_type].items():
                ours_m = rec['p_aln'] * rec['len_aln']
                total_predictions += 1

                # NOTE we're duplicating this logic to allow reevaluating
                # without redo-ing everything; eventually merge it into the
                # sim_* function. Note that we need to distinguish between
                # literal agreements between ours and igblast and those cases
                # where we argue our match has comparable quality to those of
                # igblast.
                color = 'r'
                if gene in reads[read]['igblast'][gene_type]:
                    num_agreements_strict += 1
                    num_agreements_forgiving += 1
                    color = 'g'
                else:
                    p_aln = rec['p_aln']
                    ours_m = rec['alignment'].transcript.count('M')
                    higher_p = p_aln >= min_igblast_p and \
                        ours_m >= min_igblast_m - 2
                    higher_m = ours_m >= min_igblast_m and \
                        p_aln >= min_igblast_p - .01
                    if higher_p or higher_m:
                        num_agreements_forgiving += 1
                comparison['p'][gene_type].append(
                    (rec['p_aln'], min_igblast_p, color)
                )
                comparison['K'][gene_type].append(
                    (rec['len_aln'], min_igblast_L, color)
                )
                comparison['num_match'][gene_type].append(
                    (ours_m, min_igblast_m, color)
                )

        accuracy['strict'][gene_type] = \
            100. * num_agreements_strict / total_predictions
        accuracy['forgiving'][gene_type] = \
            100. * num_agreements_forgiving / total_predictions
        sys.stderr.write('%s %.2f %.2f' % (gene_type,
                                           accuracy['strict'][gene_type],
                                           accuracy['forgiving'][gene_type]))

    avg_time = sum(elapsed_times) / len(sim_data['mappings'])
    print 't', avg_time

    def _extract_with_noise(mode, gene_type_):
        mag = {'p': .005, 'K': .5, 'num_match': .5}[mode]
        xs = [rec[0] for rec in comparison[mode][gene_type_]]
        xs += np.random.randn(len(comparison[mode][gene_type_])) * mag
        ys = [rec[1] for rec in comparison[mode][gene_type_]]
        ys += np.random.randn(len(comparison[mode][gene_type_])) * mag
        colors = [rec[2] for rec in comparison[mode][gene_type_]]
        return xs, ys, colors

    # probability of ours vs igblast
    fig = plt.figure(figsize=(14, 5))
    ax_V = fig.add_subplot(1, 3, 1)
    ax_D = fig.add_subplot(1, 3, 2)
    ax_J = fig.add_subplot(1, 3, 3)
    for gene_type, ax in zip('VDJ', [ax_V, ax_D, ax_J]):
        xs, ys, colors = _extract_with_noise('p', gene_type)
        ax.scatter(xs, ys, c=colors, alpha=.6, s=20, lw=0)
        ax.set_title('%s genes: \\%%%.2f (\\%%%.2f)' %
                     (gene_type,
                      accuracy['strict'][gene_type],
                      accuracy['forgiving'][gene_type]))
        ax.set_xlabel('WordBlot similarity')
        ax.set_ylabel('IgBlast similarity')
        # NOTE this can mislead if output is not checked already
        ax.set_xlim(.7, 1.05)
        ax.set_ylim(.7, 1.05)
        ax.set_aspect('equal')
        ax.plot([0, 1], [0, 1], c='k', lw=3, ls='--', alpha=.2)
    fig.suptitle('time per read: %.2f s' % avg_time, fontsize=8)
    savefig(fig, 'ig_genotyping[p-hat]%s.png' % suffix)

    # aligned length of ours vs igblast
    fig = plt.figure(figsize=(14, 5))
    ax_V = fig.add_subplot(1, 3, 1)
    ax_D = fig.add_subplot(1, 3, 2)
    ax_J = fig.add_subplot(1, 3, 3)
    for gene_type, ax in zip('VDJ', [ax_V, ax_D, ax_J]):
        xs, ys, colors = _extract_with_noise('K', gene_type)
        ax.scatter(xs, ys, c=colors, alpha=.4, s=20, lw=0)
        ax.set_title('%s genes: \\%%%.2f (\\%%%.2f)' %
                     (gene_type,
                      accuracy['strict'][gene_type],
                      accuracy['forgiving'][gene_type]))
        ax.set_xlabel('WordBlot aligned length')
        ax.set_ylabel('IgBlast aligned length')
        ax.set_aspect('equal')
        x_range, y_range = ax.get_xlim(), ax.get_ylim()
        xy_range = (min(x_range[0], y_range[0]), max(x_range[1], y_range[1]))
        ax.set_xlim(*xy_range)
        ax.set_ylim(*xy_range)
        ax.plot(ax.get_xlim(), ax.get_ylim(), c='k', lw=3, ls='--', alpha=.2)
    fig.suptitle('time per read: %.2f s' % avg_time, fontsize=8)
    savefig(fig, 'ig_genotyping[K-hat]%s.png' % suffix)

    # nt agreement of ours vs igblast
    fig = plt.figure(figsize=(14, 5))
    ax_V = fig.add_subplot(1, 3, 1)
    ax_D = fig.add_subplot(1, 3, 2)
    ax_J = fig.add_subplot(1, 3, 3)
    x_range, y_range = None, None
    for gene_type, ax in zip('VDJ', [ax_V, ax_D, ax_J]):
        xs, ys, colors = _extract_with_noise('num_match', gene_type)
        ax.scatter(xs, ys, c=colors, alpha=.4, s=20, lw=0)
        ax.set_title('%s genes: \\%%%.2f (\\%%%.2f)' %
                     (gene_type,
                      accuracy['strict'][gene_type],
                      accuracy['forgiving'][gene_type]))
        ax.set_xlabel('WordBlot matched nucleotides')
        ax.set_ylabel('IgBlast matched nucleotides')
        ax.set_aspect('equal')
        x_range, y_range = ax.get_xlim(), ax.get_ylim()
        xy_range = (min(x_range[0], y_range[0]), max(x_range[1], y_range[1]))
        ax.set_xlim(*xy_range)
        ax.set_ylim(*xy_range)
        ax.plot(ax.get_xlim(), ax.get_ylim(), c='k', lw=3, ls='--', alpha=.2)
    fig.suptitle('time per read: %.2f s' % avg_time, fontsize=8)
    savefig(fig, 'ig_genotyping[matches]%s.png' % suffix)


def exp_ig_genotyping():
    p_min = .8
    wordlens = {'J': 6, 'V': 8, 'D': 4}
    # cf. https://ncbiinsights.ncbi.nlm.nih.gov/tag/igblast/
    # https://www.ncbi.nlm.nih.gov/books/NBK279684/
    # note: I'm forcing these scores on blast as well:
    mismatch_scores = {'V': -1, 'D': -3, 'J': -2}
    minlens = {'V': 100, 'D': 10, 'J': 10}
    gap_open_score = -5
    gap_extend_score = -2
    suffix = '_first_1000'
    dumpfile = 'igh_s22%s.txt' % suffix

    A = Alphabet('ACGT')

    reads_file = 'data/igh-s22/s22%s.fa' % suffix
    igblast_file = 'data/igh-s22/igblast%s_clean.out' % suffix

    log('loading reads')
    reads = {}
    with open(reads_file) as f:
        for raw_seq, name, _ in load_fasta(f, num_seqs=-1):
            read = A.parse(fill_in_unknown(raw_seq, A))
            reads[name] = {
                'seq': read,
                'igblast': {'V': {}, 'D': {}, 'J': {}},
            }

    with open(igblast_file) as f:
        for line in f.readlines():
            rec = dict(zip(['gene_type', 'read', 'gene', 'p', 'length'],
                           line.strip().split()))
            gene, gene_type, name = rec['gene'], rec['gene_type'], rec['read']
            num_m, length = rec['length'].split('/')
            reads[name]['igblast'][gene_type][gene] = {
                'p': float(rec['p']) / 100,
                'length': int(length),
                'num_matches': int(num_m),
            }

    genes = {'V': {}, 'D': {}, 'J': {}}
    repertoire_prefix = 'data/igh-s22/imgt/'
    for key in genes:
        with open(repertoire_prefix + key + '.fa') as f:
            for raw_seq, name, _ in load_fasta(f):
                seq = A.parse(fill_in_unknown(raw_seq.upper(), A))
                genes[key][name] = {
                    'seq': seq,
                }
    log('running experiment')
    sim_data = sim_ig_genotyping(
        reads, genes, wordlens=wordlens, p_min=p_min, minlens=minlens,
        mismatch_scores=mismatch_scores,
        gap_open_score=gap_open_score,
        gap_extend_score=gap_extend_score,
        dumpfile=dumpfile
    )
    plot_ig_genotyping(sim_data, suffix)


def exp_rearrangement():
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
    exp_ig_genotyping()
