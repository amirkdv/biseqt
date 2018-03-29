import os
from itertools import combinations
import logging
import numpy as np
from itertools import product
from matplotlib import pyplot as plt

from biseqt.util import ProgressIndicator
from biseqt.blot import WordBlotOverlap
from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess, rand_read

from util import plot_classifier, log, with_dumpfile, plot_cdf, savefig
from util import DATA_DIR
from util import plot_with_sd


# M is the mutation process to insert noise to each read, note that when
# comparing two reads later they have both suffered mutations at the rate
# dictated by M and hence are "twice" furhter apart as they are from the
# original genome.
def create_simulated_reads(M, path, gap=None, subst=None,
                           genome_size=10000, read_len=2000, num_reads=40):
    log('simulating %d sequencing reads of %d bp for genome of %d bp.' %
        (num_reads, read_len, genome_size))
    A = Alphabet('ACGT')
    M = MutationProcess(A, subst_probs=subst, ge_prob=gap, go_prob=gap)
    genome = rand_seq(A, genome_size)
    indic = ProgressIndicator(num_total=num_reads)
    indic.start()
    reads, mappings = [], []
    for read, start in rand_read(genome, len_mean=read_len, num=num_reads):
        indic.progress()
        reads.append(M.mutate(read)[0])
        mappings.append((start, start + len(read)))
    indic.finish()

    path = os.path.join(DATA_DIR, path)
    with open(path, 'w') as f:
        for m, read in zip(mappings, reads):
            f.write('> + %d:%d\n%s\n' % (m[0], m[1], str(read)))


def load_mapped_reads(path, max_num=100):
    A = Alphabet('ACGT')
    reads, mappings = [], []
    log('loading sequencing reads from %s' % path)
    count = 0
    reads, mappings = [], []
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                break
            rc, coordinates = line[2:].split()
            assert rc in '+-'
            m0, m1 = coordinates.split(':')
            mappings.append((rc, int(m0), int(m1)))
            reads.append(A.parse(f.readline().strip()))
            count += 1
            if count == max_num:
                break
    return reads, mappings


@with_dumpfile
def sim_overlap_simulations(**kw):
    A = Alphabet('ACGT')
    wordlen, gap, subst, Ks = kw['wordlen'], kw['gap'], kw['subst'], kw['Ks']
    n_samples = kw['n_samples']
    WB_kw = {
        'g_max': .6,
        'sensitivity': .99,
        'alphabet': A,
        'wordlen': wordlen,
        'log_level': logging.WARN,
    }
    p_match = (1 - gap) * (1 - subst)

    def _zero():
        return {key: np.zeros((len(Ks), n_samples))
                for key in ['d_true', 'p_true', 'p', 'score', 'd', 'r']}

    sim_data = {
        'results': {
            'overlap': _zero(),
            'noisy overlap': _zero(),
            'short overlap': _zero(),
            'unrelated': _zero(),
        },
        'Ks': Ks,
        'gap': gap,
        'subst': subst,
        'n_samples': n_samples,
        'WB_kw': WB_kw,
    }
    M = MutationProcess(A, subst_probs=subst, ge_prob=gap, go_prob=gap)
    M_noisy = MutationProcess(
        A, subst_probs=subst * 2, ge_prob=gap * 2, go_prob=gap * 2
    )
    for (K_idx, K), idx in product(enumerate(Ks), range(n_samples)):
        d_true = np.random.randint(max(Ks))
        log('K = %d (sample %d / %d)' % (K, idx + 1, n_samples))

        def _flank(length): return rand_seq(A, length)

        overlap_seq = rand_seq(A, K)
        seq_pairs = {
            'overlap': (
                _flank(d_true) + overlap_seq,
                M.mutate(overlap_seq)[0] + _flank(min(Ks))
            ),
            'noisy overlap': (
                _flank(d_true) + overlap_seq,
                M_noisy.mutate(overlap_seq)[0] + _flank(min(Ks))
            ),
            'short overlap': (
                _flank(K / 2 + d_true) + overlap_seq + _flank(K / 2),
                _flank(K / 2) + M.mutate(overlap_seq)[0] + _flank(K / 2)
            ),
            'unrelated': (rand_seq(A, K), rand_seq(A, K))
        }
        p_true = {
            'overlap': p_match,
            'noisy overlap': (1 - 2 * gap) * (1 - 2 * subst),
            'short overlap': None,
            'unrelated': None,
        }

        for key, (S, T) in seq_pairs.items():
            WB = WordBlotOverlap(S, T, **WB_kw)
            res = WB.highest_scoring_overlap_band2(p_match * .95)
            sim_data['results'][key]['d_true'][K_idx, idx] = d_true
            sim_data['results'][key]['p'][K_idx, idx] = res['p']
            sim_data['results'][key]['p_true'][K_idx, idx] = p_true[key]
            sim_data['results'][key]['d'][K_idx, idx] = sum(res['d_band']) / 2
            sim_data['results'][key]['score'][K_idx, idx] = res['score']
            sim_data['results'][key]['r'][K_idx, idx] = \
                (res['d_band'][1] - res['d_band'][0]) / 2
    return sim_data


def plot_overlap_simulations(sim_data):
    Ks = sim_data['Ks']
    colors = {
        'overlap': 'g',
        'unrelated': 'r',
        'noisy overlap': 'b',
        'short overlap': 'm'
    }

    fig = plt.figure(figsize=(14, 5))
    ax_p_hat = fig.add_subplot(1, 3, 1)
    ax_d_hat = fig.add_subplot(1, 3, 2)
    ax_cdf = fig.add_subplot(1, 3, 3)

    for key, res in sim_data['results'].items():
        color = colors[key]
        plot_with_sd(ax_p_hat, Ks, res['p'], axis=1, alpha=.8, marker='o',
                     markersize=3, color=color, label=key.replace('_', ' '))

        if res['p_true'][0, 0] is not None:
            p_true = res['p_true'][0, 0]
            ax_p_hat.plot([min(Ks), max(Ks)], [p_true, p_true], lw=5,
                          alpha=.3, color=color, ls='--')

        if key != 'unrelated':
            color = colors[key]
            ax_d_hat.scatter(res['d_true'], res['d'], color=color, s=5, lw=0,
                             alpha=.6, label=key)

        plot_cdf(ax_cdf, res['p'].flatten(), color=color, smooth_radius=10,
                 label=key)

    for ax in [ax_p_hat, ax_d_hat, ax_cdf]:
        ax.legend(loc='best', fontsize=8)
    ax_cdf.set_title('Estimated match probability CDF')
    ax_cdf.set_xlabel('estimated match probability')

    ax_p_hat.set_ylim(-.1, 1.1)
    ax_p_hat.set_xlabel('overlap length')
    ax_p_hat.set_ylabel('estimated match probability')

    ax_d_hat.set_aspect('equal')
    ax_d_hat.plot([0, max(Ks)], [0, max(Ks)], lw=6, color='k', ls='--',
                  alpha=.2)
    ax_d_hat.set_xlabel('overlap diagonal position')
    ax_d_hat.set_ylabel('estimated overlap diagonal position')

    fig.tight_layout()
    savefig(fig, 'overlaps_simulations[d,p-hat].png')


def exp_overlap_simulations(**kw):
    wordlen = 6
    gap = .1
    subst = .1
    Ks = [i * 500 for i in range(1, 11)]

    n_samples = 20
    dumpfile = 'overlap_simulations.txt'

    sim_data = sim_overlap_simulations(
        wordlen=wordlen, gap=gap, subst=subst, Ks=Ks, n_samples=n_samples,
        dumpfile=dumpfile
    )
    plot_overlap_simulations(sim_data)


@with_dumpfile
def sim_sequencing_reads_overlap(**kw):
    reads_path, p_min, wordlen = kw['reads_path'], kw['wordlen'], kw['p_min']
    sim_data = {
        'overlap_band': {},
        'p_min': p_min,
        'wordlen': wordlen,
        'reads_path': reads_path,
    }
    A = Alphabet('ACGT')

    WB_kw = {
        'g_max': .6,
        'sensitivity': .9,
        'alphabet': A,
        'wordlen': wordlen,
        'path': ':memory:',
        'log_level': logging.WARN,
    }

    reads, mappings = load_mapped_reads(reads_path)

    num_reads = len(reads)
    assert len(reads) == len(mappings)
    num_total = (num_reads * (num_reads - 1)) / 2
    log('finding all overlapping pairs of reads')
    indic = ProgressIndicator(num_total=num_total, percentage=False)
    indic.start()
    for i, j in combinations(range(num_reads), 2):
        indic.progress()
        WB = WordBlotOverlap(reads[i], reads[j], **WB_kw)

        res = WB.highest_scoring_overlap_band2(p_min=p_min)
        sim_data['overlap_band'][(i, j)] = res

    indic.finish()
    return sim_data


def plot_sequencing_reads_overlap(sim_data, suffix=''):
    reads, mappings = load_mapped_reads(sim_data['reads_path'])

    def _overlaps(map0, map1):
        rc0, from0, to0 = map0
        rc1, from1, to1 = map1
        if rc0 != rc1:
            return False
        overlap_len = min(to0, to1) - max(from0, from1)
        if overlap_len > 10:
            return True
        else:
            return False

    pos, neg = [], []
    for i, j in combinations(range(len(reads)), 2):
        if sim_data['overlap_band'][(i, j)] is None:
            p_hat = 0
        else:
            p_hat = sim_data['overlap_band'][(i, j)]['p']

        if _overlaps(mappings[i], mappings[j]):
            pos.append(p_hat)
        else:
            neg.append(p_hat)

    labels = ['overlapping', 'non-overlapping']
    plot_classifier('overlaps[roc]%s.png' % suffix, pos, neg, labels=labels)


def exp_sequencing_reads_overlap():
    hf_gap = .15
    hf_subst = .15
    hf_match = (1 - hf_gap) * (1 - hf_subst)

    wordlen = 12
    suffix = '[w=12]'

    reads_path = 'data/leishmenia/reads_mapped.fa'
    dumpfile = 'overlaps%s.txt' % suffix

    sim_data = sim_sequencing_reads_overlap(
        reads_path=reads_path, wordlen=wordlen, p_min=hf_match,
        dumpfile=dumpfile
    )

    # HACK
    sim_data['reads_path'] = reads_path
    from util import pickle_dump
    pickle_dump(sim_data, dumpfile)

    plot_sequencing_reads_overlap(sim_data, suffix=suffix)


if __name__ == '__main__':
    # mut_gap = .08
    # mut_subst = .05
    # create_simulated_reads(M, 'reads_sim.fa', gap=mut_gap, subst=mut_subst)

    exp_overlap_simulations()
    exp_sequencing_reads_overlap()
