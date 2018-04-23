from itertools import combinations
import logging
import numpy as np
import sys
from matplotlib import pyplot as plt
from time import time

from biseqt.util import ProgressIndicator
from biseqt.blot import WordBlotOverlap, WordBlotOverlapRef
from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess

from util import plot_classifier, log, with_dumpfile, plot_cdf, savefig
from util import plot_with_sd


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
        'g_max': .3,
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
            'incomplete overlap': _zero(),
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
    for K_idx, K in enumerate(Ks):
        log('K = %d' % K, newline=False)
        for idx in range(n_samples):
            sys.stderr.write('.')
            d_true = np.random.randint(max(Ks))

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
                'incomplete overlap': (
                    _flank(K / 2 + d_true) + overlap_seq + _flank(K / 2),
                    _flank(K / 2) + M.mutate(overlap_seq)[0] + _flank(K / 2)
                ),
                'unrelated': (rand_seq(A, K), rand_seq(A, K))
            }
            p_true = {
                'overlap': p_match,
                'noisy overlap': (1 - 2 * gap) * (1 - 2 * subst),
                'incomplete overlap': None,
                'unrelated': None,
            }

            for key, (S, T) in seq_pairs.items():
                WB = WordBlotOverlap(S, T, **WB_kw)
                res = WB.highest_scoring_overlap_band(p_match * .95)
                sim_data['results'][key]['d_true'][K_idx, idx] = d_true
                sim_data['results'][key]['p'][K_idx, idx] = res['p']
                sim_data['results'][key]['p_true'][K_idx, idx] = p_true[key]
                sim_data['results'][key]['d'][K_idx, idx] = \
                    sum(res['d_band']) / 2
                sim_data['results'][key]['score'][K_idx, idx] = res['score']
                sim_data['results'][key]['r'][K_idx, idx] = \
                    (res['d_band'][1] - res['d_band'][0]) / 2
        sys.stderr.write('\n')
    return sim_data


def plot_overlap_simulations(sim_data):
    Ks = sim_data['Ks']
    colors = {
        'overlap': 'g',
        'unrelated': 'r',
        'noisy overlap': 'b',
        'incomplete overlap': 'm'
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
    ax_cdf.set_xlabel('estimated match probability')
    ax_cdf.set_ylabel('cumulative distribution')

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


def exp_overlap_simulations():
    """Performance of Word-Blot in overlap discovery *simulations*. In each
    trial for overlap length :math:`K`, a pair of sequences are presented to
    :class:`biseqt.blot.WordBlotOverlap` which have either of:

    * an overlap similarity of length :math:`K` with gap and substitution
      probabilities both equal to 0.1 (*overlap* pair),
    * an overlap similarity of length :math:`K` with twice as much mutation
      probabilities (*noisy overlap* pair),
    * a similarity of length :math:`K` with 0.1 substitution and gap
      probabilities, flanked by two unrelated regions of length
      :math:`\\frac{K}{2}` such that it does not qualify as a proper overlap
      similarity (*incomplete overlap* pair),
    * no similarities.

    In each case the highest *overlap band match probability* reported by
    Word-Blot is used as the discriminating statistic.

    **Supported Claims**

    * In simulations Word-Blot accurately detects overlaps, estimates their
      diagonal position, and their match probability and thus can be used for
      quality-sensitive overlap detection and as a guide for further banded
      alignment.

    .. figure::
        https://www.dropbox.com/s/91amhsoyep99204/
        overlaps_simulations%5Bd%2Cp-hat%5D.png?raw=1
       :target:
        https://www.dropbox.com/s/91amhsoyep99204/
        overlaps_simulations%5Bd%2Cp-hat%5D.png?raw=1
       :alt: lightbox

       Highest overlap band match probability reported by Word-Blot
       (*left*) for each of the four cases; word length 6, n=20, shaded regions
       indicate one standard deviation. For overlap and noisy overlap pairs the
       known truth is shown in thick dashed lines of corresponding color.
       Estimated diagonal position in each trial is shown in each of the three
       cases where similarities exist (*middle*), true diagonal position is the
       dashed grey diagonal. For each of case the cumulative probability
       distribution of the match probability reported by Word-Blot is shown
       (*right*), the clear separation of which indicates that Word-Blot is
       capable of detecting overlaps of a certain quality with high accuracy.
    """
    wordlen = 6
    gap = .1
    subst = .1
    Ks = [i * 500 for i in range(1, 11)]

    n_samples = 20
    dumpfile = 'overlap_simulations.txt'

    sim_data = sim_overlap_simulations(
        wordlen=wordlen, gap=gap, subst=subst, Ks=Ks, n_samples=n_samples,
        dumpfile=dumpfile, ignore_existing=False,
    )
    plot_overlap_simulations(sim_data)


@with_dumpfile
def sim_overlap_sequencing_reads(**kw):
    reads_path, wordlen, g_max = kw['reads_path'], kw['wordlen'], kw['g_max']
    sim_data = {
        'overlap_band': {},
        'wordlen': wordlen,
        'reads_path': reads_path,
        'avg_time': 0,
        'avg_len': 0,
    }
    A = Alphabet('ACGT')

    WB_kw = {
        'g_max': g_max,
        'sensitivity': .9,
        'alphabet': A,
        'wordlen': wordlen,
        'path': ':memory:',
        'log_level': logging.WARN,
    }

    reads, mappings = load_mapped_reads(reads_path)
    sim_data['avg_len'] = int(
        1. * sum(len(read) for read in reads) / len(reads)
    )

    num_reads = len(reads)
    assert len(reads) == len(mappings)
    num_total = (num_reads * (num_reads - 1)) / 2
    log('finding all overlapping pairs of reads')
    indic = ProgressIndicator(num_total=num_total, percentage=False)
    indic.start()
    t_start = time()
    for i in range(num_reads - 1):
        WB = WordBlotOverlapRef(reads[i], **WB_kw)
        for j in range(i + 1, num_reads):
            indic.progress()
            res = WB.highest_scoring_overlap_band(reads[j])
            sim_data['overlap_band'][(i, j)] = res
    sim_data['avg_time'] = (time() - t_start) / num_total
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
        if overlap_len > 50:
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
    fig, _ = plot_classifier(pos, neg, labels=labels, mark_threshold=.75)
    print 'avg_time', sim_data['avg_time'], 'avg_len', sim_data['avg_len']
    savefig(fig, 'overlaps[roc]%s.png' % suffix, comment='classifier plot')


def exp_overlap_sequencing_reads():
    """Performance of Word-Blot in overlap discovery on *Pac Bio sequencing*
    reads. A thousand Pac Bio reads from chromosome 1 of Leishmenia Donovani
    were mapped to a reference genome using BLASR which provides a ground truth
    for classification of reads into overlapping and non-overlapping pairs.
    Word-Blot overlap discovery is then applied to all pairs of reads and the
    performance of the highest overlap band match probability as a classifier
    is shown.

    **Supported Claims**

    * Word-Blot overlap discovery is a good classifier for
      overlapping/non-overlapping pairs among Pac Bio sequencing reads.

    .. figure::
        https://www.dropbox.com/s/8cmsnr5x0oe5bbt/
        overlaps%5Broc%5D%5Bw%3D10%5D.png?raw=1
       :target:
        https://www.dropbox.com/s/8cmsnr5x0oe5bbt/
        overlaps%5Broc%5D%5Bw%3D10%5D.png?raw=1
       :alt: lightbox

       ROC curve for Word-Blot (word length = 10) overlap band match
       probability estimates as a classifier statistic for
       overlapping/non-overlapping pairs in Pac Bio sequencing reads. Average
       comparison time for reads of average length 7kbp is 60 ms. With match
       probability threshold 0.75 (indicated in all subplots in blue) we
       obtains %99 specificity and %98 negative predictive value (the crucial
       factors in overlap detection) and %45 sensitivity and %46 positive
       predictive value.
    """
    wordlen = 9
    suffix = '[w=9]'

    reads_path = 'data/leishmenia/reads_mapped.fa'
    dumpfile = 'overlaps%s.txt' % suffix

    sim_data = sim_overlap_sequencing_reads(
        reads_path=reads_path, wordlen=wordlen, g_max=.2,
        dumpfile=dumpfile, ignore_existing=False
    )
    plot_sequencing_reads_overlap(sim_data, suffix=suffix)


if __name__ == '__main__':
    exp_overlap_simulations()
    exp_overlap_sequencing_reads()
