#!/usr/env/bin python
import numpy as np
from util import log, savefig, with_dumpfile
from biseqt.sequence import Alphabet, Sequence
from biseqt.kmers import KmerIndex
from matplotlib import pyplot as plt


# choose maxlen according to maximum wordlen of interest, so that under
# uniform assumption we get at least one copy of each kmer, i.e maxlen/4^w > 1
def load_seq(path, maxlen=1e7):
    log('loading sequence at %s' % path)
    with open(path) as f:
        seq = ''.join(s.strip().upper().replace('N', '')
                      for s in f.readlines() if s[0] != '>')
        contents = [{'A': 0, 'C': 1, 'G': 2, 'T': 3}[l] for l in seq[:maxlen]]
        return Sequence(Alphabet('ACGT'), contents=contents)


@with_dumpfile
def count_kmers(source_path, ws, **kw):
    A = Alphabet('ACGT')
    seq = load_seq(source_path)
    sim_data = {
        'counts': {w: [] for w in ws},
        'ws': ws,
        'len': len(seq)
    }
    for w in ws:
        path = 'dumpfiles/w=%d-kmer-freqs.db' % w
        kmer_index = KmerIndex(alphabet=A, wordlen=w, path=path)
        kmer_index.index_kmers(seq)
        with kmer_index.connection() as conn:
            cursor = conn.cursor()
            q = 'SELECT COUNT(*), kmer FROM %s GROUP BY kmer' % \
                kmer_index.kmers_table
            cursor.execute(q)
            sim_data['counts'][w] = list(row[0] for row in cursor)
    return sim_data


# the ns are indices for kmers (to some extent arbitrary, but we are sorting by
# probability to get a monotonic distribution).
# Pheno model: Gamma Distribution with alpha = 1-1/4 and beta = p0 = 1/4^w
def pheno_kmer_probs_neg_log(ns, w):
    p0 = .25 ** w
    return w + p0 * ns + np.log(ns) / 4


def plot_count_kmers(sim_data):
    N = sim_data['len']
    ws = sim_data['ws']
    fig_spectra = plt.figure(figsize=(6 * len(ws), 5))
    fig_prob = plt.figure(figsize=(6 * len(ws), 5))
    for idx, w in enumerate(ws):
        log('plotting kmer spectrum for w = %d' % w)
        ax_spectra = fig_spectra.add_subplot(1, len(ws), idx + 1)
        ax_prob = fig_prob.add_subplot(1, len(ws), idx + 1)

        p = .25 ** w
        counts = sorted(sim_data['counts'][w], reverse=True)
        probs = [1. * c / N for c in counts]
        ns = np.arange(1, len(probs) + 1)
        ax_prob.plot(ns, -np.log(probs), color='k', alpha=.8,
                     label='empirical')
        ax_prob.plot(ns, pheno_kmer_probs_neg_log(ns, w), color='b', ls='--',
                     lw=2, alpha=.7, label='$\Gamma$ approximation')
        ax_prob.axhline(y=-np.log(p), color='k', lw=2, alpha=.2, ls='--',
                        label='uniform')
        ax_prob.legend(loc='upper right', fontsize=10)

        # max count to consider (for clarity of spectrum plots)
        max_count = int(np.ceil(4 * p * N))
        counts = [1. * c for c in counts if c <= max_count]
        ax_spectra.hist(counts, bins=max_count, color='k', lw=0, normed=True,
                        alpha=.5)
        ax_spectra.axvline(x=p * N, color='k', lw=2, alpha=.5, ls='--')

        ax_spectra.set_xlabel('kmer abundance')
        ax_spectra.set_ylabel('frequency')
        ax_prob.set_xlabel('kmer')
        ax_prob.set_ylabel('- log(probability)')
        ax_prob.set_xlim(-.1 * len(probs), None)
        for ax in [ax_prob, ax_spectra]:
            ax.grid(True)
            ax.set_title(label='w=%d' % w, fontsize=10)
            ax.ticklabel_format(style='sci', scilimits=(-2, 2))

        # NOTE manually chosen xlims for w=9,11
        if w == 9:
            ax_prob.set_xlim(-2e3, 2e4)
        elif w == 11:
            ax_prob.set_xlim(-5e3, 5e4)

    savefig(fig_spectra, 'kmer_spectra.png')
    savefig(fig_prob, 'kmer_probs.png')


def exp_kmer_freqs():
    dumpfile = 'kmer_freq.txt'
    source_path = 'data/hg38/chr1.fa'
    ws = [7, 9, 11]
    sim_data = count_kmers(source_path, ws, dumpfile=dumpfile)
    plot_count_kmers(sim_data)


if __name__ == '__main__':
    exp_kmer_freqs()
