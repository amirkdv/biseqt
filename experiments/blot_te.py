#!/usr/bin/env python
import sys
import os
import logging
import numpy as np
import subprocess
from time import time
from matplotlib import pyplot as plt

from biseqt.blot import WordBlotLocalRef
from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess

from util import plot_with_sd, savefig
from util import with_dumpfile, log, DATA_DIR, load_fasta


def gen_data_set(**kw):
    TE_file = 'TE/TE.fa'
    genome_file = 'TE/genome[%.2f].fa'
    n_TE, n_seqs, = kw['n_TE'], kw['n_seqs']
    len_TE, n_TE_per_seq = kw['len_TE'], kw['n_TE_per_seq']
    ps = kw['ps']
    A = Alphabet('ACGT')
    TEs = [rand_seq(A, len_TE) for _ in range(n_TE)]

    with open(os.path.join(DATA_DIR, TE_file), 'w') as f:
        for idx, seq in enumerate(TEs):
            f.write('> %s\n%s\n' % (str(seq.content_id[:8]), seq))
        log('wrote %d transposable elements (%dnt) to %s' %
            (n_TE, len_TE, os.path.relpath(f.name, os.getcwd())))

    def _rand_genome(mutation_process=None):
        def _junk(): return rand_seq(A, 200)

        genome = A.parse('')
        included_TEs = []
        for TE_idx in np.random.choice(n_TE, n_TE_per_seq, replace=False):
            included_TEs.append(TEs[TE_idx].content_id[:8])
            genome += _junk() + TEs[TE_idx]
        genome += _junk()
        return mutation_process.mutate(genome)[0], included_TEs

    for p_match in ps:
        # distribute p_match evenly over gap and subst
        subst = gap = 1 - np.sqrt(p_match)
        M = MutationProcess(A, subst_probs=subst, ge_prob=gap, go_prob=gap)
        genomes = [_rand_genome(M) for _ in range(n_seqs)]
        with open(os.path.join(DATA_DIR, genome_file % p_match), 'w') as f:
            for genome, included_TEs in genomes:
                f.write('> %s\n%s\n' %
                        ('+'.join(TE_id for TE_id in included_TEs), genome))
            log('wrote %d genomes (%dnt) with %d TEs each (%.2f sim) to %s' %
                (n_seqs, sum(len(x[0]) for x in genomes) / len(genomes),
                 n_TE_per_seq, p_match, os.path.relpath(f.name, os.getcwd())))


@with_dumpfile
def sim_transposable_elements(**kw):
    TE_file = os.path.join(DATA_DIR, 'TE/TE.fa')
    genome_file = os.path.join(DATA_DIR, 'TE/genome[%.2f].fa')
    blast_db = os.path.join(DATA_DIR, 'TE/blast/te.db')
    n_seqs, ps, wordlen = kw['n_seqs'], kw['ps'], kw['wordlen']
    K_min, p_min = kw['K_min'], kw['p_min']
    A = Alphabet('ACGT')
    WB_kw = {'g_max': .2, 'sensitivity': .9, 'alphabet': A, 'wordlen': wordlen,
             'log_level': logging.WARN}

    kw = {'stdout': subprocess.PIPE, 'stderr': subprocess.PIPE}
    args = ['makeblastdb', '-dbtype', 'nucl', '-in', TE_file, '-out', blast_db]
    try:
        proc = subprocess.Popen(args, **kw)
        out, err = proc.communicate()
        assert proc.returncode == 0
    except (OSError, KeyError, AssertionError) as e:
        log('failed on running: ' + ' '.join(args))
        raise e

    with open(TE_file) as f:
        TEs = {name: A.parse(seq) for seq, name, _ in load_fasta(f)}

    def _zero(): return {key: np.zeros((len(ps), n_seqs))
                         for key in ['tp', 'fp', 'tn', 'fn']}

    def _seglen(seg): return (seg[1][1] - seg[1][0]) / 2.

    def _te_calls_stats(calls):
        stats = {}
        for genome, calls in calls.items():
            trues = genome.split('+')
            tp = sum(1 for true in trues if true in calls)
            fn = sum(1 for true in trues if true not in calls)
            fp = sum(1 for call in calls if call not in trues)
            tn = sum(1 for TE in TEs if TE not in calls and TE not in trues)
            stats[genome] = {'tp': tp, 'tn': tn, 'fp': fp, 'fn': fn}
        return stats

    sim_data = {
        'K_min': K_min,
        'p_min': p_min,
        'n_seqs': n_seqs,
        'ps': ps,
        'WB_kw': WB_kw,
        'te_calls': [{'blast': {}, 'wordblot': {}, 'nhmmer': {}} for _ in ps],
        'stats': {'blast': _zero(), 'wordblot': _zero(), 'nhmmer': _zero()},
        'times': {'blast': np.zeros(len(ps)),
                  'wordblot': np.zeros(len(ps)),
                  'nhmmer': np.zeros(len(ps))},
    }
    # ps = [.9]
    for p_idx, p_match in enumerate(ps):
        log('p = %.2f' % p_match)
        with open(genome_file % p_match) as f:
            genomes = {name: A.parse(seq) for seq, name, _ in load_fasta(f)}
            assert len(genomes) == n_seqs
        for alg in ['blast', 'wordblot', 'nhmmer']:
            sim_data['te_calls'][p_idx][alg] = {genome: set()
                                                for genome in genomes}

        # =============== NHMMER =================
        kw = {'stdout': subprocess.PIPE, 'stderr': subprocess.PIPE}
        # obtained manually:
        cols = ['target name', 'accession', 'query name', 'accession',
                'hmmfrom', 'hmmto', 'alifrom', 'alito', 'envfrom', 'envto',
                'sqlen', 'strand', 'E-value', 'score', 'bias',
                'description of target']
        args = ['nhmmer',
                '--qformat', 'fasta',
                '--tblout', '/dev/stdout',
                '-o', '/dev/null',
                TE_file, genome_file % p_match]
        t_nhmmer = time()
        try:
            proc = subprocess.Popen(args, **kw)
            out, err = proc.communicate()
            assert proc.returncode == 0
        except (OSError, KeyError, AssertionError) as e:
            log('failed on running: ' + ' '.join(args))
            raise e

        for line_no, rec in enumerate(out.strip().split('\n')):
            if rec[0] == '#':
                continue
            tokens = rec.split()
            assert len(tokens) == len(cols)
            TE_name = tokens[cols.index('query name')]
            genome = tokens[cols.index('target name')]
            to_pos = int(tokens[cols.index('hmmto')])
            from_pos = int(tokens[cols.index('hmmfrom')])

            if to_pos - from_pos < K_min:
                continue
            sim_data['te_calls'][p_idx]['nhmmer'][genome] |= set([TE_name])
        sim_data['times']['nhmmer'][p_idx] = time() - t_nhmmer
        call_stats = _te_calls_stats(sim_data['te_calls'][p_idx]['nhmmer'])
        for idx, (genome, counts) in enumerate(call_stats.items()):
            for key in ['tp', 'fp', 'tn', 'fn']:
                sim_data['stats']['nhmmer'][key][p_idx, idx] = counts[key]
        tp = sim_data['stats']['nhmmer']['tp'][p_idx, :]
        fp = sim_data['stats']['nhmmer']['fp'][p_idx, :]
        tn = sim_data['stats']['nhmmer']['tn'][p_idx, :]
        fn = sim_data['stats']['nhmmer']['fn'][p_idx, :]
        tpr = 1. * tp / (tp + fn)
        fpr = 1. * fp / (tn + fp)
        log('   [nhmmer] tpr = %.2f, fpr = %.2f' % (tpr.mean(), fpr.mean()))

        # ============= BLAST ==================
        kw = {'stdout': subprocess.PIPE, 'stderr': subprocess.PIPE}
        # args = ['blastn',
        # rmblast from
        #   -> ftp://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/2.2.28
        args = ['ncbi-rmblastn-2.2.28/bin/rmblastn',
                '-gapextend', '2',  # necessary for rmblast to not crash!
                '-db', blast_db, '-word_size', str(wordlen),
                '-perc_identity', str(int(p_min * 100)),
                '-evalue', str(1e3),
                '-outfmt', '6 qseqid sseqid length pident',
                '-query', genome_file % p_match]
        t_blast = time()
        try:
            proc = subprocess.Popen(args, **kw)
            out, err = proc.communicate()
            assert proc.returncode == 0
        except (OSError, KeyError, AssertionError) as e:
            log('failed on running: ' + ' '.join(args))
            sys.stderr.write('STDOUT:\n  ' + '\n  '.join(out.split('\n')))
            sys.stderr.write('STDERR:\n  ' + '\n  '.join(err.split('\n')))
            raise e

        for rec in out.strip().split('\n'):
            genome, TE_name, length, p_hat = rec.split()
            # we already asking blastn to filter by percent identity, just
            # check length:
            if int(length) < K_min:
                continue
            sim_data['te_calls'][p_idx]['blast'][genome] |= set([TE_name])
        sim_data['times']['blast'][p_idx] = time() - t_blast
        call_stats = _te_calls_stats(sim_data['te_calls'][p_idx]['blast'])
        for idx, (genome, counts) in enumerate(call_stats.items()):
            for key in ['tp', 'fp', 'tn', 'fn']:
                sim_data['stats']['blast'][key][p_idx, idx] = counts[key]
        tp = sim_data['stats']['blast']['tp'][p_idx, :]
        fp = sim_data['stats']['blast']['fp'][p_idx, :]
        tn = sim_data['stats']['blast']['tn'][p_idx, :]
        fn = sim_data['stats']['blast']['fn'][p_idx, :]
        tpr = 1. * tp / (tp + fn)
        fpr = 1. * fp / (tn + fp)
        log('   [blast] tpr = %.2f, fpr = %.2f' % (tpr.mean(), fpr.mean()))

        # =============== WORD-BLOT ===============
        t_wordblot = time()
        for seq_idx, (genome, seq) in enumerate(genomes.items()):
            WB = WordBlotLocalRef(seq, **WB_kw)
            for TE_name, TE in TEs.items():
                # Word blot output is guaranteed to satisfy minimum length
                # and probability requirements by construction
                res = list(WB.similar_segments(TE, K_min, p_min))
                if res:
                    sim_data['te_calls'][p_idx]['wordblot'][genome] |= \
                        set([TE_name])

        sim_data['times']['wordblot'][p_idx] = time() - t_wordblot

        call_stats = _te_calls_stats(sim_data['te_calls'][p_idx]['wordblot'])
        for idx, (genome, counts) in enumerate(call_stats.items()):
            for key in ['tp', 'fp', 'tn', 'fn']:
                sim_data['stats']['wordblot'][key][p_idx, idx] = counts[key]
        tp = sim_data['stats']['wordblot']['tp'][p_idx, :]
        fp = sim_data['stats']['wordblot']['fp'][p_idx, :]
        tn = sim_data['stats']['wordblot']['tn'][p_idx, :]
        fn = sim_data['stats']['wordblot']['fn'][p_idx, :]
        tpr = 1. * tp / (tp + fn)
        fpr = 1. * fp / (tn + fp)
        log('   [wordblot] tpr = %.2f, fpr = %.2f' % (tpr.mean(), fpr.mean()))

    return sim_data


def plot_transposable_elements(sim_data, suffix=''):
    ps = sim_data['ps']
    fig = plt.figure(figsize=(12, 4))
    ax_tpr = fig.add_subplot(1, 3, 1)
    ax_fpr = fig.add_subplot(1, 3, 2)
    ax_time = fig.add_subplot(1, 3, 3)

    for key, color in zip(['blast', 'nhmmer', 'wordblot'], 'bgk'):
        tp = sim_data['stats'][key]['tp']
        fp = sim_data['stats'][key]['fp']
        tn = sim_data['stats'][key]['tn']
        fn = sim_data['stats'][key]['fn']
        tpr = tp / (tp + fn)
        fpr = fp / (tn + fp)

        label = 'rmblast' if key == 'blast' else key
        kw = {'alpha': .6, 'lw': 1, 'color': color, 'marker': 'o',
              'markersize': 3, 'label': label}
        plot_with_sd(ax_tpr, ps, tpr, axis=1, **kw)
        plot_with_sd(ax_fpr, ps, 1 - fpr, axis=1, **kw)
        ax_time.plot(ps, sim_data['times'][key], **kw)

    for ax in [ax_fpr, ax_tpr, ax_time]:
        ax.set_xlabel('match probability')
    for ax in [ax_fpr, ax_tpr]:
        ax.legend(loc='lower right')
    ax_time.legend(loc='upper left')
    ax_tpr.set_ylabel('Sensitivity (TPR)')
    ax_fpr.set_ylabel('Specificity (1 - FPR)')
    ax_time.set_ylabel('time (s)')
    ax_fpr.set_ylim(-.1, 1.2)
    ax_tpr.set_ylim(-.1, 1.2)

    savefig(fig, 'transposable_elements%s.png' % suffix)


def exp_transposable_elements():
    """Comparison between Word-Blot, Blast, and Hmmer in identifying
    *simulated* transposable elements (TE). TEs of length 1000nt are simulated
    (20 TEs), 20 genomes are simulated, each containing mutated versions of 10
    TEs with varying match probabilities. Each algorithm (Blast or Word-Blot)
    is used with word length 8, minimum percent identity %60, and minimum
    length 500nt to identify local similarities between all known TEs and the
    genomic sequences and a TE is called if a local similarity of above
    properties is found between a genome and a TE. Example blastn command is::

        $ makeblastdb -dbtype nucl -in TE.fa -out blast/te.db
        $ rmblast -gapextend 2 \\  # necessary for rmblast to not crash!
                 -db blast/te.db \\
                 -word_size 8 \\
                 -evalue 1e3 \\
                 -outfmt '6 qseqid sseqid length pident' \\
                 -query 'genome[0.69].fa'

    and example nhmmer command is::

        $ nhmmer --qformat fasta \\
                 --tblout /dev/stdout -o /dev/null \\
                 TE.fa 'genome[0.69.fa]'


    **Supported Claims**

    * Word-Blot greatly outperforms rmblast (blastn tuned for RepeatMasker) in
      sensitivity (with roughly 5x time) while maintaining reasonable
      specificity (>%98) specifically when identifying *distant* TE homologos
      sequences (between %40 and %70 identity).
    * Word-Blot outperforms nhmmer in speed (more than 2x) and to some extent
      in sensitivity (>%10) when identifying *distant* TEs.

    .. figure::
        https://www.dropbox.com/s/x4j1gbbpddpmz5k/
        transposable_elements.png?raw=1
       :target:
        https://www.dropbox.com/s/x4j1gbbpddpmz5k/
        transposable_elements.png?raw=1
       :alt: lightbox

       Sensitivity (*left*), specificity (*middle*) and computation time
       (*right*) for Blast (blue), nhmmer (green), and Word-Blot (black) for
       varying similarities (horizontal axes) between simulated TEs and their
       occurences in simulated genomes.
    """
    ps = [.3 + .03 * i for i in range(21)]
    n_seqs = 20
    # n_TE = 20
    # len_TE = 1000
    # n_TE_per_seq = 10
    # gen_data_set(n_TE=n_TE, n_seqs=n_seqs, len_TE=len_TE, ps=ps,
    #              n_TE_per_seq=n_TE_per_seq)

    K_min, p_min = 500, .6
    suffix = ''
    wordlen = 8  # kept in memory; don't go too high up
    dumpfile = 'transposable_elements%s.txt' % suffix

    sim_data = sim_transposable_elements(
        K_min=K_min, p_min=p_min, ps=ps, wordlen=wordlen, n_seqs=n_seqs,
        dumpfile=dumpfile, ignore_existing=False,
    )
    plot_transposable_elements(sim_data, suffix=suffix)


if __name__ == '__main__':
    exp_transposable_elements()
