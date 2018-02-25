import os
import sys
import re
import pickle
import numpy as np
from bisect import bisect_left
from textwrap import TextWrapper
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.gridspec as gridspec
from Bio import AlignIO
from itertools import combinations
import pysam
from biseqt.stochastics import rand_seq


plt.rc('text', usetex=True)


def log(msg, newline=True):
    wrapper = TextWrapper(initial_indent='* ', subsequent_indent='  ',
                          width=100)
    sys.stderr.write('\n'.join(wrapper.wrap(msg)) + '\n' if newline else '')


def _make_absolute(path):
    if not os.path.isabs(path):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), path)
    return path


DATA_DIR = _make_absolute(os.getenv('DATA_DIR', 'data/'))
PLOT_DIR = _make_absolute(os.getenv('PLOT_DIR', 'plots/'))
DUMP_DIR = _make_absolute(os.getenv('DUMP_DIR', 'dumpfiles/'))


# =============================================================================
# Generic IO Helpers
# =============================================================================
def pickle_dump(obj, path, comment=None):
    """Pickle dumps the given object to given relative path in DUMP_DIR."""
    try:
        with open(os.path.join(DUMP_DIR, path), 'w') as f:
            pickle.dump(obj, f)
    except IOError, e:
        log('Error writing to %s: %s' % (path, str(e)))
        n = np.random.choice(list('0123456789abcdef'), size=10)
        path = os.join('/tmp', ''.join(str(i) for i in n))
        with open(path, 'w') as f:
            pickle.dump(obj, f)
    log('dumped %s to %s.' % (comment if comment else 'object', path))


def pickle_load(path):
    log('loading dumpfile %s' % path)
    path = os.path.join(DUMP_DIR, path)
    with open(path, 'rb') as f:
        return pickle.load(f)


def with_dumpfile(func):
    """Decorator for storing simulation results, usage:

        @with_dumpfile
        foo(..., **kw):
            ...
            return results # this is dumped to dumpfile

        # if dumpfile exists, foo() will not be executed
        res = foo(..., dumpfile='foo.txt')
        # force a re-run of foo(), dumpfile will be overwritten
        res = foo(..., dumpfile='foo.txt', ignore_existing=True)
    """
    def wrapped_func(*args, **kwargs):
        dumpfile = kwargs['dumpfile']
        path = os.path.join(DUMP_DIR, dumpfile)
        if os.path.exists(path) and not kwargs.get('ignore_existing', False):
            return pickle_load(dumpfile)
        out = func(*args, **kwargs)
        pickle_dump(out, dumpfile)
        return out

    return wrapped_func


# =============================================================================
# Real bio sequence and alignment helpers
# =============================================================================
# TODO this is really slow on large sequences
def load_fasta(f):
    """Given file handle reads fasta sequences and yields tuples of (seq, name,
       pos) containing the raw sequence, its FASTA name, and its starting
       position in f.
    """
    cur_name = cur_seq = ''
    cur_pos = 0

    # NOTE `for raw_line in f` uses a read-ahead buffer which makes `f.tell()`
    # useless for remembering where a sequence begins.
    # cf. https://docs.python.org/2/library/stdtypes.html#file.next
    for raw_line in iter(f.readline, ''):
        line = raw_line.strip()
        if not line:
            continue
        if line[0] == '>':
            if cur_seq:
                assert cur_name
                yield cur_seq, cur_name, cur_pos
                cur_seq, cur_name = '', ''
            cur_name = line[1:].strip()
            cur_pos = f.tell() - len(raw_line)
            cur_seq = ''
        else:
            assert cur_name
            cur_seq += line
    if cur_seq:
        assert cur_name
        yield cur_seq, cur_name, cur_pos


# NOTE sam alignments don't have substitutions, so outputs here are MID,
# whereas outputs of MAF are in {SMID}*.
def sam_to_opseqs(sam_path, opseq_path):
    """ Usage:
        >>> sam_to_opseqs('foo.sam', 'foo_opseqs.fa')

    The opseqs path will contain a single string of mutations obtained from SAM
    alignments.
    """
    def _opseq_from_cigar_string(cigarstring):
        cigar_ops = [x for x in re.split('(I|M|D)', cigarstring) if x]
        assert len(cigar_ops) % 2 == 0, 'incomplete cigar string'
        ops = []

        # clip flanking insertion/deletions
        if cigar_ops[1] in 'DI':
            cigar_ops = cigar_ops[2:]
        if cigar_ops[-1] in 'DI':
            cigar_ops = cigar_ops[:-2]

        for i in range(len(cigar_ops) / 2):
            ops += [cigar_ops[2 * i + 1]] * int(cigar_ops[2 * i])
        assert all(op in 'IMD' for op in ops)
        return ''.join(ops)

    samfile = pysam.AlignmentFile(sam_path)
    log('loading SAM records from %s.' % sam_path)
    with open(opseq_path, 'w') as f:
        for rec in samfile.fetch():
            opseq = _opseq_from_cigar_string(rec.cigarstring)
            assert opseq.count('I') + opseq.count('M') + opseq.count('S') == \
                rec.query_alignment_length
            f.write(opseq)
        f.write('\n')


def seq_in_maf(alignment, id_):
    for rec in alignment:
        if rec.id == id_:
            return np.array(list(rec))
    return np.array([])


def get_seqs_from_maf(maf_path):
    alignments = list(AlignIO.parse(maf_path, 'maf'))
    # NOTE trust the first alignment to have all the ids
    ids = list(rec.id for rec in alignments[0])
    for i in ids:
        seq = np.concatenate([seq_in_maf(alignment, i)
                              for alignment in alignments])
        assert seq[0] != '-', \
            'Alignments starting with gaps not allowed (pw vs indiv madness)'
        yield i, ''.join(x for x in seq if x != '-')


def get_pws_from_maf(maf_path):
    alignments = list(AlignIO.parse(maf_path, 'maf'))
    # NOTE trust the first alignment to have all the ids
    ids = list(rec.id for rec in alignments[0])

    for i, j in combinations(ids, 2):
        seqs_i, seqs_j = [], []
        for alignment in alignments:
            seq_i_in_aln = seq_in_maf(alignment, i)
            seq_j_in_aln = seq_in_maf(alignment, j)
            if len(seq_i_in_aln) and len(seq_j_in_aln):
                seqs_i.append(seq_i_in_aln)
                seqs_j.append(seq_j_in_aln)
        full_alignment = np.array([np.concatenate(seqs_i),
                                   np.concatenate(seqs_j)])

        # remove double gaps, e.g.
        # AA-T
        # AA-T
        double_gaps = []
        for idx in range(len(full_alignment[0, :])):
            if full_alignment[0, idx] == full_alignment[1, idx] == '-':
                double_gaps.append(idx)
        full_alignment = np.delete(full_alignment, double_gaps, axis=1)

        # translate alignment to SMDI opseq
        def _let_pair_to_op(let0, let1):
            assert let0 != '-' or let1 != '-'
            if let0 == '-':
                return 'I'
            if let1 == '-':
                return 'D'
            if let0 == let1:
                return 'M'
            else:
                return 'S'
        opseq = ''.join(_let_pair_to_op(*full_alignment[:, idx])
                        for idx in range(len(full_alignment[0, :])))
        yield i, j, opseq


# maf_to_individual_fasta('data/actb/actb-7vet.maf', 'data/actb/actb-7vet.fa')
def maf_to_individual_fasta(maf_path, out_path):
    """ Usage:
        >>> maf_to_individual_fasta('foo.maf', 'foo_indiv.fa')

    Each entry in output is in FASTA format where each record is:

        > [id] (as per maf)
        [seq]  (concatenated and - removed as per maf)
    """
    with open(out_path, 'w') as f:
        for id_, seq in get_seqs_from_maf(maf_path):
            f.write('> %s\n%s\n' % (id_, seq))


# maf_to_pw_fasta('data/actb/actb-7vet.maf', 'data/actb/actb-7vet-pws.fa')
def maf_to_pw_fasta(maf_path, out_path):
    """ Usage:
        >>> maf_to_pw_fasta('foo.maf', 'foo_pw.fa')

    Each entry in output is in FASTA format where each record is:

        > [id0]:[id1] (of pair of sequnces as per maf)
        [opsseq]  (in SMID, concatenated and double - removed)
    """
    with open(out_path, 'w') as f:
        for i, j, opseq in get_pws_from_maf(maf_path):
            log('producing pairwise alignment for %s, %s' % (i, j))
            f.write('> %s:%s\n%s\n' % (i, j, opseq))


# =============================================================================
# Random Generation Helpers
# =============================================================================
def seq_pair(n, alphabet, mutation_process=None):
    S = rand_seq(alphabet, n)
    T, _ = mutation_process.mutate(S)
    return S, T


def sample_bio_seqs(ns, source_seq):
    def _intersect(l0, l1):
        r0, r1 = l0 + n, l1 + n
        return min(r0, r1) >= max(l0, l1)
    assert len(source_seq) > sum(ns)
    starts = []
    for n in ns:
        start = None
        while start is None or any(_intersect(start, l) for l in starts):
            start = np.random.randint(0, len(source_seq) - max(ns))
        starts.append(start)
    return [source_seq[s:s + n] for s in starts]


def sample_bio_opseqs(opseq, K, n_samples, gap=None, match=None,
                      resolution=1e-3):
    radius = K / 2
    cands = []
    if match is None:
        matches = []
    else:
        matches = estimate_match_probs_in_opseq(opseq, radius)
    if gap is None:
        gaps = []
    else:
        gaps = estimate_gap_probs_in_opseq(opseq, radius)

    for idx in range(len(opseq) - radius):
        if match is not None and abs(matches[idx] - match) >= resolution:
            continue
        if gap is not None and abs(gaps[idx] - gap) >= resolution:
            continue
        cands.append(idx)
    assert len(cands) > n_samples, 'not enough samples found'
    samples = np.random.choice(cands, size=n_samples)
    for pos in samples:
        yield opseq[pos - radius:pos + radius]


def apply_opseq(S, opseq):
    i_S = 0
    out = ''
    letters = S.alphabet
    for op in opseq:
        if op == 'M':
            out += letters[S[i_S]]
            i_S += 1
        if op == 'S':
            choices = [let for let in letters if let != letters[S[i_S]]]
            out += np.random.choice(choices)
            i_S += 1
        if op == 'D':
            i_S += 1
        if op == 'I':
            out += np.random.choice(letters)
    return S.alphabet.parse(out)


def estimate_gap_probs_in_opseq(opseq, radius):
    assert isinstance(radius, int)
    assert len(opseq) > 2 * radius
    gaps = np.zeros(len(opseq))
    start = opseq[:2 * radius + 1]
    n_MS = start.count('M') + start.count('S')
    for i in range(radius, len(opseq) - radius):
        if opseq[i - radius] in 'MS':
            n_MS -= 1
        if opseq[i + radius] in 'MS':
            n_MS += 1
        gaps[i] = 1 - n_MS / (2. * radius)
    return gaps


def estimate_match_probs_in_opseq(opseq, radius):
    assert isinstance(radius, int)
    assert len(opseq) > 2 * radius
    matches = np.zeros(len(opseq))
    start = opseq[:2 * radius + 1]
    n_M = start.count('M')
    for i in range(radius, len(opseq) - radius):
        if opseq[i - radius] == 'M':
            n_M -= 1
        if opseq[i + radius] == 'M':
            n_M += 1
        matches[i] = n_M / (2. * radius)
    return matches


# =============================================================================
# MATPLOTLIB HELPERS
# =============================================================================
def color_code(values, cmap='nipy_spectral'):
    """Returns an array of RGB colors of the same length as values."""
    cmap = plt.cm.get_cmap(cmap)
    m, M = 1. * min(values), 1. * max(values)
    assert len(values) == len(set(values)), 'all values should be distinct'
    if len(values) > 1:
        return [cmap(.1 + .8 * (v - m) / (M - m)) for v in values]
    else:
        return [cmap(.5)]


def plot_with_sd(ax, xs, ys, axis=None, n_sds=1, y_max=None, color='k', **kw):
    """Plots the mean of given data set (with specified axis) and add a shaded
    region around specifying a given number of standard deviations."""
    assert axis is not None
    means = ys.mean(axis=axis)
    sds = np.sqrt(ys.var(axis=axis))

    ys_l = means - n_sds * sds
    ys_h = means + n_sds * sds

    if y_max is not None:
        assert np.all(means <= y_max), 'all means must be bounded above by max'
        ys_h = np.minimum(ys_h, y_max)

    ax.plot(xs, means, color=color, **kw)
    ax.fill_between(xs, ys_l, ys_h, facecolor=color, edgecolor=color, alpha=.2)


def savefig(fig, path, grid_all=True, comment=None):
    if grid_all:
        for ax in fig.get_axes():
            ax.grid(True)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(os.path.join(PLOT_DIR, path), dpi=270)
    log('created plot%s at %s\n' %
        ('' if comment is None else ' "%s"' % comment, path))


# Call this after all work (e.g. plot_seeds, plot_global_alignment, etc) is
# done on the Axes object.
def adjust_pw_plot(ax, N0, N1):
    ax.set_ylim(-N0 * .1, N0 * 1.1)
    ax.set_xlim(-N1 * .1, N1 * 1.1)
    ax.set_xlabel('Pos. in T', fontsize=12)
    ax.set_ylabel('Pos. in S', fontsize=12)
    ax.set_aspect('equal')
    line_kw = {'c': 'k', 'alpha': .1, 'lw': '5'}
    ax.axvline(x=-5,  **line_kw)
    ax.axvline(x=N1+5, **line_kw)
    ax.axhline(y=-5,  **line_kw)
    ax.axhline(y=N0+5, **line_kw)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.tick_params(labelsize=10)


def plot_seeds(ax, seeds, c='k'):
    idx_S, idx_T = [], []
    for seed in seeds:
        idx_S.append(seed[0])
        idx_T.append(seed[1])

    # x and y are flipped when going from matrix notation to plotting.
    ax.scatter(idx_T, idx_S, facecolor=c, lw=0, s=2, alpha=.3)
    ax.grid(True)


def plot_scored_seeds(ax, ax_colorbar, scored_seeds):
    idx_S, idx_T, cs = [], [], []
    cmap = plt.cm.get_cmap('jet')
    max_score = max(score for _, score in scored_seeds)
    for (i, j), score in scored_seeds:
        idx_S.append(i)
        idx_T.append(j)
        cs.append(cmap(score/max_score)[:3])
    # x and y are flipped when going from matrix notation to plotting.
    ax.scatter(idx_T, idx_S, facecolor=cs, lw=0, s=2, alpha=.9)
    norm = matplotlib.colors.Normalize(vmin=0, vmax=max_score)
    matplotlib.colorbar.ColorbarBase(ax_colorbar, cmap=cmap, norm=norm,
                                     orientation='vertical')


def plot_similar_segment(ax, segment, color='k', **kw):
    (d_min, d_max), (a_min, a_max) = segment
    seg_ds = [d_min, d_min, d_max, d_max]
    seg_as = [a_min, a_max, a_max, a_min]

    seg_xs = [0, 0, 0, 0]
    seg_ys = [0, 0, 0, 0]
    for i in range(4):
        d, a = seg_ds[i], seg_as[i]
        x, y = (a + max(d, 0), a - min(d, 0))
        # x and y is flipped between matrix notation and plotting
        seg_xs[i], seg_ys[i] = y, x

    ax.fill(seg_xs, seg_ys, facecolor=color, edgecolor=color, **kw)


def plot_global_alignment(ax, opseq, **kw):
    xs, ys = [0], [0]
    for op in opseq:
        if op in 'MS':
            xs.append(xs[-1] + 1)
            ys.append(ys[-1] + 1)
        elif op == 'I':
            xs.append(xs[-1])
            ys.append(ys[-1] + 1)
        elif op == 'D':
            xs.append(xs[-1] + 1)
            ys.append(ys[-1])

    ax.plot(xs, ys, **kw)


# =============================================================================
# STATISTICAL HELPERS
#
# Cumulative distribution handling and plots for binary classifiers
# =============================================================================
def empirical_cdf(sample):
    """Empirical cumulative distribution function (returns a callable)."""
    m, M = min(sample), max(sample)
    sorted_ = np.sort(sample)
    yvals = np.arange(len(sorted_) + 1)/float(len(sorted_))

    def _empirical_cdf(xs):
        cdfs = []
        for x in xs:
            if x >= M:
                cdfs.append(1)
            elif x <= m:
                cdfs.append(0)
            else:
                cdfs.append(yvals[bisect_left(sorted_, x)])
        return np.array(cdfs)
    return _empirical_cdf


def plot_cdf(ax, sample, num_points=1000, **kwargs):
    sample = np.sort(sample)
    F = empirical_cdf(sample)
    m, M = sample[0], sample[-1]
    w = .05 * (M - m)
    xs = np.linspace(m - w, M + w, num_points)
    ax.plot(xs, F(xs), **kwargs)
    ax.set_ylim([-.1, 1.1])
    ax.set_title('Cumulative Distribution')


def roc(pos, neg, classifier='>', num_points=1000):
    """Samples the ROC curve for given positive and negative
    samples."""
    assert classifier in '><'
    F_x, F_y = empirical_cdf(pos), empirical_cdf(neg)
    all_ = np.sort(np.concatenate([pos, neg]))
    m, M = all_[0], all_[-1]
    w = .25 * (M - m)
    xs = np.linspace(m - w, M + w, num_points)
    if classifier == '>':
        return (xs, 1 - F_y(xs), 1 - F_x(xs))
    else:
        return (xs, F_y(xs), F_x(xs))


def plot_roc(ax, pos, neg, classifier='>', **kwargs):
    ax.set_xlim([-.1, 1.1])
    ax.set_ylim([-.1, 1.1])
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    _, x, y = roc(pos, neg, classifier=classifier)

    ax.plot(x, y, **kwargs)
    ax.plot([0, 1], [0, 1], lw=2, ls='--', alpha=.2, c='k')


def ppv(pos, neg, classifier='>', num_points=1000):
    """Positive predictive values for given positive and negative samples. The
    classifier argument indicates the direction of positive classification;
    i.e. '>' ('<') means positive samples are expected to have higher (lower)
    scores than negative samples."""
    assert classifier in '><'
    _all = list(pos) + list(neg)
    params = np.linspace(min(_all), max(_all), num_points)
    F_pos, F_neg = empirical_cdf(pos)(params), empirical_cdf(neg)(params)

    if classifier == '>':
        Nc_pos, Nc_neg = (1 - F_pos) * len(pos), (1 - F_neg) * len(neg)
        with np.errstate(divide='ignore', invalid='ignore'):
            # don't complain about points where denominator is 0
            return params, Nc_pos / (Nc_pos + Nc_neg)
    else:
        N_pos, N_neg = F_pos * len(pos), F_neg * len(neg)
        with np.errstate(divide='ignore', invalid='ignore'):
            return params, N_pos / (N_pos + N_neg)


def npv(pos, neg, classifier='>', num_points=1000):
    """Same as ppv but for negative predictive values."""
    assert classifier in '><'
    _all = list(pos) + list(neg)
    params = np.linspace(min(_all), max(_all), num_points)
    F_pos, F_neg = empirical_cdf(pos)(params), empirical_cdf(neg)(params)

    if classifier == '>':
        N_pos, N_neg = F_pos * len(pos), F_neg * len(neg)
        with np.errstate(divide='ignore', invalid='ignore'):
            return params, N_neg / (N_pos + N_neg)
    else:
        Nc_pos, Nc_neg = (1 - F_pos) * len(pos), (1 - F_neg) * len(neg)
        with np.errstate(divide='ignore', invalid='ignore'):
            return params, Nc_neg / (Nc_pos + Nc_pos)


def plot_roc_pv(ax, pos, neg, classifier='>', num_points=1000, **kwargs):
    params, npvs = npv(pos, neg, classifier=classifier, num_points=num_points)
    params, ppvs = ppv(pos, neg, classifier=classifier, num_points=num_points)
    ax.plot(ppvs, npvs, **kwargs)
    ax.set_xlabel('Positive Predictive Value')
    ax.set_ylabel('Negative Predictive Value')
    ax.set_title('Predictive ROC')


def plot_classifier(path, pos, neg, labels=None, classifier='>'):
    if labels is None:
        labels = ['positive', 'negative']
    assert len(labels) == 2

    fig = plt.figure(figsize=(12, 8))

    grids = gridspec.GridSpec(2, 2, height_ratios=[2.5, 1])
    ax_roc = fig.add_subplot(grids[0, 0])
    ax_p_roc = fig.add_subplot(grids[0, 1])
    ax_cdf = fig.add_subplot(grids[1, :])

    plot_cdf(ax_cdf, pos, color='g', label=labels[0], lw=1.5)
    plot_cdf(ax_cdf, neg, color='r', label=labels[1], lw=1.5)
    ax_cdf.legend(fontsize=12, loc='lower right')

    kw = {'lw': 1.5, 'c': 'k'}
    plot_roc(ax_roc, pos, neg, classifier=classifier, **kw)
    roc_ticks = [i * .1 for i in range(11)]
    ax_roc.set_xticks(roc_ticks)
    ax_roc.set_yticks(roc_ticks)

    plot_roc_pv(ax_p_roc, pos, neg, classifier=classifier, **kw)

    for [ax, title] in zip([ax_cdf, ax_roc, ax_p_roc],
                           ['Score CDF', 'ROC', 'Predictive ROC']):
        ax.set_title(title, fontsize=12)
        ax.tick_params(labelsize=12)

    savefig(fig, path, comment='classifier plot')
