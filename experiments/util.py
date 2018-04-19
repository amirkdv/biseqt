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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from Bio import AlignIO
from itertools import combinations
import pysam
from biseqt.stochastics import rand_seq
from scipy.ndimage.filters import gaussian_filter1d


plt.rc('text', usetex=True)


def log(msg, newline=True):
    wrapper = TextWrapper(initial_indent='* ', subsequent_indent='  ',
                          width=120)
    sys.stderr.write('\n'.join(wrapper.wrap(msg)) + ('\n' if newline else ''))
    sys.stderr.flush()


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
def load_fasta(f, num_seqs=-1):
    """Given file handle reads fasta sequences and yields tuples of (seq, name,
       pos) containing the raw sequence, its FASTA name, and its starting
       position in f.
    """
    cur_name = cur_seq = ''
    cur_pos = 0
    cnt = 0

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
                cnt += 1
                if num_seqs > 0 and cnt == num_seqs:
                    break
                cur_seq, cur_name = '', ''
            cur_name = line[1:].strip()
            cur_pos = f.tell() - len(raw_line)
            cur_seq = ''
        else:
            assert cur_name
            cur_seq += line
    if cur_seq:
        if num_seqs < 0 or cnt < num_seqs:
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


def seq_in_mse(alignment, id_):
    for rec in alignment:
        if rec.id == id_:
            return np.array(list(rec))
    return np.array([])


def get_seqs_from_mse(mse_path, fmt='maf'):
    alignments = list(AlignIO.parse(mse_path, fmt))
    ids = set()
    for alignment in alignments:
        ids = ids.union(set(rec.id for rec in alignment))
    for i in ids:
        seq = np.concatenate([seq_in_mse(alignment, i)
                              for alignment in alignments])
        yield i, ''.join(x for x in seq if x != '-')


def pws_from_mse(mse_path, fmt='maf'):
    alignments = list(AlignIO.parse(mse_path, fmt))
    # NOTE trust the first alignment to have all the ids
    ids = list(rec.id for rec in alignments[0])

    for i, j in combinations(ids, 2):
        seqs_i, seqs_j = [], []
        for alignment in alignments:
            seq_i_in_aln = seq_in_mse(alignment, i)
            seq_j_in_aln = seq_in_mse(alignment, j)
            if len(seq_i_in_aln) == 0:
                seq_i_in_aln = ['-'] * len(seq_j_in_aln)
            elif len(seq_j_in_aln) == 0:
                seq_j_in_aln = ['-'] * len(seq_i_in_aln)
            assert len(seq_i_in_aln) == len(seq_j_in_aln)
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
                if let0 == 'N':
                    # unknown nucleotides shouldn't match with each other.
                    return 'S'
                else:
                    return 'M'
            else:
                return 'S'
        opseq = ''.join(_let_pair_to_op(*full_alignment[:, idx])
                        for idx in range(len(full_alignment[0, :])))
        yield i, j, opseq


# some opseqs have long chains of insertion or deletion. This leads to weird
# similarity stats (e.g. segment of length 150 with 98% match followed by 100
# bases of deletion gives a total match probability of ~ 50%). For the purpose
# of identifying the correct length of alignment we remove any stretch of N
# deletions or insertions from the opseq
def pws_to_opseq(pw_path, opseq_path, max_gap_length=100):
    def _remove_long_gaps(opseqs, max_gap_length):
        p = re.compile('D{%d,}|I{%d,}' % (max_gap_length, max_gap_length))
        for opseq in opseqs:
            yield p.sub('', opseq)

    with open(pw_path) as f_pw:
        log('loading pairwise opseqs from %s' % pw_path)
        opseqs = [(name, opseq) for opseq, name, _ in load_fasta(f_pw)]

    opseqs_clean = list(_remove_long_gaps([opseq for name, opseq in opseqs],
                                          max_gap_length))
    with open(opseq_path, 'w') as f_ops:
        log('writing cleaned up opseqs to %s' % pw_path)
        for (name, opseq), opseq_clean in zip(opseqs, opseqs_clean):
            log('[%s] removed %d edit ops (out of %d) from opseq for %s' %
                (opseq_path, len(opseq) - len(opseq_clean), len(opseq), name))
            f_ops.write('> (cleaned, max_gap = %d) %s\n%s\n' %
                        (max_gap_length, name, opseq_clean))


def mse_to_individual_fasta(mse_path, out_path, fmt='maf'):
    """ Usage:
        >>> mse_to_individual_fasta('foo.maf', 'foo_indiv.fa')

    Each entry in output is in FASTA format where each record is:

        > [id] (as per maf)
        [seq]  (concatenated and - removed as per maf)
    """
    with open(out_path, 'w') as f:
        for id_, seq in get_seqs_from_mse(mse_path, fmt=fmt):
            log('[%s] extracting "%s" (%d) from MSE' %
                (out_path, id_, len(seq)))
            f.write('> %s\n%s\n' % (id_, seq))


def mse_to_pw_fasta(mse_path, out_path, fmt='maf'):
    """ Usage:
        >>> mse_to_pw_fasta('foo.maf', 'foo_pw.fa')

    Each entry in output is in FASTA format where each record is:

        > [id0]:[id1] (of pair of sequnces as per MSE)
        [opsseq]  (in SMID, concatenated and double - removed)
    """
    with open(out_path, 'w') as f:
        for i, j, opseq in pws_from_mse(mse_path, fmt=fmt):
            log('[%s] producing pairwise alignment for %s, %s' %
                (out_path, i, j))
            f.write('> %s:%s\n%s\n' % (i, j, opseq))


# =============================================================================
# Random Generation Helpers
# =============================================================================
def seq_pair(n, alphabet, mutation_process=None):
    S = rand_seq(alphabet, n)
    T, _ = mutation_process.mutate(S)
    return S, T


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


# if projection is either of 1 or 2, the positions are projected on either of
# the first or second sequence (origin or mutatnt); if None positions are
# against the opseq.
def estimate_match_probs_in_opseq(opseq, radius, projection=None):
    assert isinstance(radius, int)
    assert all(op in 'MSID' for op in opseq)
    assert len(opseq) > 2 * radius
    probs = np.zeros(len(opseq))
    start = opseq[:2 * radius]
    n_M = start.count('M')
    for i in range(radius, len(opseq) - radius):
        probs[i] = n_M / (2. * radius)
        if opseq[i - radius] == 'M':
            n_M -= 1
        if opseq[i + radius] == 'M':
            n_M += 1

    assert len(probs) == len(opseq)

    if projection is not None:
        probs_proj = np.zeros(len(opseq))  # guaranteed upper bound for length
        assert projection in [1, 2]
        i, j = 0, 0
        for op, prob in zip(opseq, probs):
            probs_proj[{1: i, 2: j}[projection]] = prob
            if op in 'MSD':
                i += 1
            if op in 'MSI':
                j += 1
        probs = probs_proj[:{1: i, 2: j}[projection] + 1]

    return probs


def fill_in_unknown(seq, alphabet):
    def _filter(s): return s if s != 'N' else np.random.choice(alphabet)

    return ''.join(_filter(s) for s in seq)


def seeds_from_opseq(opseq, wordlen):
    i, j = 0, 0
    for idx, op in enumerate(opseq):
        if opseq[idx:idx + wordlen] == 'M' * wordlen:
            yield idx, (i, j)
        if op in 'MS':
            i, j = i + 1, j + 1
        elif op in 'D':
            i = i + 1
        elif op in 'I':
            j = j + 1
        else:
            raise ValueError('op %s not understood' % op)


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
    ax.fill_between(xs, ys_l, ys_h, facecolor=color, edgecolor=color, alpha=.1)


def savefig(fig, path, dpi=270, grid_all=True, comment=None):
    if grid_all:
        for ax in fig.get_axes():
            ax.grid(True)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(os.path.join(PLOT_DIR, path), dpi=dpi)
    log('created plot%s at %s\n' %
        ('' if comment is None else ' "%s"' % comment, path))
    plt.close(fig)


# Call this after all work (e.g. plot_seeds, plot_global_alignment, etc) is
# done on the Axes object.
def adjust_pw_plot(ax, N0, N1):
    ax.set_ylim(-N0 * .1, N0 * 1.1)
    ax.set_xlim(-N1 * .1, N1 * 1.1)
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
    xs, ys = [], []
    for seed in seeds:
        ys.append(seed[0])
        xs.append(seed[1])

    ax.scatter(xs, ys, facecolor=c, lw=0, s=2, alpha=.3)
    ax.grid(True)


def plot_scored_seeds(ax, scored_seeds, threshold=5, extent=None, **kw):
    idx_S, idx_T, cs, ss = [], [], [], []
    cmap = plt.cm.get_cmap('Greys')
    if extent is None:
        max_score = max(score for _, score in scored_seeds)
        min_score = min(score for _, score in scored_seeds)
    else:
        min_score, max_score = extent
    for (i, j), score in scored_seeds:
        idx_S.append(i)
        idx_T.append(j)
        cs.append(cmap(score/max_score)[:3])
        ss.append(10 if threshold and score > threshold else 1)
    # x and y are flipped when going from matrix notation to plotting.
    ax.scatter(idx_T, idx_S, facecolor=cs, s=ss, lw=0, **kw)

    # colorbar Axes
    ax_c = make_axes_locatable(ax).append_axes('right', size='4%', pad=0.05)

    norm = matplotlib.colors.Normalize(vmin=min_score, vmax=max_score)
    matplotlib.colorbar.ColorbarBase(ax_c, cmap=cmap, norm=norm,
                                     orientation='vertical')


# NOTE hardcodes coordinate conversion logic
def plot_similar_segment(ax, segment, color='k', **kw):
    (d_min, d_max), (a_min, a_max) = segment
    seg_ds = [d_min, d_min, d_max, d_max]
    seg_as = [a_min, a_max, a_max, a_min]

    seg_xs = [0, 0, 0, 0]
    seg_ys = [0, 0, 0, 0]
    for i in range(4):
        d, a = seg_ds[i], seg_as[i]
        x, y = (a + max(d, 0), a - min(d, 0))
        x, y = (a + d) / 2, (a - d) / 2
        # x and y is flipped between matrix notation and plotting
        seg_xs[i], seg_ys[i] = y, x

    ax.fill(seg_xs, seg_ys, facecolor=color, edgecolor=color, **kw)


def opseq_path(opseq, x0=0, y0=0):
    xs, ys = [x0], [y0]
    for op in opseq:
        if op in 'MS':
            xs.append(xs[-1] + 1)
            ys.append(ys[-1] + 1)
        elif op == 'D':
            xs.append(xs[-1] + 1)
            ys.append(ys[-1])
        elif op == 'I':
            xs.append(xs[-1])
            ys.append(ys[-1] + 1)
    return xs, ys


def plot_global_alignment(ax, opseq, **kw):
    ys, xs = opseq_path(opseq)
    ax.plot(xs, ys, **kw)


def plot_local_alignment(ax, opseq, i, j, **kw):
    ys, xs = opseq_path(opseq, x0=i, y0=j)
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


def plot_cdf(ax, sample, **kw):
    smooth_radius = kw.pop('smooth_radius', None)
    num_points = kw.pop('num_points', 1000)
    mark_threshold = kw.pop('mark_threshold', None)

    sample = np.sort(sample)
    F = empirical_cdf(sample)
    m, M = sample[0], sample[-1]
    w = .05 * (M - m)
    xs = np.linspace(m - w, M + w, num_points)
    if smooth_radius is None:
        ys = F(xs)
    else:
        ys = gaussian_filter1d(F(xs), smooth_radius)
    ax.plot(xs, ys, **kw)
    ax.set_ylim([-.1, 1.1])
    ax.set_title('Cumulative Distribution')
    if mark_threshold is not None:
        ax.axvline(mark_threshold, c='b', ls='--', lw=3, alpha=.2)


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


def plot_roc(ax, pos, neg, classifier='>', mark_threshold=None, **kw):
    ax.set_xlim([-.1, 1.1])
    ax.set_ylim([-.1, 1.1])
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    params, x, y = roc(pos, neg, classifier=classifier)

    ax.plot(x, y, **kw)
    ax.plot([0, 1], [0, 1], lw=2, ls='--', alpha=.2, c='k')

    if mark_threshold is not None:
        idx = bisect_left(params, mark_threshold)
        ax.scatter([x[idx]], [y[idx]], c='b', s=100, lw=0, alpha=.3)


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


def plot_roc_pv(ax, pos, neg, classifier='>', **kw):
    mark_threshold = kw.pop('mark_threshold', None)
    num_points = kw.pop('num_points', 1000)

    params, npvs = npv(pos, neg, classifier=classifier, num_points=num_points)
    params, ppvs = ppv(pos, neg, classifier=classifier, num_points=num_points)
    ax.plot(ppvs, npvs, **kw)
    ax.set_xlabel('Positive Predictive Value')
    ax.set_ylabel('Negative Predictive Value')
    ax.set_title('Predictive ROC')
    if mark_threshold is not None:
        idx = bisect_left(params, mark_threshold)
        ax.scatter([ppvs[idx]], [npvs[idx]], c='b', s=100, lw=0, alpha=.3)


def plot_classifier(pos, neg, labels=None, classifier='>', title='',
                    mark_threshold=None):
    if labels is None:
        labels = ['positive', 'negative']
    assert len(labels) == 2

    fig = plt.figure(figsize=(12, 8))

    grids = gridspec.GridSpec(2, 2, height_ratios=[2.5, 1])
    ax_roc = fig.add_subplot(grids[0, 0])
    ax_p_roc = fig.add_subplot(grids[0, 1])
    ax_cdf = fig.add_subplot(grids[1, :])

    plot_cdf(ax_cdf, pos, mark_threshold=mark_threshold, color='g',
             label=labels[0], lw=1.5)
    plot_cdf(ax_cdf, neg, mark_threshold=mark_threshold, color='r',
             label=labels[1], lw=1.5)
    ax_cdf.legend(fontsize=12, loc='lower right')

    kw = {'lw': 1.5, 'c': 'k'}
    plot_roc(ax_roc, pos, neg, classifier=classifier,
             mark_threshold=mark_threshold, **kw)
    roc_ticks = [i * .1 for i in range(11)]
    ax_roc.set_xticks(roc_ticks)
    ax_roc.set_yticks(roc_ticks)

    plot_roc_pv(ax_p_roc, pos, neg, classifier=classifier,
                mark_threshold=mark_threshold, **kw)

    for [ax, ax_title] in zip([ax_cdf, ax_roc, ax_p_roc],
                              ['Classifier statistic CDF',
                               'ROC for %d(+), %d(-) samples' %
                               (len(pos), len(neg)),
                               'Predictive ROC']):
        ax.set_title(ax_title, fontsize=12)
        ax.tick_params(labelsize=12)

    return fig, [ax_cdf, ax_roc, ax_p_roc]
