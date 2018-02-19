import os
import sys
import pickle
import numpy as np
from bisect import bisect_left
from textwrap import TextWrapper
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec


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
# IO Helpers
# =============================================================================
# TODO if we don't need file position use biopython instead
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


# TODO change syntax to pickle_dump(obj, path)
def pickle_dump(path, obj, comment=None):
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
        pickle_dump(dumpfile, out)
        return out

    return wrapped_func


# =============================================================================
# MATPLOTLIB HELPERS
# =============================================================================
def color_code(values, cmap='rainbow'):
    """Returns an array of RGB colors of the same length as values."""
    cmap = plt.cm.get_cmap(cmap)
    m, M = 1. * min(values), 1. * max(values)
    return [cmap((v - m) / (M - m)) for v in values]


def plot_with_sd(ax, xs, ys, axis=None, n_sds=1, y_max=None, color='k', **kw):
    """Plots the mean of given data set (with specified axis) and add a shaded
    region around specifying a given number of standard deviations."""
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
    log('recreated plot%s at %s\n' %
        ('' if comment is None else ' "%s"' % comment, path))


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


def plot_classifier(pos, neg, labels=None, name=None, classifier='>'):
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

    savefig(fig, 'classifier_%s.png' % name, comment='classifier plot')
