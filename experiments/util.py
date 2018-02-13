import os
import sys
import pickle
import numpy as np
from matplotlib import pyplot as plt
from bisect import bisect_left
from textwrap import TextWrapper


plt.rc('text', usetex=True)


def log(msg, newline=True):
    wrapper = TextWrapper(initial_indent='* ', subsequent_indent='  ',
                          width=100)
    sys.stderr.write('\n'.join(wrapper.wrap(msg)) + '\n' if newline else '')


# =============================================================================
# File IO Helpers
# =============================================================================
# TODO change syntax to pickle_dump(obj, path)
def pickle_dump(path, obj, comment=None):
    try:
        with open(path, 'w') as f:
            pickle.dump(obj, f)
    except IOError, e:
        log('Error writing to %s: %s' % (path, str(e)))
        n = np.random.choice(list('0123456789abcdef'), size=10)
        path = '/tmp/' + ''.join(str(i) for i in n)
        with open(path, 'w') as f:
            pickle.dump(obj, f)
    log('dumped %s to %s.' % (comment if comment else 'object', path))


def pickle_load(path):
    with open(path, 'rb') as f:
        return pickle.load(f)


def with_dumpfile(func):
    """Decorator for storing simulation results, usage:

        @with_dumpfile
        sim_foo(..., **kw):
            ...
            return results # this is dumped to dumpfile

        # load from dumpfile if it exists, sim_foo() will not be executed
        res = sim_foo(..., dumpfile='foo.txt', ignore_existing=False)
        # force a re-run of sim_foo(), dumpfile will be overwritten
        res = sim_foo(..., dumpfile='foo.txt', ignore_existing=True)
    """
    def wrapped_func(*args, **kwargs):
        dumpfile = kwargs['dumpfile']
        if os.path.exists(dumpfile) and not kwargs.get('repeat', False):
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
    return [cmap((v - m) / M) for v in values]


def plot_with_sd(ax, xs, ys, num_sds=1, axis=None, color='k', **plot_kw):
    """Plots the mean of given data set (with specified axis) and add a shaded
    region around specifying a given number of standard deviations."""
    means = ys.mean(axis=axis)
    sds = np.sqrt(ys.var(axis=axis))

    ys_l = means - num_sds * sds
    ys_h = means + num_sds * sds

    ax.plot(xs, means, color=color, **plot_kw)
    ax.fill_between(xs, ys_l, ys_h, facecolor=color, edgecolor=color, alpha=.2)


def savefig(fig, path, grid_all=True, comment=None):
    if grid_all:
        for ax in fig.get_axes():
            ax.grid(True)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(path, dpi=270)
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

    import matplotlib.gridspec as gridspec
    gs = gridspec.GridSpec(2, 2, height_ratios=[2.5, 1])
    ax_roc = fig.add_subplot(gs[0, 0])
    ax_p_roc = fig.add_subplot(gs[0, 1])
    ax_cdf = fig.add_subplot(gs[1, :])

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
