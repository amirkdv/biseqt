# -*- coding: utf-8 -*-
import numpy as np
import scipy.stats
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvas


matplotlib.rc('font', style='normal', size=8)

_default_scatter_kwargs = {
    's': 1,
    'color': 'k',
}
_default_savefig_kwargs = {
    'dpi': 500,
}


def _scatter(ax, *args, **kwargs):
    _kwargs = _default_scatter_kwargs.copy()
    _kwargs.update(kwargs)
    ax.scatter(*args, **_kwargs)


# Initializes new figure with the AGG backend_ and adds a canvas to it.
# cf. http://matplotlib.org/faq/usage_faq.html#what-is-a-backend
def figure():
    fig = Figure()
    FigureCanvas(fig)  # add a canvas to figure
    return fig


# kwargs are passed as is to matplotlib.figure.Figure.savefig
def save(fig, path, **kwargs):
    for ax in fig.get_axes():
        ax.grid(True)
        ax.legend()
    fig.tight_layout()
    _kwargs = _default_savefig_kwargs
    _kwargs.update(kwargs)
    fig.savefig(path, **_kwargs)


# Gaussian density estimate of given sample. kwargs are passed as is to
# matplotlib.axes.Axes.plot.
def plot_density(ax, sample, num_points=1000, **kwargs):
    density = scipy.stats.gaussian_kde(sample)
    xs = np.linspace(min(sample), max(sample), num=num_points)
    ax.plot(xs, density(xs), **kwargs)


# Cumulative distribution obtained by integrating estimated density. Value
# range enforces range of values over which CDF is calculated.
# Returns xs, cdf
def cdf(sample, num_points=1000, value_range=None):
    density = scipy.stats.gaussian_kde(sample)
    if value_range is None:
        value_range = min(sample), max(sample)
    xs = np.linspace(*value_range, num=num_points)
    cdf = scipy.integrate.cumtrapz(density(xs), xs)
    return xs, np.insert(cdf, [0], 0)


def plot_cdf(ax, sample, num_points=1000, **kwargs):
    xs, F = cdf(sample, num_points=num_points)
    ax.set_title('Cumulative Distribution')
    _scatter(ax, xs, F, **kwargs)


# Samples the ROC curve for given positive (X) and negative (Y) samples.
def roc(X, Y, classifier='>', num_points=1000):
    assert classifier in '><'
    _all = list(X) + list(Y)
    value_range = min(_all), max(_all)
    xs, F_X = cdf(X, num_points=num_points, value_range=value_range)
    _, F_Y = cdf(Y, num_points=num_points, value_range=value_range)
    return (xs, 1 - F_Y, 1 - F_X) if classifier == '>' else (xs, F_Y, F_X)


def plot_roc(ax, X, Y, **kwargs):
    ax.set_xlim([-.1, 1.1])
    ax.set_ylim([-.1, 1.1])
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    _, x, y = roc(X, Y)
    _scatter(ax, x, y, **kwargs)


# Positive predictive value for given positive (X) and negative (Y) samples.
def ppv(X, Y, classifier='>', num_points=1000):
    assert classifier in '><'
    _all = list(X) + list(Y)
    value_range = min(_all), max(_all)
    params, F_X = cdf(X, num_points=num_points, value_range=value_range)
    _, F_Y = cdf(Y, num_points=num_points, value_range=value_range)

    if classifier == '>':
        Nc_X, Nc_Y = (1 - F_X) * len(X), (1 - F_Y) * len(Y)
        with np.errstate(divide='ignore', invalid='ignore'):
            # don't complain about points where denominator is 0
            return params, Nc_X / (Nc_X + Nc_Y)
    else:
        N_X, N_Y = F_X * len(X), F_Y * len(Y)
        with np.errstate(divide='ignore', invalid='ignore'):
            return params, N_X / (N_X + N_Y)


# Negative predictive value for given positive (X) and negative (Y) samples.
def npv(X, Y, classifier='>', num_points=1000):
    assert classifier in '><'
    _all = list(X) + list(Y)
    value_range = min(_all), max(_all)
    params, F_X = cdf(X, num_points=num_points, value_range=value_range)
    _, F_Y = cdf(Y, num_points=num_points, value_range=value_range)

    if classifier == '>':
        N_X, N_Y = F_X * len(X), F_Y * len(Y)
        with np.errstate(divide='ignore', invalid='ignore'):
            return params, N_Y / (N_X + N_Y)
    else:
        Nc_X, Nc_Y = (1 - F_X) * len(X), (1 - F_Y) * len(Y)
        with np.errstate(divide='ignore', invalid='ignore'):
            return params, Nc_Y / (Nc_X + Nc_X)


def plot_ppv(ax, X, Y, classifier='>', num_points=1000, **kwargs):
    params, ppvs = ppv(X, Y, classifier=classifier, num_points=num_points)
    ax.set_title('Positive Predictive Value')
    _scatter(ax, params, ppvs, **kwargs)


def plot_npv(ax, X, Y, classifier='>', num_points=1000, **kwargs):
    params, npvs = npv(X, Y, classifier=classifier, num_points=num_points)
    ax.set_title('Negative Predictive Value')
    _scatter(ax, params, npvs, **kwargs)
