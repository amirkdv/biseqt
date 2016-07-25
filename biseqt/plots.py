# -*- coding: utf-8 -*-
"""
.. wikisection:: overview
    :title: Plots

    The :mod:`biseqt.plots` module provides plotting utilities based on
    ``matplotlib`` functionality in addition to basic tools for inspecting
    binary classifiers.

    >>> from biseqt.plots import figure, save, cdf, plot_density
    >>> X = ... # a list of numbers (the positive data set)
    >>> Y = ... # a list of numbers (the negatvie data set)
    >>> fig = figure()
    >>> plot_density(fig.subplot(2, 1, 1), X)
    >>> plot_density(fig.subplot(2, 1, 2), Y)
    >>> save(fig, '/path/to/fig.png')

    >>> from biseqt.plots import plot_roc, plot_ppv
    >>> fig = figure()
    >>> plot_roc(fig.subplot(1, 1, 1), X, Y)

    >>> fig = figure()
    >>> plot_ppv(fig.subplot(1, 1, 1), X, Y)
"""

import numpy as np
import scipy.stats
import matplotlib


matplotlib.rc('font', style='normal', size=8)

_default_scatter_kwargs = {
    's': 1,
    'color': 'k',
}


def _scatter(ax, *args, **kwargs):
    _kwargs = _default_scatter_kwargs.copy()
    _kwargs.update(kwargs)
    ax.scatter(*args, **_kwargs)


def figure():
    """Initializes new figure with the AGG backend_ and adds a canvas to it.

    .. _backend: http://matplotlib.org/faq/usage_faq.html#what-is-a-backend

    Returns:
        matplotlib.figure.Figure
    """
    fig = matplotlib.figure.Figure()
    matplotlib.backends.backend_agg.FigureCanvas(fig)  # add a canvas to figure
    return fig


def save(fig, path, **kwargs):
    """Saves figure to path. All keyword arguments are passed as is to
    ``matplotlib.figure.Figure.savefig``.

    Args:
        fig (matplotlib.figure.Figure)
        path (str)
    """
    for ax in fig.get_axes():
        ax.grid(True)
        ax.legend()
    fig.tight_layout()
    fig.savefig(path, **kwargs)


def plot_density(ax, sample, num_points=1000, **kwargs):
    """Plots the Guassian density estimate of the given sample. All un-named
    keyword arguments are passed as is to ``matplotlib.axes.Axes.plot``.

    Args:
        ax (matplotlib.axes.Axes): The matplotlib axes object to use.
        sample (list): A list of numbers; the sample set.

    Keyword Args:
        num_points (int): The sample size of the density estimate.
    """
    density = scipy.stats.gaussian_kde(sample)
    xs = np.linspace(min(sample), max(sample), num=num_points)
    ax.plot(xs, density(xs), **kwargs)


def cdf(sample, num_points=1000, value_range=None):
    """Calculate the cumulative distribution function of a given sample set.

    Args:
        sample (list): A list of numbers; the sample set.

    Keyword Args:
        num_points (int): The sample size of the CDF.
        value_range (tuple): The minimum and maximum values where the CDF
            should be calculated; default is None in which case the bounds of
            the sample set are used.

    Returns:
        tuple:
            A tuple containing two lists: the values at which the CDF is
            reported and CDF values themselves.
    """
    density = scipy.stats.gaussian_kde(sample)
    if value_range is None:
        value_range = min(sample), max(sample)
    xs = np.linspace(*value_range, num=num_points)
    cdf = scipy.integrate.cumtrapz(density(xs), xs)
    return xs, np.insert(cdf, [0], 0)


def plot_cdf(ax, sample, num_points=1000, **kwargs):
    """Plots the cumulative distribution of a given sample set. All unnamed
    keyword arguments are passed as is to ``matplotlib.axes.Axes.scatter``.

    Args:
        ax (matplotlib.axes.Axes): The matplotlib axes object to use.
        sample (list): A list of numbers; the sample set.

    Keyword Args:
        num_points (int): The sample size of the CDF.
    """
    xs, F = cdf(sample, num_points=num_points)
    ax.set_title('Cumulative Distribution')
    _scatter(ax, xs, F, **kwargs)


def roc(X, Y, classifier='>', num_points=1000):
    """Sample the ROC curve for a pair of sample sets labeled as *positive* and
    *negative*. The ROC curve is defined by the
    :math:`(\mathrm{FPR},\mathrm{TPR})` pairs of false and true positive rates
    obtained for various values of threshold. The rates are defined to be:

    .. math::
        \mathrm{FPR} = \\frac{f.p.}{f.p. + t.n.}

        \mathrm{TPR} = \\frac{t.p.}{t.p. + f.n.}

    where :math:`t.p., f.p., t.n., f.n.` stand for true and false positives and
    negatives.  Letting :math:`F_X,F_Y` be the cumulative distributions of the
    sample sets we get:

    .. math::
        \mathrm{FPR} = 1 - F_Y, \\text{ and } \mathrm{TPR} = 1 - F_X

    when the direction for positive selection is ``>``, and:

    .. math::
        \mathrm{FPR} = F_Y, \\text{ and } \mathrm{TPR} = F_X

    when the direction for positive selection is ``<``.

    Args:
        X (list): The "positive" sample set.
        Y (list): The "negative" sample set.

    Keyword Args:
        classifier (str): Either of ``>`` or ``<`` denoting the rule for
            positive classification: the former indicates positive
            classification if observed value is *larger* than threshold and
            the latter indicates positive classification if observed value is
            *smaller* than threshold; default is ``>``, i.e positively
            classified values are thought to be mostly larger than negatively
            classified values.
        num_points (int): The sample size of the ROC curve; translates to
            sample points of the CDF of each sample set.

    Returns:
        tuple:
            A tuple containing three lists: the parameter values (i.e threshold
            values) for which the ROC curve is sampled, the *false positive*
            rates and the *true positive* rates at every sample point.
    """
    assert classifier in '><'
    _all = list(X) + list(Y)
    value_range = min(_all), max(_all)
    xs, F_X = cdf(X, num_points=num_points, value_range=value_range)
    _, F_Y = cdf(Y, num_points=num_points, value_range=value_range)
    return (xs, 1 - F_Y, 1 - F_X) if classifier == '>' else (xs, F_Y, F_X)


def plot_roc(ax, X, Y, **kwargs):
    """Plots the ROC curve of two sample sets, cf. :func:`roc`. All keyword
    arguments are passed as is to ``matplotlib.axes.Axes.scatter``.

    Args:
        ax (matplotlib.axes.Axes): The matplotlib axes object to use.
        X (list): The "positive" sample set (horizontal axis).
        Y (list): The "negative" sample set (vertical axis).
    """
    ax.set_xlim([-.1, 1.1])
    ax.set_ylim([-.1, 1.1])
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    _, x, y = roc(X, Y)
    _scatter(ax, x, y, **kwargs)


def ppv(X, Y, classifier='>', num_points=1000):
    """Samples the *positive predictive value* for a pair of sample sets
    labeled as *positive* and *negative*. The positive negative predictive
    value is defined as:

    .. math::
        \mathrm{PPV} = \\frac{t.p.}{t.p. + f.p.}

    where :math:`t.p., f.p., t.n., f.n.` stand for true and false positives and
    negatives. Depending on the classifier direction this translates to:

    .. math::
        \mathrm{PPV} = \\frac{(1-F_X)|X|}{(1-F_X)|X| + (1-F_Y)|Y|}

    when the direction for positive selection is ``>``, and:

    .. math::
        \mathrm{PPV} = \\frac{F_X|Y|}{F_X|X| + F_Y|Y|}

    when the direction for positive selection is ``<``.

    Args:
        X (list): The "positive" sample set.
        Y (list): The "negative" sample set.

    Keyword Args:
        classifier (str): Either of ``>`` or ``<`` denoting the rule for
            positive classification: the former indicates positive
            classification if observed value is *larger* than threshold and
            the latter indicates positive classification if observed value is
            *smaller* than threshold; default is ``>``, i.e positively
            classified values are thought to be mostly larger than negatively
            classified values.
        num_points (int): The number of sample points.

    Returns:
        tuple:
            A tuple containing two lists: the parameter values where the PPV is
            sampled, and PPV at each point.
    """
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


def npv(X, Y, classifier='>', num_points=1000):
    """Samples the *negative predictive value* for a pair of sample sets
    labeled as *positive* and *negative*. The positive negative predictive
    value is defined as:

    .. math::
        \mathrm{NPV} = \\frac{t.n.}{t.n. + f.n.}

    where :math:`t.p., f.p., t.n., f.n.` stand for true and false positives and
    negatives. Depending on the classifier direction this translates to:

    .. math::
        \mathrm{NPV} = \\frac{F_Y|Y|}{F_X|X| + F_Y|Y|}

    when the direction for positive selection is ``>``, and:

    .. math::
        \mathrm{NPV} = \\frac{(1-F_Y)|Y|}{(1-F_X)|X| + (1-F_Y)|Y|}


    when the direction for positive selection is ``<``.

    Args:
        X (list): The "positive" sample set.
        Y (list): The "negative" sample set.

    Keyword Args:
        classifier (str): Either of ``>`` or ``<``, cf. :func:`ppv`.
        num_points (int): The number of sample points.

    Returns:
        tuple:
            A tuple containing two lists: the parameter values where the PPV is
            sampled, and NPV at each point.
    """
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
    """Plots the positive predictive value over a range of parameters. All
    keyword arguments are passed as is to ``matplotlib.axes.Axes.scatter`.

    Args:
        X (list): The "positive" sample set.
        Y (list): The "negative" sample set.

    Keyword Args:
        classifier (str): Either of ``>`` or ``<``, cf. :func:`ppv`.
        num_points (int): The number of sample points.
    """
    params, ppvs = ppv(X, Y, classifier=classifier, num_points=num_points)
    ax.set_title('Positive Predictive Value')
    _scatter(ax, params, ppvs, **kwargs)


def plot_npv(ax, X, Y, classifier='>', num_points=1000, **kwargs):
    """Plots the negative predictive value over a range of parameters. All
    keyword arguments are passed as is to ``matplotlib.axes.Axes.scatter`.

    Args:
        X (list): The "positive" sample set.
        Y (list): The "negative" sample set.

    Keyword Args:
        classifier (str): Either of ``>`` or ``<``, cf. :func:`ppv`.
        num_points (int): The number of sample points.
    """
    params, npvs = npv(X, Y, classifier=classifier, num_points=num_points)
    ax.set_title('Negative Predictive Value')
    _scatter(ax, params, npvs, **kwargs)
