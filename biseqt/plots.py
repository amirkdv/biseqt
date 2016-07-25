# -*- coding: utf-8 -*-
# TODO figure out proper usage of style sheets so all the font, tight layout,
# dpi, ... gets out of code, cf. http://matplotlib.org/users/style_sheets.html
# and http://matplotlib.org/users/customizing.html

# TODO document http://matplotlib.org/faq/usage_faq.html#what-is-a-backend

import numpy as np
import scipy.stats
import matplotlib


matplotlib.rc('font', style='normal', size=8)


def fig():
    """Initializes new figure.

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
        ax (matplotlib.axes.Axes): The matplotlib axes object to plot things
            on.
        sample (list): A list of numbers; the sample set.

    Keyword Args:
        num_points (int): The sample size of the density estimate.
    """
    density = scipy.stats.gaussian_kde(sample)
    xs = np.linspace(min(sample), max(sample), num=num_points)
    ax.plot(xs, density(xs), **kwargs)


def cdf(sample, num_points=1000):
    """Calculate the cumulative distribution function of a given sample set.

    Args:
        sample (list): A list of numbers; the sample set.

    Keyword Args:
        num_points (int): The sample size of the CDF.

    Returns:
        tuple:
            A tuple containing two lists: the values at which the CDF is
            reported and CDF values themselves.
    """
    density = scipy.stats.gaussian_kde(sample)
    xs = np.linspace(min(sample), max(sample), num=num_points)
    cdf = scipy.integrate.cumtrapz(density(xs), xs)
    return xs, np.insert(cdf, [0], 0)


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
    xs, F_x = cdf(X, num_points=num_points)
    _, F_y = cdf(Y, num_points=num_points)
    return (xs, 1 - F_y, 1 - F_x) if classifier == '>' else (xs, F_y, F_x)


def predictive_value(X, Y, classifier='>', num_points=1000):
    """Samples the PPV-NPV curve for a pair of sample sets, say positive and
    negative label sets. The curve is similar to the ROC curve but captures
    positive and negative *predictive values* instead of sensitivity and
    specificity. The positive and negative predictive values are defined as:

    .. math::
        \mathrm{PPV} = \\frac{t.p.}{t.p. + f.p.}

        \mathrm{NPV} = \\frac{t.n.}{t.n. + f.n.}

    where :math:`t.p., f.p., t.n., f.n.` stand for true and false positives and
    negatives. Depending on the classifier direction this translates to:

    .. math::
        \mathrm{PPV} = \\frac{(1-F_X)|X|}{(1-F_X)|X| + (1-F_Y)|Y|}

        \mathrm{NPV} = \\frac{F_Y|Y|}{F_X|X| + F_Y|Y|}

    when the direction for positive selection is ``>``, and:

    .. math::
        \mathrm{PPV} = \\frac{F_X|Y|}{F_X|X| + F_Y|Y|}

        \mathrm{NPV} = \\frac{(1-F_Y)|Y|}{(1-F_X)|X| + (1-F_Y)|Y|}


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
            A tuple containing three lists: the parameter values at which the
            curve is sampled, and the x and y coordinates of the sample points
            on the ROC curve.
    """
    assert classifier in '><'
    xs, F_x = cdf(X, num_points=num_points)
    _, F_y = cdf(Y, num_points=num_points)

    _X, _Xc = F_x * len(X), (1 - F_x) * len(X)
    _Y, _Yc = F_y * len(Y), (1 - F_y) * len(Y)
    if classifier == '>':
        ppv = _Xc / (_Xc + _Yc)
        npv = _Y / (_X + _Y)
    else:
        ppv = _X / (_X + _Y)
        npv = _Yc / (_Xc + _Yc)
    return xs, ppv, npv


def plot_roc(ax, X, Y, **kwargs):
    """Plots the ROC curve of two sample sets, cf. :func:`roc`. All keyword
    arguments are passed as is to ``matplotlib.axes.Axes.scatter``.

    Args:
        ax (matplotlib.axes.Axes): The matplotlib axes object to plot things
            on.
        X (list): The "positive" sample set (horizontal axis).
        Y (list): The "negative" sample set (vertical axis).
    """
    ax.set_xlim([-.1, 1.1])
    ax.set_ylim([-.1, 1.1])
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    _, x, y = roc(X, Y)
    ax.scatter(x, y, **kwargs)
