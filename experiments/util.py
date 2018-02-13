import pickle
import os
import sys
import numpy as np
from matplotlib import pyplot as plt


# TODO change syntax to pickle_dump(obj, path)
def pickle_dump(path, obj, comment=None):
    try:
        with open(path, 'w') as f:
            pickle.dump(obj, f)
    except IOError, e:
        sys.stderr.write('Error writing to to %s: %s\n' % (path, str(e)))
        n = np.random.choice(list('0123456789abcdef'), size=10)
        path = '/tmp/' + ''.join(str(i) for i in n)
        with open(path, 'w') as f:
            pickle.dump(obj, f)
    sys.stderr.write('dumped %s to %s.\n' %
                     (comment if comment else 'object', path))


def pickle_load(path):
    with open(path, 'rb') as f:
        return pickle.load(f)


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


def savefig(fig, path, comment=None):
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(path, dpi=270)
    sys.stderr.write('recreated plot%s at %s\n' %
                     ('' if comment is None else ' "%s"' % comment, path))
