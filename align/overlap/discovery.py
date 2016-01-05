from collections import namedtuple
from matplotlib import pyplot as plt
from math import log, sqrt, ceil
import os.path
from .. import pw, lib, ffi, ProgressIndicator

# TODO make this a C struct so the C code cleans up.
class SeedExtensionParams(namedtuple('SeedExtensionParams',
    ['window', 'min_overlap_score', 'max_new_mins', 'align_params'])):
    """Represents the set of tuning parameters for seed extension.

    Attributes:
        window (int): The size of the rolling window.
        min_overlap_score (float): The minimum required score for an alignment
            to be reported
        max_new_mins (int): Maximum number of new minima observed in the score
            random walk until a segement is dropped.
        align_params (pw.AlignParams): The alignment parameters for the
            rolling alignment.
    """
    pass

def extend_segments(S, T, segments, params, rw_collect=False):
    """Wraps :c:func:`extend()`: given two sequences and a number of
    matching segments returns the first fully extended segment.

    Args:
        S (seq.Sequence): The "from" sequence.
        T (seq.Sequence): The "to" sequence.
        segments (list[pw.Segment]): The starting segments. If called
            from :func:`build`, these are seeds but no assumption is made.
        params (SeedExtensionParams): The parameters for seed extension.

    Returns:
        pw.Segment|None: A segment corresponding to an overlap alignment or
            None if none found.
    """
    if not segments:
        return None

    segs = ffi.new('segment* []', [seg.c_obj for seg in segments])
    res = lib.extend(
        segs, len(segs),
        S.c_idxseq, T.c_idxseq, len(S), len(T), params.align_params.c_obj,
        params.window, params.max_new_mins, params.min_overlap_score, int(bool(rw_collect))
    )
    return pw.Segment(c_obj=res) if res != ffi.NULL else None

# FIXME docs
def plot_seed_extension_rws(path, seqinfo, max_rws=225, draw_type='-+',
    logfile='scores.txt', true_overlaps=[]):
    with open(logfile) as f:
        data = [l.strip().split() for l in f.readlines() if l.strip()[-1] in draw_type]

    data = data[:max_rws]
    data = [[d[0], eval(d[1]), d[2]] for d in data]
    indicator = ProgressIndicator('plotting score random walks', len(data))
    indicator.start()
    dim = ceil(sqrt(len(data)))

    plt.clf()
    fig = plt.figure(figsize=(5*dim,5*dim))
    for idx, datum in enumerate(data):
        S_tok, T_tok = datum[0][1:-1].split(',')
        S_id, S_idx = [int(i) for i in S_tok.split(':')]
        T_id, T_idx = [int(i) for i in T_tok.split(':')]
        S_start, T_start = seqinfo[S_id]['start'], seqinfo[T_id]['start']
        ax = fig.add_subplot(dim, dim, idx+1)
        if set([S_id, T_id]) in true_overlaps:
            color = 'green'
            true_shift = T_start - S_start
            ax.set_title('%d, %d (%d)' % (S_id, T_id, true_shift))
        else:
            ax.set_title('%d, %d' % (S_id, T_id))
            color = 'red'

        xs = [x*50 for x in range(len(datum[1]))]
        label = ' '.join([str(S_idx - T_idx), datum[2]])
        ax.plot(xs, datum[1], color=color, label=label)
        ax.legend()
        indicator.progress()

    indicator.finish()
    plt.savefig(path)


# FIXME docs
def rolling_sum(data, width):
    if len(data) < width:
        return
    cur = 0
    for idx in range(len(data)):
        if idx >= width:
            cur -= data[idx - width]
        if idx < len(data):
            cur += data[idx]
        yield cur

# FIXME docs
class ShiftWindow(namedtuple('ShiftWindow', ['S_len', 'T_len', 'width', 'shift'])):
    # FIXME docs
    def significance(self, num):
        diag_len_squared = (self.S_len - abs(self.shift))**2 + \
            (self.T_len - abs(self.shift))**2
        log_pvalue = - log(self.S_len) - log(self.T_len) \
            + log(self.width) \
            + 0.5 * log(diag_len_squared)
        # we have num observations (each a seed) with the same probability
        # of being matched by the null hypothesis:
        log_pvalue *= num
        # we are testing S_len+T_len simultaneous hypotheses; apply a Bonferroni
        # correction:
        log_pvalue += log(self.S_len + self.T_len)
        return log_pvalue

# FIXME find the strip with max number of shifts *per unit area*.
def most_signitifcant_shift(S_len, T_len, seeds, rolling_sum_width):
    """Builds a smoothed distribution (via a rolling sum of known width) of
    shifts for the provided seeds of a pair of sequences of known lengths and
    returns the shift with most number of seeds and its p-value.

    The shift window corresponding to shift :math:`s` is the interval
    :math:`[s-w, s]` of shift values where :math:`w` is the rolling sum width.

    The p-value for a given shift is calculated as follows: consider the
    distribution of seeds in the :math:`(i_S,i_T)` plane where :math:`i_S,i_T`
    are the starting coordinates of each seed. In this plane, all seeds reside
    within the rectangle :math:`M = [0, |S|] \\times [0, |T|]`. We assume that
    under the null hypothesis (the two sequences are not overlapping) seeds
    are uniformly distributed over :math:`M`. Further, in this plane a shift
    window :math:`[s-w, s]` corresponds to a diagonal strip at horizontal
    and vertical distance :math:`s` of the main quadrant diagonal. The p-value
    is then calculated as the probability that as many seeds would be observed
    in the corresponding strip:

        :math:`\\Pr(X_s\\ge n) \\simeq \\left(\\frac{A_s}{|S|\cdot|T|}\\right)^n`

    where :math:`A_s` is the area of the diagonal :math:`w`-strip corresponding
    to shift :math:`s`:

        :math:`A_s = w\\sqrt{(|S|-|s|)^2 + (|T|-|s|)^2}`

    The final formula for the returned pvalue logarithm is:

        :math:`\\log(|S|+|T|) + n \\left[\\log(A_s) - \\log(|S|) - \\log(|T|)\\right]`

    where the first term is a Bonferroni correction since the process involves
    testing the same hypothesis for :math:`|S|+|T|` different strips of shift
    windows.

    Args:
        seeds (list[Segment]): Exactly matching seeds.
        S_len (int): The length of the "from" sequence.
        T_len (int): The length of the "to" sequence.
        rolling_sum_width(int): The parameter :math:`w` in the above.

    Returns:
        (mode_shift, log_pvalue): A tuple of the form ``(int, float)``.

    """
    shift_range = range(-T_len, S_len)
    shift_coverage = {shift:0 for shift in shift_range}
    for seed in seeds:
        # all seeds are the same length at this stage (and they
        # are potentially overlapping):
        shift_coverage[seed.tx.S_idx - seed.tx.T_idx] += 1

    shift_coverage = [x[1] for x in sorted(shift_coverage.items())]
    shift_distrib = rolling_sum(shift_coverage, rolling_sum_width)
    mode_idx, mode = max(enumerate(shift_distrib), key=lambda x: x[1])
    most_dense_window = ShiftWindow(S_len=S_len, T_len=T_len,
        width=rolling_sum_width,
        shift=shift_range[min(mode_idx, len(shift_range)-1)]
    )

    return most_dense_window.shift, most_dense_window.significance(mode)

# FIXME docs
def plot_shift_signifiance_discrimination(path, index, rolling_sum_width,
    true_overlaps, num_bins=500):
    seqinfo = index.tuplesdb.seqinfo()
    ids = seqinfo.keys()
    pos_pvalues = []
    neg_pvalues = []
    msg = 'Finding most significant shift for all pairs of sequences'
    indicator = ProgressIndicator(msg,
        len(seqinfo) * (len(seqinfo)-1) / 2.0, percentage=False)
    indicator.start()
    for S_id_idx in range(len(ids)):
        for T_id_idx in range(S_id_idx+1, len(ids)):
            indicator.progress()
            S_id, T_id = ids[S_id_idx], ids[T_id_idx]
            seeds = index.seeds(S_id, T_id)
            if not seeds:
                continue
            S_len, T_len = seqinfo[S_id]['length'], seqinfo[T_id]['length']
            _, log_pvalue = most_signitifcant_shift(S_len, T_len,
                seeds, rolling_sum_width)
            if set([S_id, T_id]) in true_overlaps:
                pos_pvalues += [log_pvalue]
            else:
                neg_pvalues += [log_pvalue]

    indicator.finish()

    plt.clf()
    n_neg, bins_neg, hist_neg = plt.hist(neg_pvalues, num_bins, color='red',
        histtype='step', cumulative=True, normed=True, label='Non-overlapping reads')
    n_pos, bins_pos, hist_pos = plt.hist(pos_pvalues, num_bins, color='green',
        histtype='step', cumulative=True, normed=True, label='Overlapping reads')
    xmin = max(
        bins_neg[len(filter(lambda x: n_neg[x] < 0.005, range(len(bins_neg)-1)))],
        bins_pos[len(filter(lambda x: n_pos[x] < 0.005, range(len(bins_pos)-1)))]
    )
    xmax = int(min(bins_neg)*-0.1)
    plt.grid(True)
    plt.xlim(xmin, xmax)
    plt.ylim(-0.1, 1.2)
    plt.axvline(x=0, ymin=-0.1, ymax=1.2, color='k')
    plt.axhline(y=0, xmin=xmin, xmax=xmax, color='k')
    plt.xticks([int(xmin) + i*100 for i in range(int(abs(xmin)/100) + 1)])
    plt.yticks([i*0.05 for i in range(21)])
    plt.xlabel('smallest log(p-values) for a shift window on %d-mers' % index.wordlen)
    plt.ylabel('Proportion of read-pairs (cumulative)')
    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.savefig(path, dpi=300)

def plot_all_seeds(index, rolling_sum_width, basedir='', true_overlaps=[]):
    seqinfo = index.tuplesdb.seqinfo()
    ids = seqinfo.keys()
    indicator = ProgressIndicator('Plotting all seeds',
        len(ids) * (len(ids) - 1) / 2.0, percentage=False)
    indicator.start()
    for S_id_idx in range(len(ids)):
        for T_id_idx in range(S_id_idx+1, len(ids)):
            indicator.progress()
            S_id, T_id = ids[S_id_idx], ids[T_id_idx]
            S_start, T_start = seqinfo[S_id]['start'], seqinfo[T_id]['start']
            S_len, T_len = seqinfo[S_id]['length'], seqinfo[T_id]['length']
            seeds = index.seeds(S_id, T_id)
            if not seeds:
                continue
            best_shift, log_pvalue = most_signitifcant_shift(S_len, T_len,
                seeds, rolling_sum_width)
            label = 'Most significant shift = %d\nlog(p-value)=%.2f' % (best_shift, log_pvalue)
            overlay = [(best_shift, '#333333', label)]

            path = os.path.join(basedir, '%d_%d' % (S_id, T_id))
            if true_overlaps:
                if set([S_id, T_id]) in true_overlaps:
                    color = 'green'
                    path += '.p.png'
                    true_shift = T_start - S_start
                    overlay += [(true_shift, 'green', 'True shift = %d' % true_shift)]
                else:
                    path += '.n.png'
                    color = 'red'

                plot_seeds(path, seeds, seqinfo, color=color, shift_overlay=overlay)
            else:
                path += '.png'
                plot_seeds(path, seeds, seqinfo, shift_overlay=overlay)

    indicator.finish()

def plot_seeds(path, seeds, seqinfo, color='k', shift_overlay=[]):
    plt.clf()
    plt.gca().set_aspect('equal')
    plt.scatter([x.tx.S_idx for x in seeds], [x.tx.T_idx for x in seeds],
        marker='o', s=5, color=color)
    plt.ylim(-1000)
    plt.xlim(-1000)
    plt.grid(True)

    S_id, T_id = seeds[0].S_id, seeds[0].T_id
    S_name, T_name = seqinfo[S_id]['name'], seqinfo[T_id]['name']
    S_len, T_len = seqinfo[S_id]['length'], seqinfo[T_id]['length']

    for shift, color, label in shift_overlay:
        xrange = (max(0, shift), S_len)
        yrange = (max(0, -shift), S_len - shift)
        plt.plot(xrange, yrange, alpha=0.4, linewidth=5, color=color, label=label)

    plt.title('%s vs. %s\n%d total seeds' % (S_name, T_name, len(seeds)), fontsize=8)
    plt.axvline(x=0, ymin=plt.ylim()[0], ymax=plt.ylim()[1], color='k')
    plt.axhline(y=0, xmin=plt.xlim()[0], xmax=plt.xlim()[1], color='k')

    plt.xlabel('Position in %s' % S_name)
    plt.ylabel('Position in %s' % T_name)
    plt.legend(prop={'size':8})
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    plt.savefig(path, dpi=300)
