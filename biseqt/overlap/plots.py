import os.path
from math import sqrt, ceil
from matplotlib import pyplot as plt
import igraph
from scipy.stats import gaussian_kde
from scipy.integrate import cumtrapz
import numpy as np

from .discovery import most_significant_shift
from .. import ProgressIndicator
from ..mapping import Mapping

# FIXME docs
# FIXME does this still work?
def plot_num_seeds_discrimination(path, index, true_overlaps, num_bins=500, min_overlap=-1):
    plt.clf()
    seqinfo = index.seqdb.seqinfo()
    ids = seqinfo.keys()
    msg = 'Counting the number of seeds for all pairs of sequences'
    indicator = ProgressIndicator(msg, len(ids) * (len(ids)-1) / 2.0)
    indicator.start()
    pos = []
    neg = []
    for S_id_idx in range(len(ids)):
        for T_id_idx in range(S_id_idx+1, len(ids)):
            indicator.progress()
            S_id, T_id = ids[S_id_idx], ids[T_id_idx]
            S_name, T_name = seqinfo[S_id]['name'], seqinfo[T_id]['name']
            count = len(index.seeds(S_id, T_id))
            if count == 0:
                continue
            if (S_name, T_name) in true_overlaps or (T_name, S_name) in true_overlaps:
                pos += [count]
            else:
                neg += [count]

    indicator.finish()

    n_neg, bins_neg, hist_neg = plt.hist(neg, num_bins, color='red',
        histtype='step', cumulative=True, normed=True, label='Non-overlapping reads')
    n_pos, bins_pos, hist_pos = plt.hist(pos, num_bins, color='green',
        histtype='step', cumulative=True, normed=True, label='Overlapping reads')
    xmax = max(
        bins_neg[len(filter(lambda x: n_neg[x]<0.999, range(len(bins_neg)-1)))],
        bins_pos[len(filter(lambda x: n_pos[x]<0.999, range(len(bins_pos)-1)))]
    )
    plt.grid(True)
    plt.xlim(-xmax/10, xmax)
    plt.ylim(-0.1, 1.2)
    plt.axvline(x=0, ymin=-0.1, ymax=1.2, color='k')
    plt.axhline(y=0, xmin=-100, xmax=xmax, color='k')
    plt.xticks([i*100 for i in range(int(xmax/100) + 1)], rotation='vertical')
    plt.yticks([i*0.1 for i in range(11)], rotation='vertical')
    plt.tick_params(axis='x', labelsize=8, direction='vertical')
    plt.xlabel('Number of matching %d-mers' % index.wordlen)
    plt.ylabel('Proportion of read-pairs (cumulative)')
    plt.legend(loc='right')
    plt.savefig(path)

def plot_shift_consistency(path, index, true_overlaps, min_shift_significance=50, num_bins=500, min_overlap=-1):
    seqinfo = index.seqdb.seqinfo()
    ids_by_name = {seqinfo[i]['name']: i for i in seqinfo}
    errs = []
    offbysign = []
    msg = 'Finding shift error for overlapping sequences'
    indicator = ProgressIndicator(msg, len(true_overlaps))
    indicator.start()
    safety = 400 # FIXME expose
    for edge in true_overlaps:
        indicator.progress()
        S_name, T_name = tuple(edge)
        S_id, T_id = ids_by_name[S_name], ids_by_name[T_name]
        true_shift = seqinfo[S_id]['length'] - true_overlaps[edge]
        shift, sig = most_significant_shift(S_id, T_id, index, min_overlap=min_overlap)
        if abs(shift) < safety or sig < min_shift_significance:
            continue
        if shift < 0:
            offbysign += [abs(true_shift)]
        else:
            errs += [abs(shift - true_shift)]

    indicator.finish()
    errs = sorted(errs)
    offbysign = sorted(offbysign)

    plt.clf()

    density = gaussian_kde(errs)
    density.covariance_factor = lambda : .2
    markratio = 0.9
    plt.plot(errs, density(errs), antialiased=True, color='g',
        label='Absolue shift error when signs match (shaded up to %%%.0f)' % (100 * markratio))
    errs = errs[:-int(len(errs)*(1-markratio))]
    plt.fill_between(errs, density(errs), color='g', alpha=0.2)

    plt.grid(True)
    plt.xlabel('Error in estimated shift (sgn safety= %d, sig. cutoff=%.1f, word len.=%d, min. overlap=%d)' %
        (safety, min_shift_significance, index.wordlen, min_overlap))
    plt.ylabel('Density of read-pairs')
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig(path, dpi=300)


# FIXME docs
def plot_shift_pvalues(path, index, true_overlaps, num_bins=500, min_overlap=-1):
    seqinfo = index.seqdb.seqinfo()
    ids = seqinfo.keys()
    pos_pvalues = []
    neg_pvalues = []
    msg = 'Finding most significant shift for all pairs of sequences'
    indicator = ProgressIndicator(msg, len(ids) * (len(ids)-1) / 2.0)
    indicator.start()
    for S_id_idx in range(len(ids)):
        for T_id_idx in range(S_id_idx+1, len(ids)):
            S_id, T_id = ids[S_id_idx], ids[T_id_idx]
            if seqinfo[S_id]['name'][:-1] == seqinfo[T_id]['name'][:-1]:
                continue
            indicator.progress()
            S_name = seqinfo[S_id]['name']
            T_name = seqinfo[T_id]['name']
            _, significance = most_significant_shift(S_id, T_id, index, min_overlap=min_overlap)
            if significance is None:
                continue
            if (S_name, T_name) in true_overlaps or (T_name, S_name) in true_overlaps:
                pos_pvalues += [significance]
            else:
                neg_pvalues += [significance]

    indicator.finish()

    # FIXME refactor all this out (repeated in plot_shift_consistency; but
    # direction of marking is different and one is cumulative one is not)
    pos_pvalues = sorted(pos_pvalues)
    density = gaussian_kde(pos_pvalues)
    cdf = np.insert(cumtrapz(density(pos_pvalues), pos_pvalues), 0, [0])
    density.covariance_factor = lambda : .2
    markratio = 0.9
    plt.plot(pos_pvalues, cdf, antialiased=True, color='g', label='Overlapping reads')
    cdf = cdf[int(len(pos_pvalues)*(1-markratio)):]
    pos_pvalues = pos_pvalues[int(len(pos_pvalues)*(1-markratio)):]
    plt.fill_between(pos_pvalues, cdf, color='g', alpha=0.2)

    neg_pvalues = sorted(neg_pvalues)
    density = gaussian_kde(neg_pvalues)
    cdf = np.insert(cumtrapz(density(neg_pvalues), neg_pvalues), 0, [0])
    density.covariance_factor = lambda : .2
    markratio = 0.9
    plt.plot(neg_pvalues, cdf, antialiased=True, color='r', label='Non-overlapping reads')
    cdf = cdf[int(len(neg_pvalues)*(1-markratio)):]
    neg_pvalues = neg_pvalues[int(len(neg_pvalues)*(1-markratio)):]
    plt.fill_between(neg_pvalues, cdf, color='r', alpha=0.2)

    plt.grid(True)
    plt.xlabel('largest significance for a shift window on %d-mers' % index.wordlen)
    plt.ylabel('Proportion of read-pairs (cumulative)')
    plt.legend(loc='upper right', fontsize=10)
    plt.tight_layout()
    plt.savefig(path, dpi=300)

def plot_all_seeds(index, basedir='', true_overlaps={}, mappings={}, min_overlap=-1):
    seqinfo = index.seqdb.seqinfo()
    ids = seqinfo.keys()
    indicator = ProgressIndicator('Plotting all seeds',
        len(ids) * (len(ids) - 1) / 2.0, percentage=False)
    indicator.start()
    for S_id_idx in range(len(ids)):
        for T_id_idx in range(S_id_idx+1, len(ids)):
            indicator.progress()
            S_id, T_id = ids[S_id_idx], ids[T_id_idx]
            S_name, T_name = seqinfo[S_id]['name'], seqinfo[T_id]['name']
            S_hname, T_hname = seqinfo[S_id]['hname'], seqinfo[T_id]['hname']
            # FIXME debug
            if (S_name, T_name) not in true_overlaps and (T_name, S_name) not in true_overlaps:
                continue

            seeds = index.seeds(S_id, T_id)
            if not seeds:
                continue
            best_shift, significance = most_significant_shift(S_id, T_id, index, min_overlap=-1)
            if significance is None:
                continue
            label = 'Most significant shift = %d\nlog(p-value)=%.2f' % (best_shift, significance)
            overlay = [(best_shift, '#333333', label)]

            path = os.path.join(basedir, '%s_%s' % (S_hname, T_hname))
            if true_overlaps:
                # find the true shift:
                color = 'green'
                if (S_name, T_name) in true_overlaps:
                    # if S -> T then true shift for (S, T) is |S| - |overlap|.
                    true_shift = seqinfo[S_id]['length'] - true_overlaps[(S_name,T_name)]
                    overlay += [(true_shift, 'green', 'True shift = %d' % true_shift)]
                    path += '.p.png'
                elif (T_name, S_name) in true_overlaps:
                    # if T -> S then true shift for (S, T) is |overlap| - |T|
                    true_shift = - seqinfo[T_id]['length'] + true_overlaps[(T_name, S_name)]
                    overlay += [(true_shift, 'green', 'True shift = %d' % true_shift)]
                    path += '.p.png'
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
    plt.scatter([x.tx.T_idx for x in seeds], [x.tx.S_idx for x in seeds],
        marker='o', s=5, color=color)
    plt.ylim(0, None)
    plt.xlim(0, None)
    plt.rc('text', usetex=True)
    plt.grid(True)

    S_id, T_id = seeds[0].S_id, seeds[0].T_id
    S_hname, T_hname = seqinfo[S_id]['hname'], seqinfo[T_id]['hname']
    S_len, T_len = seqinfo[S_id]['length'], seqinfo[T_id]['length']

    for shift, color, label in shift_overlay:
        _yrange = (max(0, shift), S_len)
        _xrange = (max(0, -shift), S_len - shift)
        plt.plot(_xrange, _yrange, alpha=0.2, linewidth=5, color=color, label=label)

    plt.title('\\texttt{%s} vs. \\texttt{%s} (%d seeds)' % (S_hname, T_hname, len(seeds)), fontsize=8)
    plt.axvline(x=0, ymin=plt.ylim()[0], ymax=plt.ylim()[1], color='k')
    plt.axhline(y=0, xmin=plt.xlim()[0], xmax=plt.xlim()[1], color='k')

    plt.xlabel('\\texttt{%s}' % T_hname)
    plt.ylabel('\\texttt{%s}' % S_hname)
    plt.legend(prop={'size':8}, loc='lower left')
    plt.gca().invert_yaxis()
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    plt.savefig(path, dpi=300)
