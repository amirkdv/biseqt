from .discovery import most_significant_shift
import os.path
from math import sqrt, ceil
from matplotlib import pyplot as plt
from .. import ProgressIndicator

# FIXME docs
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
            count = len(index.seeds(S_id, T_id))
            if count == 0:
                continue
            if set([S_id,T_id]) in true_overlaps:
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

# FIXME docs
def plot_shift_signifiance_discrimination(path, index, true_overlaps, num_bins=500, min_overlap=-1):
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
            if set([S_name, T_name]) in true_overlaps:
                pos_pvalues += [significance]
            else:
                neg_pvalues += [significance]

    indicator.finish()

    plt.clf()
    # hist returns 3 lists (n, bins, _): n is values at bins, bins is edges.
    n_neg, bins_neg, _ = plt.hist(neg_pvalues, num_bins, color='red',
        histtype='step', cumulative=True, normed=True, label='Non-overlapping reads')
    n_pos, bins_pos, _= plt.hist(pos_pvalues, num_bins, color='green',
        histtype='step', cumulative=True, normed=True, label='Overlapping reads')
    xmax = min(
        bins_neg[len(filter(lambda x: n_neg[x] < 1 - 0.001, range(len(bins_neg)-1)))],
        bins_pos[len(filter(lambda x: n_pos[x] < 1 - 0.001, range(len(bins_pos)-1)))]
    )
    plt.grid(True)
    plt.xlim(-50, xmax)
    ymax = max(max(n_neg), max(n_pos))*1.1
    plt.ylim(ymax*-0.1, ymax)
    plt.axvline(x=0, ymin=plt.ylim()[0], ymax=plt.ylim()[1], color='k')
    plt.axhline(y=0, xmin=plt.xlim()[0], xmax=plt.xlim()[1], color='k')
    x_step = 100
    y_step = 0.05
    plt.xticks([int(plt.xlim()[0]) + i*x_step for i in range(int((plt.xlim()[1]-plt.xlim()[0])/x_step) + 1)], rotation=90)
    plt.yticks([i*y_step for i in range(int(plt.ylim()[1]/y_step))])
    plt.xlabel('largest significance for a shift window on %d-mers' % index.wordlen)
    plt.ylabel('Proportion of read-pairs (cumulative)')
    plt.legend(loc='upper left', fontsize=10)
    plt.tight_layout()
    plt.savefig(path, dpi=300)

def plot_all_seeds(index, basedir='', true_overlaps=[], mappings={}, min_overlap=-1):
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
            seeds = index.seeds(S_id, T_id)
            if not seeds:
                continue
            best_shift, significance = most_significant_shift(S_id, T_id, index, min_overlap=-1)
            if significance is None:
                continue
            # FIXME debug
            #if significance >= 50 or set([S_name, T_name]) not in true_overlaps:
                #continue
            label = 'Most significant shift = %d\nlog(p-value)=%.2f' % (best_shift, significance)
            overlay = [(best_shift, '#333333', label)]

            path = os.path.join(basedir, '%s_%s' % (S_hname, T_hname))
            if true_overlaps:
                if set([S_name, T_name]) in true_overlaps:
                    color = 'green'
                    path += '.p.png'
                    true_shift = mappings[T_name[:-1]].ref_from - mappings[S_name[:-1]].ref_from
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
    plt.ylim(0, None)
    plt.xlim(0, None)
    plt.rc('text', usetex=True)
    plt.grid(True)

    S_id, T_id = seeds[0].S_id, seeds[0].T_id
    S_hname, T_hname = seqinfo[S_id]['hname'], seqinfo[T_id]['hname']
    S_len, T_len = seqinfo[S_id]['length'], seqinfo[T_id]['length']

    for shift, color, label in shift_overlay:
        xrange = (max(0, shift), S_len)
        yrange = (max(0, -shift), S_len - shift)
        plt.plot(xrange, yrange, alpha=0.2, linewidth=5, color=color, label=label)

    plt.title('\\texttt{%s} vs. \\texttt{%s} (%d seeds)' % (S_hname, T_hname, len(seeds)), fontsize=8)
    plt.axvline(x=0, ymin=plt.ylim()[0], ymax=plt.ylim()[1], color='k')
    plt.axhline(y=0, xmin=plt.xlim()[0], xmax=plt.xlim()[1], color='k')

    plt.xlabel('\\texttt{%s}' % S_hname)
    plt.ylabel('\\texttt{%s}' % T_hname)
    plt.legend(prop={'size':8}, loc='lower left')
    plt.gca().invert_yaxis()
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    plt.savefig(path, dpi=300)
