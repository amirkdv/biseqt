import sqlite3
import igraph
import os
from matplotlib import pyplot as plt
from align import seq, tuples, ProgressIndicator

plt.rc('font', family='Cardo')


db = 'genome.leishmania.hp_assembly.db'
wordlen = int(os.environ['WORDLEN'])
base_dir = 'spectra-%d' % wordlen
pos_dir = os.path.join(base_dir, 'positive')
neg_dir = os.path.join(base_dir, 'negative')
shift_distrib_window = 40
min_threshold_coeff = 6
max_threshold_coeff = 15
min_ymax = 300
num_bins = 100
G = igraph.read('leishmania_true.gml')
db_id_from_graph_id = lambda vid: int(G.vs[vid]['name'].split('#')[1])
true_overlaps = [set([db_id_from_graph_id(u), db_id_from_graph_id(v)]) for u,v in G.get_edgelist()]
vertices = {db_id_from_graph_id(int(v['id'])):v['name'] for v in G.vs}
start_pos_from_db_id = lambda dbid: int(vertices[dbid].split('#')[0].split()[1].split('-')[0])

red = '#ffe7e7'
green = '#6bdb6b'

B = tuples.TuplesDB(db, alphabet=seq.Alphabet('ACGT'))
I = tuples.Index(B, wordlen=wordlen)
def rolling_sum(src, width):
    if len(src) < width:
        return
    cur = 0
    for i in range(0, len(src)):
        if i >= width:
            cur -= src[i-width]
        if i < len(src):
            cur += src[i]
        yield cur

from math import sqrt,log
pvalue = lambda L, S_len, T_len, shift, num: log(2*min(S_len,T_len)) + num*(-log(S_len) - log(T_len) + log(L) + 0.5*log(((S_len-abs(shift))**2 + (T_len-abs(shift))**2)/2))

pos_pvalues = []
neg_pvalues = []

pos_mode_ratios = []
neg_mode_ratios = []
with sqlite3.connect(db) as conn:
    c = conn.cursor()
    c.execute('SELECT id, LENGTH(seq) FROM seq ORDER BY id ASC')
    lengths_by_id = {row[0]: row[1] for row in c}
    ids = lengths_by_id.keys()
    # ids = [3,7]
    if not os.path.exists(base_dir):
        os.mkdir(base_dir)
    if not os.path.exists(pos_dir):
        os.mkdir(pos_dir)
    if not os.path.exists(neg_dir):
        os.mkdir(neg_dir)
    indic = ProgressIndicator('building spectra', len(ids)*(len(ids)-1) / 2, percentage=False)
    indic.start()
    for S_id_idx in range(len(ids)):
        for T_id_idx in range(S_id_idx, len(ids)):
            S_id, T_id = ids[S_id_idx], ids[T_id_idx]
            S_len, T_len = lengths_by_id[S_id], lengths_by_id[T_id]
            seeds = I.seeds(S_id, T_id)
            if not seeds:
                continue
            indic.progress()
            drange = range(-T_len, S_len)
            shifts = {d:0 for d in drange}
            for seed in seeds:
                d = seed.tx.S_idx - seed.tx.T_idx
                shifts[d] += len(seed.tx.opseq)
            shifts = [x for x in rolling_sum([x[1] for x in sorted(shifts.items())], shift_distrib_window)]
            #uniform_prob = float(sum(shifts))/len(shifts)

            plt.clf()
            mode_idx, mode = max(enumerate(shifts), key=lambda x: x[1])
            mode_shift = drange[min(mode_idx, len(drange)-1)]

            color = 'red'
            if set([S_id, T_id]) in true_overlaps:
                color = 'green'
                pos_pvalues += [pvalue(shift_distrib_window, S_len, T_len, mode_shift, mode)]
                true_shift = start_pos_from_db_id(T_id) - start_pos_from_db_id(S_id)
            #     plt.axvline(x=true_shift, color='#333333', linewidth=20, alpha=0.4)
            #     pos_mode_ratios += [max(shifts)/uniform_prob]
            else:
            #     neg_mode_ratios += [max(shifts)/uniform_prob]
                neg_pvalues += [pvalue(shift_distrib_window, S_len, T_len, mode_shift, mode)]
            #print '\n%s %s' % ('-' if color == 'red' else '+', (mode, mode_shift))

            # plt.axhline(y=0, xmin=-1000, xmax=1000, color='k')
            # plt.axvline(x=0, ymin=-1000, ymax=1000, color='k')
            # plt.axhline(y=max_threshold_coeff*uniform_prob, xmin=drange[0], xmax=drange[-1], linestyle='--', color='k')
            # plt.plot(drange, shifts, color=color)
            # plt.tick_params(axis='both', which='major', labelsize=8)
            # plt.ylim(0, max(min_ymax, max(max_threshold_coeff*uniform_prob, max(shifts))) + 100)


            plt.scatter([x.tx.S_idx for x in seeds], [x.tx.T_idx for x in seeds], marker='o', s=5, color=color)
            plt.gca().set_aspect('equal')
            plt.ylim(-1000)
            plt.xlim(-1000)
            plt.grid(True)
            if color == 'green':
                plt.plot((0, S_len), (-true_shift, S_len-true_shift), linewidth=10, alpha=0.4, color='#333333')
            plt.title('%s ----> %s (%d total seeds)' % (vertices[S_id], vertices[T_id], len(seeds)), fontsize=8, fontname='Inconsolata')
            plt.savefig(os.path.join(
                pos_dir if color == 'green' else neg_dir,
                '%d_%d.png' % (S_id, T_id)
            ))

    indic.finish()

    n_neg, bins_neg, hist_neg = plt.hist(neg_pvalues, num_bins, color='red',
        histtype='step', cumulative=True, normed=True, label='Non-overlapping reads')
    n_pos, bins_pos, hist_pos = plt.hist(pos_pvalues, num_bins, color='green',
        histtype='step', cumulative=True, normed=True, label='Overlapping reads')
    xmin,xmax = -10000, 1000
    plt.grid(True)
    plt.xlim(xmin, xmax)
    plt.ylim(-0.1, 1.2)
    plt.axvline(x=0, ymin=-0.1, ymax=1.2, color='k')
    plt.axhline(y=0, xmin=xmin, xmax=xmax, color='k')
    plt.xticks([xmin + i*1000 for i in range(int(abs(xmin)/1000) + 1)])
    plt.yticks([i*0.05 for i in range(21)])
    #plt.xlabel(r'$\frac{\mathrm{max. freq.}}{\mathrm{unif. freq.}}$ of %d-mer shifts' % wordlen)
    #plt.ylabel('Proportion of read-pairs (cumulative)')
    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.savefig('pvalues.%d.png' % wordlen)
