import sqlite3
import igraph
import os
from matplotlib import pyplot as plt
from align import seq, tuples, ProgressIndicator


db = 'genome.leishmania.hp_assembly.db'
wordlen = int(os.environ['WORDLEN'])
base_dir = 'spectra-%d' % wordlen
pos_dir = os.path.join(base_dir, 'positive')
neg_dir = os.path.join(base_dir, 'negative')
shift_distrib_window = 500
min_threshold_coeff = 6
max_threshold_coeff = 15
min_ymax = 300
num_bins = 500
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

pos_mode_ratios = []
neg_mode_ratios = []
with sqlite3.connect(db) as conn:
    c = conn.cursor()
    c.execute('SELECT id, LENGTH(seq) FROM seq ORDER BY id ASC')
    lengths_by_id = {row[0]: row[1] for row in c}
    ids = lengths_by_id.keys()
    #os.mkdir(base_dir)
    #os.mkdir(pos_dir)
    #os.mkdir(neg_dir)
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
            uniform_prob = float(sum(shifts))/len(shifts)

            #plt.clf()
            #color = 'red'
            if set([S_id, T_id]) in true_overlaps:
                #color = 'green'
                #true_shift = start_pos_from_db_id(T_id) - start_pos_from_db_id(S_id)
                #plt.axvline(x=true_shift, color='#333333', linewidth=20, alpha=0.4)
                pos_mode_ratios += [max(shifts)/uniform_prob]
            else:
                neg_mode_ratios += [max(shifts)/uniform_prob]

            #plt.axhline(y=min_threshold_coeff*uniform_prob, xmin=drange[0], xmax=drange[-1], linestyle='--', color='k')
            #plt.axhline(y=max_threshold_coeff*uniform_prob, xmin=drange[0], xmax=drange[-1], linestyle='--', color='k')
            #plt.plot(drange, shifts, color=color)
            #plt.tick_params(axis='both', which='major', labelsize=8)
            #plt.ylim(0, max(min_ymax, max(max_threshold_coeff*uniform_prob, max(shifts))) + 100)
            #plt.grid(True)
            #plt.title('%s ----> %s (%d total seeds)' % (vertices[S_id], vertices[T_id], len(seeds)), fontsize=8, fontname='Inconsolata')
            #plt.savefig(os.path.join(
                #pos_dir if color == 'green' else neg_dir,
                #'%d_%d.png' % (S_id, T_id)
            #))

    indic.finish()

    plt.rc('font', family='Cardo')
    plt.rc('text', usetex=True) # requires dvipng installed
    n_neg, bins_neg, hist_neg = plt.hist(neg_mode_ratios, num_bins, color='red',
        histtype='step', cumulative=True, normed=True, label='Non-overlapping reads')
    n_pos, bins_pos, hist_pos = plt.hist(pos_mode_ratios, num_bins, color='green',
        histtype='step', cumulative=True, normed=True, label='Overlapping reads')
    xmax = max(
        bins_neg[len(filter(lambda x: n_neg[x]<0.999, range(len(bins_neg)-1)))],
        bins_pos[len(filter(lambda x: n_pos[x]<0.9, range(len(bins_pos)-1)))]
    )
    plt.grid(True)
    plt.xlim(-xmax/10, xmax)
    plt.ylim(-0.1, 1.2)
    plt.axvline(x=0, ymin=-0.1, ymax=1.2, color='k')
    plt.axhline(y=0, xmin=-100, xmax=xmax, color='k')
    plt.xticks([i*5 for i in range(int(xmax/5) + 1)], rotation='vertical')
    plt.yticks([i*0.1 for i in range(11)], rotation='vertical')
    plt.xlabel(r'$\frac{\mathrm{max. freq.}}{\mathrm{unif. freq.}}$ of %d-mer shifts' % wordlen)
    plt.ylabel('Proportion of read-pairs (cumulative)')
    plt.legend(loc='upper left')
    plt.savefig('ratios.%d.png' % wordlen)
