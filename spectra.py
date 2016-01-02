import sqlite3
import igraph
import os
from matplotlib import pyplot as plt
from align import seq, tuples, ProgressIndicator
from Bio import Seq, SeqIO, SeqRecord

plt.rc('font', family='Cardo')


db = 'genome.leishmania.hp_assembly.db'
wordlen = int(os.environ['WORDLEN'])
base_dir = 'spectra-%d' % wordlen
pos_dir = os.path.join(base_dir, 'positive')
neg_dir = os.path.join(base_dir, 'negative')
shift_rolling_sum_width = 100
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
        yield sum(src)
        return
    cur = 0
    for i in range(0, len(src)):
        if i >= width:
            cur -= src[i-width]
        if i < len(src):
            cur += src[i]
        yield cur

from math import sqrt,log

def _shift_log_pvalue(S_len, T_len, shift, num, num_reads):
    if num == 0:
        return 0
    L = shift_rolling_sum_width
    log_prob = -log(S_len) - log(T_len) + log(L) + 0.5 * log(
        (S_len - abs(shift))**2 + (T_len - abs(shift))**2
    )
    # 1- we have num observations (each a seed) with the same probability
    #    of being matched by the null hypothesis
    # 2- we are testing S_len+T_len simultaneous hypotheses;
    #    apply a Bonferroni correction:
    return 2*log(num_reads) + log(S_len + T_len) + num * log_prob

pos_pvalues = []
neg_pvalues = []

pos_mode_ratios = []
neg_mode_ratios = []
with sqlite3.connect(db) as conn:
    c = conn.cursor()
    c.execute('SELECT id, LENGTH(seq) FROM seq ORDER BY id ASC')
    lengths_by_id = {row[0]: row[1] for row in c}
    num_reads = len(lengths_by_id)
    ids = lengths_by_id.keys()
    #ids = [1]
    if not os.path.exists(base_dir):
        os.mkdir(base_dir)
    if not os.path.exists(pos_dir):
        os.mkdir(pos_dir)
    if not os.path.exists(neg_dir):
        os.mkdir(neg_dir)
    indic = ProgressIndicator('building spectra', len(ids)*(len(ids)-1) / 2, percentage=False)
    indic.start()

    #annotated = []
    for S_id_idx in range(len(ids)):
        for T_id_idx in range(S_id_idx+1, len(ids)):
            S_id, T_id = ids[S_id_idx], ids[T_id_idx]
            S_len, T_len = lengths_by_id[S_id], lengths_by_id[T_id]
            #if S_id != 1:
                #continue

            seeds = I.seeds(S_id, T_id)
            if not seeds:
                continue
            indic.progress()
            drange = range(-T_len, S_len)
            shifts = {d:[] for d in drange}
            for seed in seeds:
                d = seed.tx.S_idx - seed.tx.T_idx
                shifts[d] += [seed]
            shift_distrib = [x for x in rolling_sum([len(x[1]) for x in sorted(shifts.items())], shift_rolling_sum_width)]

            plt.clf()
            max_idx, num = max(
                enumerate(shift_distrib),
                key=lambda x: x[1]
            )
            best_shift = drange[min(max_idx, len(drange)-1)]

            #c.execute('select name, description, seq from seq where id = ?', (T_id,))
            #name, descr, seqstr = c.next()
            #assert(set(seqstr).issubset('ACGT'))
            #newid = '%s_P%d' % (name.split('_')[0], best_shift)
            #annotated += [SeqRecord.SeqRecord(
                #Seq.Seq(seqstr), id=newid, description='%s %s' % (newid, descr.split()[1])
            #)]
            # seeds = sum([shifts[d] for d in range(best_shift-shift_rolling_sum_width,best_shift+1)], [])
            # S_idx_range_min = min([s.tx.S_idx for s in seeds])
            # S_idx_range_max = max([s.tx.S_idx for s in seeds]) + 1
            # S_idx_range = range(S_idx_range_min, S_idx_range_max)
            # shifts_in_slice = {idx:[] for idx in S_idx_range}
            # for seed in seeds:
            #     shifts_in_slice[seed.tx.S_idx] += [seed]
            # in_slice_rolling_sum_width = 50 * shift_rolling_sum_width
            # shifts_in_slice_distrib = [
            #     x for x in rolling_sum(
            #         [len(x[1]) for x in sorted(shifts_in_slice.items())], in_slice_rolling_sum_width
            #     )
            # ]
            # max_idx, num = max(
            #     enumerate(shifts_in_slice_distrib),
            #     key=lambda x: x[1]
            # )
            # best_S_idx = S_idx_range[max_idx]
            log_pvalue = _shift_log_pvalue(S_len, T_len, best_shift, num, num_reads)

            color = 'red'
            if set([S_id, T_id]) in true_overlaps:
                color = 'green'
                pos_pvalues += [log_pvalue]
                true_shift = start_pos_from_db_id(T_id) - start_pos_from_db_id(S_id)
            #     plt.axvline(x=true_shift, color='#333333', linewidth=20, alpha=0.4)
            #     pos_mode_ratios += [max(shifts)/uniform_prob]
            else:
            #     neg_mode_ratios += [max(shifts)/uniform_prob]
                neg_pvalues += [log_pvalue]

            plt.axhline(y=0, xmin=-1000, xmax=1000, color='k')
            plt.axvline(x=0, ymin=-1000, ymax=1000, color='k')
            plt.tick_params(axis='both', which='major', labelsize=8)

            #plt.gca().set_aspect('equal')
            #plt.scatter([x.tx.S_idx for x in seeds], [x.tx.T_idx for x in seeds], marker='o', s=5, color=color)
            #plt.ylim(-1000)
            #plt.xlim(-1000)
            #plt.grid(True)

            #def plot_shift(shift, color, label):
                #xrange = (max(0, shift), S_len)
                #yrange = (max(0, -shift), S_len - shift)
                #plt.plot(xrange, yrange, alpha=0.4, linewidth=5, color=color, label=label)

            #plot_shift(best_shift, color, 'Lowest p-value shift = %d\nlog(p-value)=%.2f' % (best_shift, log_pvalue))
            #if color == 'green':
                #plot_shift(true_shift, '#333333', 'True shift = %d' % true_shift)

            #plt.title('%s vs. %s\n%d total seeds' % (vertices[S_id], vertices[T_id], len(seeds)), fontsize=8)
            #plt.axvline(x=0, ymin=plt.ylim()[0], ymax=plt.ylim()[1], color='k')
            #plt.axhline(y=0, xmin=plt.xlim()[0], xmax=plt.xlim()[1], color='k')

            #plt.xlabel('Position in %s' % vertices[S_id].split()[-1])
            #plt.ylabel('Position in %s' % vertices[T_id].split()[-1])
            #plt.legend(prop={'size':10})
            #plt.savefig(os.path.join(
                #pos_dir if color == 'green' else neg_dir,
                #'%d.%d.png' % (S_id,T_id)
            #))

    indic.finish()
    #SeqIO.write(annotated, 'reads.annotated.fixed.fa', 'fasta')

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
    plt.xlabel('smallest (Bonferroni-corrected) log(p-values) for a shift window on %d-mers' % wordlen)
    plt.ylabel('Proportion of read-pairs (cumulative)')
    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.savefig('pvalues.%d.png' % wordlen)
