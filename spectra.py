import sqlite3
import igraph
import os
from matplotlib import pyplot as plt

db = 'genome.leishmania.hp_assembly.db'
wordlen = 10
base_dir = 'spectra-%d' % wordlen
pos_dir = os.path.join(base_dir, 'positive')
neg_dir = os.path.join(base_dir, 'negative')
num_bins = 1000
ylim = (0,100)
xlim = (-15000,15000)
G = igraph.read('leishmania_true.gml')
db_id_from_graph_id = lambda vid: int(G.vs[vid]['name'].split('#')[1])
true_overlaps = [set([db_id_from_graph_id(u), db_id_from_graph_id(v)]) for u,v in G.get_edgelist()]
vertices = {db_id_from_graph_id(int(v['id'])):v['name'] for v in G.vs}
start_pos_from_db_id = lambda dbid: int(vertices[dbid].split('#')[0].split()[1].split('-')[0])

red = '#ffe7e7'
green = '#6bdb6b'

from align import ProgressIndicator

with sqlite3.connect(db) as conn:
    c = conn.cursor()
    c.execute('SELECT id FROM seq ORDER BY id ASC')
    ids = [row[0] for row in c]
    N = len(ids)*(len(ids)-1) / 2
    os.mkdir(base_dir)
    os.mkdir(pos_dir)
    os.mkdir(neg_dir)
    indic = ProgressIndicator('building spectra', N, percentage=False)
    indic.start()
    for S_id_idx in range(len(ids)):
        for T_id_idx in range(S_id_idx, len(ids)):
            S_id = ids[S_id_idx]
            T_id = ids[T_id_idx]
            q = """
                SELECT S_idx - T_idx FROM seeds_%d
                WHERE S_id = %d AND T_id = %d
            """ % (wordlen, S_id, T_id)
            c.execute(q)
            shifts = [row[0] for row in c]
            if len(shifts) < 5:
                continue
            indic.progress()

            plt.clf()
            color = 'red'
            if set([S_id, T_id]) in true_overlaps:
                color = 'green'
                true_shift = start_pos_from_db_id(T_id) - start_pos_from_db_id(S_id)
                plt.axvline(x=true_shift, ymin=ylim[0], ymax=ylim[1], color='#333333', linewidth=20, alpha=0.4)
            plt.hist(shifts, num_bins, histtype='stepfilled', color=color, edgecolor=color)
            plt.tick_params(axis='both', which='major', labelsize=8)

            plt.xlim(*xlim)
            plt.ylim(*ylim)
            plt.grid(True)
            plt.title('%s ----> %s (%d total seeds)' % (vertices[S_id], vertices[T_id], len(shifts)), fontsize=8, fontname='Inconsolata')
            plt.savefig(os.path.join(
                pos_dir if color == 'green' else neg_dir,
                '%d_%d.png' % (S_id, T_id)
            ))

    indic.finish()
