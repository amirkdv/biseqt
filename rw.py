from matplotlib import pyplot as plt
from math import sqrt,ceil
from align import ProgressIndicator
import igraph

G = igraph.read('leishmania_true.gml')
db_id_from_graph_id = lambda vid: int(G.vs[vid]['name'].split('#')[1])
true_overlaps = [set([db_id_from_graph_id(u), db_id_from_graph_id(v)]) for u,v in G.get_edgelist()]
vertices = {db_id_from_graph_id(int(v['id'])):v['name'] for v in G.vs}
start_pos_from_db_id = lambda dbid: int(vertices[dbid].split('#')[0].split()[1].split('-')[0])

max_plots = 225
min_rw_length = 5
draw_type = ['-', '+']

with open('scores.txt') as f:
    data = [l.strip().split() for l in f.readlines() if l.strip()[-1] in draw_type]

data = data[:max_plots]
data = [[d[0], eval(d[1]), d[2]] for d in data]
data = filter(lambda x: len(x[1]) >= min_rw_length, data)
indicator = ProgressIndicator('plotting score random walks', len(data))
indicator.start()
dim = ceil(sqrt(len(data)))
fig = plt.figure(figsize=(5*dim,5*dim))
for idx,datum in enumerate(data):
    S_tok, T_tok = datum[0][1:-1].split(',')
    S_id, S_idx = [int(i) for i in S_tok.split(':')]
    T_id, T_idx = [int(i) for i in T_tok.split(':')]
    ax = fig.add_subplot(dim, dim, idx+1)
    if set([S_id, T_id]) in true_overlaps:
        color = 'green'
        true_shift = start_pos_from_db_id(T_id) - start_pos_from_db_id(S_id)
        ax.set_title('%d, %d (%d)' % (S_id, T_id, true_shift))
    else:
        ax.set_title('%d, %d' % (S_id, T_id))
        color = 'red'

    ax.plot([x*50 for x in range(len(datum[1]))], datum[1], color=color, label=' '.join([str(S_idx-T_idx), datum[2]]))
    ax.legend()
    #ax.set_xticks([])
    #ax.set_yticks([])
    indicator.progress()

indicator.finish()

plt.savefig('rw.png')
