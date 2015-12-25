import sqlite3
from itertools import product
from matplotlib import pyplot as plt
import igraph
import os

wordlen = int(os.environ['WORDLEN'])
num_bins = 500
max_seed_count = 3000
min_seeds_for_homology = 140 # just marks the appropriate line in the graph
db = 'genome.leishmania.hp_assembly.db'
G = igraph.read('leishmania_true.gml')
name_to_id = lambda name: int(name.split()[-1][1:])
def _endpoint_ids(eid):
    if isinstance(eid, igraph.Edge):
        eid = eid.index
    uid, vid = G.es[eid].tuple
    return name_to_id(G.vs[uid]['name']), name_to_id(G.vs[vid]['name'])

sE = set([_endpoint_ids(e) for e in G.es])
q = 'select count(*) from seeds_%d where S_id = ? and T_id = ?;' % wordlen
# Get false positives via:
# make -f leishmania.mk assembly ASSEMBLY_TARGET=diff \
#   | grep '^+' | grep '\[[0-9]\{4\}\.' \
#   | python -c 'import re,sys; print "[(" + "), (".join(", ".join(re.findall("(?<=\#)\d+", l.strip())) for l in sys.stdin.readlines()) + ")]"'
fp = []
pos = []
neg = []
with sqlite3.connect(db) as conn:
    c = conn.cursor()
    for S_id, T_id in product(range(1,301), repeat=2):
        if S_id >= T_id:
            continue
        c.execute(q, (S_id, T_id))
        count = int(c.next()[0])
        # mark the false positive data points we want to investigate:
        if (S_id,T_id) in fp or (T_id,S_id) in fp:
            plt.axvline(x=count, color='black', linewidth=0.5, alpha=0.4)
        if count == 0:
            continue
        if (S_id,T_id) in sE or (T_id,S_id) in sE:
            pos += [count]
        else:
            neg += [count]

#plt.axvline(x=min_seeds_for_homology, ymin=0, ymax=1, color='#333333', linewidth=2)
plt.rc('font', family='Cardo')
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
plt.xticks([i*1000 for i in range(int(xmax/1000) + 1)], rotation='vertical')
plt.yticks([i*0.1 for i in range(11)], rotation='vertical')
plt.tick_params(axis='x', labelsize=8, direction='vertical')
plt.xlabel('Number of matching %d-mers' % wordlen)
plt.ylabel('Proportion of read-pairs (cumulative)')
plt.legend(loc='right')
plt.savefig('num_seeds.%d.png' % wordlen)
