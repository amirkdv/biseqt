import sqlite3
from itertools import product
from matplotlib import pyplot as plt
import igraph

wordlen = 15
num_bins = 50
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
# Get this via:
# make -f leishmania.mk assembly ASSEMBLY_TARGET=diff | grep '^+' | grep '\[[0-9]\{4\}\.' | python -c 'import re,sys; print "), (".join(", ".join(re.findall("(?<=\#)\d+", l.strip())) for l in sys.stdin.readlines())'
fp = [(145, 177), (9, 243), (55, 145), (55, 125), (49, 57), (104, 223), (15, 23), (243, 223), (263, 227), (263, 93), (192, 217), (258, 32), (223, 259), (191, 259), (191, 104), (191, 243), (191, 36), (191, 85)]
pos = []
neg = []
with sqlite3.connect(db) as conn:
    c = conn.cursor()
    for S_id, T_id in product(range(1,301), repeat=2):
        if S_id >= T_id:
            continue
        c.execute(q, (S_id, T_id))
        count = int(c.next()[0])
        if (S_id,T_id) in fp or (T_id,S_id) in fp:
            plt.axvline(x=count, color='black', linewidth=0.5, alpha=0.4)
        if count == 0:
            continue
        if (S_id,T_id) in sE or (T_id,S_id) in sE:
            pos += [count]
        else:
            neg += [count]

plt.hist(neg, num_bins, histtype='stepfilled', color='red', alpha=0.3, cumulative=True)
plt.hist(pos, num_bins, histtype='stepfilled', color='green', alpha=0.3, cumulative=True)
plt.grid(True)
plt.xlabel('Seed count')
plt.ylabel('Number of read pairs')
plt.savefig('num_seeds.png',dpi=200)
