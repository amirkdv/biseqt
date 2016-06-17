#!/usr/bin/env python

from biseqt.mapping import Mapping
import sqlite3
from matplotlib import pyplot as plt

def load_mappings(path):
    with open(path) as f:
        return eval(f.read())

first = load_mappings('leishmania/blasr.mappings.txt')
second = load_mappings('leishmania/bwa.mappings.txt')
names = {'first': 'Our', 'second': 'BLASR'}
#DB = 'leishmania/genome.mapping.db'
DB = 'leishmania/genome.db'

with sqlite3.connect(DB) as conn:
    lengths = dict(
        [x for x in conn.cursor().execute('SELECT name, LENGTH(seq) FROM seq;')]
    )

first_mapped = set(first.keys())
second_mapped = set(second.keys())
tpl = '%s: from=%s, to=%s, strands=%s\n'

with open('diff.txt', 'w') as f:
    for read in first_mapped - second_mapped:
        f.write('%s: not mapped by second\n' % read)

    for read in second_mapped - first_mapped:
        f.write('%s: not mapped by first\n' % read)

    from_diffs, to_diffs, colors, sizes = [], [], [], []
    for read in first_mapped.intersection(second_mapped).intersection(set(lengths.keys())):
        from_diff = first[read].ref_from - second[read].ref_from
        to_diff = first[read].ref_to - second[read].ref_to
        if first[read].rc == second[read].rc:
            from_diffs += [from_diff]
            to_diffs += [to_diff]
            colors += ['g']
            f.write(tpl % (read, ('%+d' % from_diff).rjust(6), ('%+d' % to_diff).rjust(6), 'same'))
        else:
            colors += ['r']
            f.write(tpl % (read, ' '*6, ' '*6, 'opposite'))
        sizes += [int(lengths[read]/1000)]

# histogram of avg (from/to) diffs for same strand pairs, otherwise won't make
# sense to compare.
n, bins, _ = plt.hist(
    [abs((from_diffs[i] + to_diffs[i])/2) for i in range(len(from_diffs)) if colors[i] == 'g'],
    2000, cumulative=True, normed=True, histtype='step', color='k')
xmax = bins[len(filter(lambda x: n[x]<0.95, range(len(bins)-1)))]
plt.xlim(-xmax/10, xmax)
plt.grid(True)
plt.xlabel('absolute mapping distance b.w. %s and %s' % (names['first'], names['second']))
plt.ylabel('cumulative distribution')
plt.xticks([x*200 for x in range(30)], rotation='vertical')
plt.savefig('diff_hist.png', dpi=300)

plt.clf()
plt.plot([-150000,200000], [-150000,200000], 'k-', lw=2, alpha=0.5)
plt.scatter(from_diffs, to_diffs, marker='o', color=colors, s=sizes, alpha=0.5)
plt.grid(True)
plt.xlabel('starting position difference')
plt.ylabel('ending position difference')
plt.savefig('diffs.png', dpi=300)
