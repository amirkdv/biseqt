#!/usr/bin/env python

from biseqt.mapping import Mapping
import sqlite3
from matplotlib import pyplot as plt

def load_mappings(path):
    with open(path) as f:
        return eval(f.read())

bwa = load_mappings('leishmania/bwa.mappings.txt')
blasr = load_mappings('leishmania/blasr.mappings.txt')

with sqlite3.connect('genome.leishmania.db') as conn:
    lengths = dict(
        [x for x in conn.cursor().execute('SELECT name, LENGTH(seq) FROM seq;')]
    )

title = lambda read: '%s (len=%s)' % (read, str(lengths[read]).rjust(5))

bwa_mapped = set(bwa.keys())
blasr_mapped = set(blasr.keys())

with open('diff.txt', 'w') as f:
    for read in bwa_mapped - blasr_mapped:
        f.write('%s: not mapped by blasr\n' % title(read))

    for read in blasr_mapped - bwa_mapped:
        f.write('%s: not mapped by bwa\n' % title(read))

    from_diffs, to_diffs, colors, sizes = [], [], [], []
    for read in bwa_mapped.intersection(blasr_mapped):
        from_diff = bwa[read].ref_from - blasr[read].ref_from
        to_diff = bwa[read].ref_to - blasr[read].ref_to
        if bwa[read].strand == blasr[read].strand:
            from_diffs += [from_diff]
            to_diffs += [to_diff]
        strand_diff = bwa[read].strand != blasr[read].strand
        colors += ['r' if strand_diff else 'g']
        sizes += [int(lengths[read]/1000)]
        tpl = '%s: from=%s, to=%s, strands=%s\n'
        if strand_diff:
            f.write(tpl % (title(read), ' '*6, ' '*6, 'opposite'))
        else:
            f.write(tpl % (title(read), ('%+d' % from_diff).rjust(6), ('%+d' % to_diff).rjust(6), 'same'))

# histogram of avg (from/to) diffs for same strand pairs, otherwise won't make
# sense to compare.
n, bins, _ = plt.hist(
    [abs((from_diffs[i] + to_diffs[i])/2) for i in range(len(from_diffs)) if colors[i] == 'g'],
    2000, cumulative=True, normed=True, histtype='step', color='k')
xmax = bins[len(filter(lambda x: n[x]<0.95, range(len(bins)-1)))]
plt.xlim(-xmax/10, xmax)
plt.grid(True)
plt.xlabel('absolute mapping distance b.w. BWA and BLASR')
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
