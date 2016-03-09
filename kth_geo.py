#!/usr/bin/env python

from matplotlib import pyplot as plt
from math import ceil, sqrt
from bisect import bisect_right


extra_time = 20
w = 10  # word length
p = 0.85
e = 0.001 # f.n. rate

# u[k][n] is u(n,k)
u = [[0.0] * (w+1)]
u[0][0] = 1.0
from pprint import pprint as pp

extra = 0
while u[-1][w] < 1-e:
    u.append([0.0] * (w+1))
    for n in range(1, w):
        u[-1][n] = p*u[-2][n-1]
    u[-1][w] = p*u[-2][w-1] + u[-2][w]
    u[-1][0] = 1 - sum(u[-1])

plt.clf()
plt.rc('text', usetex=1)
plt.scatter(range(len(u)), [u[k][w] for k in range(len(u))], c='k', s=2)
plt.axvline(x=len(u), ymin=0, ymax=1, color='r', label='$u(n,w)\ge 1-({:.1e})$'.format(e))
plt.xlabel('sub-alignment length')
plt.ylabel('Pr[seed observed]')
plt.xticks(list(plt.xticks()[0]) + [len(u)], rotation=90, fontsize=8)
plt.grid(True)
plt.legend()
plt.savefig('chain_gap.png', dpi=300)
