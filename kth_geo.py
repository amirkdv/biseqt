#!/usr/bin/env python

from matplotlib import pyplot as plt
from math import ceil, sqrt
from bisect import bisect_right


K = 300 # maximum time
w = 10  # word length
p = 0.85
e = 0.001 # f.n. rate

# u[k][n] is u(n,k)
u = [[0.0] * (w+1) for _ in range(K+1)]

#fig = plt.figure(figsize=(5*w, 5*int(ceil(K/w))))
#ax1 = None
for k in range(0,K+1):
    #ax = fig.add_subplot(int(ceil(K/w))+1, w+1, k+1, sharey=ax1)
    if k == 0:
        u[k][0] = 1.0
        #ax1 = ax
    else:
        for n in range(1, w):
            u[k][n] = p*u[k-1][n-1]
        u[k][w] = p*u[k-1][w-1] + u[k-1][w]
        u[k][0] = 1 - sum(u[k])
    #ax.set_title('k = %d' % k)
    #ax.grid(True)
    #ax.scatter(range(w+1), u[k], s=5)

#plt.savefig('kth-geo.png', dpi=300)

plt.clf()
es = [u[k][w] for k in range(K+1)]
plt.scatter(range(K+1), es, c='k', s=2)
plt.axvline(x=bisect_right(es, 1-e), ymin=0, ymax=1, color='r', label='{:.1e}'.format(e))
plt.grid(True)
plt.legend()
plt.savefig('chain_gap.png', dpi=300)
