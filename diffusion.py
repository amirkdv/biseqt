#!/usr/bin/env python

from matplotlib import pyplot as plt
from scipy.special import erfinv
from math import sqrt, pi, exp, cos, sin, erf, ceil

K = 10000
p = 0.1
p0 = 0.001
w = int(ceil(2*sqrt(p*K)*erfinv(1-p0)))
print 'max ops: %d, band radius: %d' % (K, w)

def radii():
    radius = lambda x: int(ceil(2*sqrt(p*x)*erfinv(1-p0)))
    xs = [x*100 for x in range(100)]
    ys = [radius(x) for x in xs]
    return xs, ys

#plt.plot(xs, ys, '-k')
#plt.grid(True)
#plt.savefig('radius.png', dpi=300)

plt.rc('text', usetex=1)
ticks_every = 20
nrange = range(-2*w, 2*w +1)
plt.xticks([x * ticks_every + ticks_every*(nrange[0]/ticks_every) for x in range(len(nrange)/ticks_every+1)], rotation='vertical', fontsize=10)
plt.yticks([0.1*n for n in range(11)], fontsize=10)
plt.grid(True)
plt.ylim(-0.1, 1.1)
plt.xlabel('$n$')
plt.ylabel('$u(n,K)$')
plt.title('band randius = %d for $p_\circ=%.4f$ after $K=%d$ steps' % (w, p0, K), size=10)


def recurrence_solution():
    u = {}
    for k in range(K+1):
        if k == 0:
            u[k] = {n: 1 if abs(n) < w + K else 0 for n in nrange}
            continue
        u[k] = {}
        for n in nrange:
            if abs(n) >= w:
                u[k][n] = 0
            else:
                u[k][n] = p*u[k-1][n-1] + p*u[k-1][n+1] + (1-2*p)*u[k-1][n]
    return u

u = recurrence_solution()
plt.scatter(nrange, [u[K][n] for n in nrange], edgecolor='g', facecolor='none', lw=0.5, s=20, label='recurrence solution')

def fs_solution(num_terms):
    cos_rl = lambda n: ((n+0.5)*pi)/w
    sin_rl = lambda n: n*pi/2
    fourier_cos = lambda x, t, n: (2*(-1)**n / ((n+0.5)*pi)) * exp(t*-p*cos_rl(n)**2) * cos(cos_rl(n)*x)
    fourier_sin = lambda x, t, n: (2*(-1)**(n+1) / (n*pi))   * exp(t*-p*sin_rl(n)**2) * sin(sin_rl(n)*x) if n else 0
    diffusion_correct = lambda N, x, t: sum(fourier_cos(x,t,n) + fourier_sin(x,t,n) for n in range(N))
    return [diffusion_correct(num_terms, n, K) if abs(n) <= w else 0 for n in nrange]

num_terms = 4
soln = fs_solution(num_terms)
plt.plot(nrange, soln, color='k', label='%d terms of F.S.: $u(0, K)=%.4f$' % (num_terms*2-1, soln[nrange.index(0)]))

def nobc_approx_solution():
    return [0.5 * ( erf((n+w)/(2*sqrt(p*K))) - erf((n-w)/(2*sqrt(p*K))) ) for n in nrange]

soln = nobc_approx_solution()
plt.plot(nrange, soln, color='r', label='ignoring B.C: $\hat{u}(0, K)=%.4f$' % soln[nrange.index(0)])
plt.legend(fontsize=8)
plt.savefig('diffusion.png', dpi=300)
