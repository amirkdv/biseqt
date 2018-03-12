#!/usr/bin/env python
from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess
from biseqt.blot import HomologyFinder
from util import savefig

import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_scored_seeds_3d(fig, ax, scored_seeds):
    idx_S, idx_T1, idx_T2, cs, ss = [], [], [], [], []
    cmap = plt.cm.get_cmap('jet')
    scores = [score for _, score in scored_seeds]
    max_score = max(scores)
    for (i, j, k), score in scored_seeds:
        idx_S.append(i)
        idx_T1.append(j)
        idx_T2.append(k)
        cs.append(cmap(score/max_score)[:3])
        ss.append(10 if score > 100 else 1)
        # ss.append(10 if score > -1 else 1)

    ax.scatter(idx_S, idx_T1, idx_T2, facecolor=cs, lw=0, s=ss, alpha=.3)
    ax.set_aspect('equal')
    ax.elev = 10
    ax.azim = 135
    norm = matplotlib.colors.Normalize(vmin=0, vmax=max_score)
    m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    m.set_array(scores)
    fig.colorbar(m)


subst, gap = .15, .05
wordlen = 6
K = 1000
path = ':memory:'
A = Alphabet('ACGT')
M = MutationProcess(A, subst_probs=subst, ge_prob=gap, go_prob=gap)
S = rand_seq(A, K)
T1, _ = M.mutate(S)
T2, _ = M.mutate(S)


def junk(): return rand_seq(A, 2*K)


S, T1, T2 = junk() + S + junk(), junk() + T1 + junk(), junk() + T2 + junk()
HF_kw = {'g_max': .2, 'sensitivity': .9, 'alphabet': A, 'wordlen': wordlen,
         'path': ':memory:'}
HF = HomologyFinder(S, T1, T2, **HF_kw)

scored_seeds = HF.score_seeds(K_min=100, p_min=(1-gap) * (1-subst))
# convert to ij coordinates and leave only the H0 score
scored_seeds = [(HF.to_ij_coordinates(ds, a), s0)
                for (ds, a), (s0, s1) in scored_seeds]

fig = plt.figure()
ax = fig.gca(projection=Axes3D.name)
plot_scored_seeds_3d(fig, ax, scored_seeds)

for axis in 'xyz':
    ax.tick_params(axis=axis, labelsize=5)
ax.set_xlabel('Sequence 1')
ax.set_ylabel('Sequence 2')
ax.set_zlabel('Sequence 3')

ax.set_title('Multiple Sequence Comparison (H0 score)')

fig.tight_layout()
savefig(fig, 'multiple-sequence.png', dpi=300)
