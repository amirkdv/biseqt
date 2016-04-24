#!/usr/env/bin python

xdata = []
ydata = []
# to generate results:
# for score in -1 -2 -3; do for i in $(seq 50); do SCORE=$score python -m biseqt.tests.pw
#   | grep '^(' | head -c 20 | cut -d\( -f2 | cut -d, -f1
#   | xargs echo $score, ; done ; done | tee results.txt
with open('mismatch_score_data.txt') as f:
    data = eval('[(' + f.read().strip().replace('\n', '), (') + ')]')

xdata = [i[0] for i in data]
ydata = [i[1] for i in data]

from matplotlib import pyplot as plt

plt.scatter(xdata, ydata, s=2, color='k')
plt.rc('text', usetex=1)
plt.xlabel('Mismatch score')
plt.ylabel('Optimal shift')
plt.axvline(x=-1, lw=3, alpha=0.4, color='b', label='Correct mismatch score\n(match score is fixed at +1)')
plt.axhline(y=3000, lw=3, alpha=0.4, color='g', label='Correct shift\n(sequence lengths are 4 and 10 kbp)')
plt.legend(fontsize=9)
plt.grid(True)
plt.savefig('mismatch_score_sensitivity.png', dpi=300)
