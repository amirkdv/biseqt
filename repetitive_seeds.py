from matplotlib import pyplot as plt
import sqlite3
import os
from math import erf, sqrt, log

normal_pvalue = lambda mu, sd, x: 0.5 * (1 - erf((x - mu) / (sd * sqrt(2))))
# Normal approximation (mean and standard deviation) of a binomial distribution
B2N = lambda n, p: (n*p, sqrt(n * p * (1-p)))
prob_word = lambda length: 0.25 ** length

WORDLEN = int(os.environ['WORDLEN'])
pvalues = {}
with sqlite3.connect('genome.leishmania.hp_assembly.db') as conn:
    c = conn.cursor()
    c.execute('select sum(length(seq)) from seq')
    L = int(c.next()[0]) # total sequence length
    mu, sd = B2N(L, prob_word(WORDLEN)) # approximating normal distribution of word counts

    c.execute('select count(*) from tuples_%d' % WORDLEN)
    N = int(c.next()[0]) # total number of words


    c.execute('select tuple, hits from tuples_%d' % WORDLEN)
    for row in c:
        # apply a bonferroni correction since we are testing N hypotheses
        pvalue = N * normal_pvalue(mu, sd, row[1].count('@'))
        if pvalue:
            pvalues[row[0]] = log(pvalue)
        #else:
            #print row[0]

plt.grid(True)
plt.hist(pvalues.values(), 100, normed=True, histtype='step', cumulative=True, color='k')
plt.savefig('repetitive_seeds.png', dpi=300)
