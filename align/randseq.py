#!/usr/bin/env python
import numpy as np
import random

# distribution is a dict of letter -> probability
def gen(length, dist={k:0.25 for k in 'ACGT'}):
    """Generates a random sequence of the specified length within the alphabet
    specified by the keys in the distribution matrix. The distribution matrix
    should look like this for nucleotide sequences:
        {'A': 0.25,
         'C': 0.25,
         'G': 0.25,
         'T': 0.25}

    :param length(int): length of generated sequence.
    :param dist(dict): keys are the letters in the alphabet and values are the
        probability of an arbitrary nucleatoride being each letter.
    """
    assert abs(1-sum(dist.values())) < 0.001
    space = []
    for k in dist:
        space.extend([k] * int(100*dist[k]))
    return ''.join([random.choice(space) for _ in range(length)])

def mutate(S, gap_open=None, gap_continue=0.5, rates=None):
    """Mutates a given sequence with specified probabilities. The sequence is
    scanned and copied to the mutated sequence where at each position:
    * the current letter will be replaced by an arbitrary letter with a
        distribution that depends on the original letter.
    * with a certain fixed probability a gap may be openned with random length
        with geometric distribution.
    Accordingly, an opseq (see `align.solve()`) is generated which corresponds
    to the performed edit sequence.

    :param S(str): original sequence.
    :param gap_open(float): probability of a gap starting at any position.
    :param gap_continue(float): Bernoulli success probability of the gap
        extension distribution
    :param rates(dict): a letter-by-letter error probability matrix. For
        nucleotides, for example, it should have the following structure:

            {'A':{'A': 0.7,
                  'C': 0.1,
                  'G': 0.1,
                  'T': 0.1},
             'C':{...
             ...
            }
    """
    T = ''
    k = 0
    transcript = ''
    while k < len(S):
        if gap_open is not None:
            assert(gap_open < 1) # if not none, gap_open is assumed to be the gap probability
            assert(gap_continue < 1)
            if random.randint(0, 1000) < gap_open * 1000:
                length = np.random.geometric(1 - gap_continue)
                if random.choice([0,1]):
                    # deletion
                    transcript += 'D' * length
                    k += length
                    continue
                else:
                    # insertion
                    transcript += 'I' * length
                    T += gen(length)
                    k += 1
                    continue
        T += gen(1, dist=rates[S[k]])[0]
        transcript += 'M' if T[-1] == S[k] else 'S'
        k += 1
    return (T, transcript)
