# -*- coding: utf-8 -*-
"""A collection of tools for simulating random processes."""

import numpy as np
import random

from . import ProgressIndicator
from .sequence import Sequence, Alphabet, EditTranscript


def rand_seq(alphabet, size, p=None):
    """Generate a random :class:`Sequence` of the given length from the
    given :class:`Alphabet`.

    Keyword Args:
        alphabet (Alphabet)
        size (int): The length of the randomly generated sequence.
        p (list): The discrete probability distribution for letters of
            alphabet to appear at each position; default is None in which
            case letters are chosen uniformly.
    """
    assert isinstance(alphabet, Alphabet)
    contents = np.random.choice(range(len(alphabet)), size=size, p=p)
    return Sequence(alphabet, contents)


class MutationProcess(object):
    """
    Attributes:
        alphabet (Alphabet): The :class:`sequence.Alphabet` this mutation
            process operates on.
        subst_probs (list): The probability matrix for the
            distribution of substitutions such that ``subst_probs[i][j]`` is
            the probability of the i-th letter being replaced by the j-th
            letter at a non-indel edit operation.  Alternatively, a single
            number can be provided in which case it is treated as the
            probability of *any* substitution to occur and all letters and all
            substitutions are considered equally. For instance, if the single
            number 0.2 is given and the alphabet has 3 letters the probability
            matrix will be:

            .. math::
                \\begin{pmatrix} 0.8 & 0.1 & 0.1 \\\\ 0.1 & 0.8 & 0.1 \\\\
                                 0.1 & 0.1 & 0.8 \\end{pmatrix}

        go_prob (float): probability of a gap starting at any
            position, default is 0 (use 1 for linear gap penalty).
        ge_prob (float): The probability of an open gap to be extended (must
            be at least as large as the gap open probability); default is 0.
        insert_dist (list): the probability distribution for inserted
            content; default is None which is taken to mean uniform.
    """
    def __init__(self, alphabet, subst_probs=None, ge_prob=0, go_prob=0,
                 insert_dist=None):
        assert isinstance(alphabet, Alphabet)
        self.alphabet = alphabet

        if not isinstance(subst_probs, list):
            L = len(self.alphabet)
            assert subst_probs < 1 and subst_probs >= 0
            any_subst = float(subst_probs)
            each_subst = any_subst / (L - 1)
            match = 1 - any_subst
            self.subst_probs = [[match if i == j else each_subst
                                for j in range(L)] for i in range(L)]

        assert go_prob < 1 and ge_prob < 1 and go_prob >= 0 and ge_prob >= 0
        assert go_prob <= ge_prob, 'Gap-open probability cannot be larger ' + \
                                   'than gap-extend probability'
        self.go_prob, self.ge_prob = go_prob, ge_prob
        # let insert_dist remain None; np.random.choice treats it as uniform.
        self.insert_dist = insert_dist

    def mutate(self, seq):
        """Returns a mutant of the given sequence by copying it while at each
        position:

        - either, the current letter is replaced by a random letter according
          to the distribution dictated by ``subst_probs``. This could be a
          match or a substitution.
        - or, a gap is opened (with probability ``go_prob``) the length of
          which follows a geometric distribution with parameter ``ge_prob``.

        Accordingly, a transcript is generated which corresponds to the
        performed edit sequence.

        Args:
            seq (Sequence): The original sequence.
        Returns:
            tuple:
                The mutant :class:`Sequence` and the corresponding
                :class:`pw.EditTranscript`.
        """
        L = len(self.alphabet)
        pos = 0
        T = []
        op, opseq = '', ''

        # np.random.choice defaults to uniform if p == None.
        def rand_let(p): return np.random.choice(L, p=p)

        def coin_toss(p=0.5): return np.random.choice([1, 0], p=[p, 1 - p])

        while pos < len(seq):
            if op and op in 'ID':
                # previous op was an indel, decide whether to extend it:
                if coin_toss(self.ge_prob):
                    if op == 'I':
                        T.append(rand_let(self.insert_dist))
                    else:
                        pos = pos + 1
                else:
                    op = ''  # force the gap to end
            else:
                # previous op is not an indel, decide whether to open a gap:
                if coin_toss(self.go_prob):
                    # It's an insertion or a deletion with equal chance:
                    if coin_toss():
                        op = 'D'
                        pos += 1
                    else:
                        op = 'I'
                        T.append(rand_let(self.insert_dist))
                else:
                    copy = rand_let(self.subst_probs[seq[pos]])
                    T.append(copy)
                    op = 'M' if copy == seq[pos] else 'S'
                    pos += 1

            opseq += op
        return Sequence(self.alphabet, T), EditTranscript(opseq)

    # FIXME
    def rand_read(self, **kw):
        """Generates a random collection of lossy reads from this sequences.
        Each read is generated by taking a substring with Gaussian length and
        uniformly chosen starting point, and mutating it with given
        probabilites.

        Keyword Args:
            subst_probs(List[List[float]]): As in :func:`~mutate`.
            go_prob(Optional[float]): As in :func:`~mutate`.
            ge_prob(Optional[float]): As in :func:`~mutate`.
            coverage (Optional[float]): The expected number of times each
                letter in the sequence appears in the entire read collection.
                Default is 5.
            len_mean (Optional[float]): The mean of the normal distribution of
                read lengths. Default is 100.
            len_sd (Optional[float]): The standard deviation of the normal
                distribution or read lengths. Default is 1.

        Yields:
            tuple: of :class:`Sequence` and ``int`` (starting position).
        """
        coverage = kw.get('coverage', 5)
        len_mean, len_sd = kw.get('len_mean', 100), kw.get('len_sd', 1)
        N = self.length
        num = int(1.0 * N * coverage/len_mean)
        msg = 'generating %d random reads with length ~ ' +  \
            '(mean=%.2f, s.d.=%.2f) for a %d nucl. long genome'
        indicator = ProgressIndicator(
            msg % (num + 2, float(len_mean), float(len_sd), N), num + 2
        )
        indicator.start()

        # include a read that reaches the begenning:
        read, _ = Sequence(self[:len_mean], self.alphabet).mutate(**kw)
        yield read, 0
        indicator.progress()

        for i in range(num):
            # minimum read leangth is 10, and max is N-1
            length = max(10, min(N - 1, int(random.gauss(len_mean, len_sd))))
            start = random.randint(0, N - length)
            read = ''.join([self[k] for k in range(start, start + length)])
            read = Sequence(read, self.alphabet)
            read, _ = read.mutate(**kw)
            yield read, start
            indicator.progress()

        # include a read that reaches the end:
        # FIXME guarantee coverage instead
        read, _ = Sequence(self[-len_mean:], self.alphabet).mutate(**kw)
        yield read, N - len_mean
        indicator.progress()
        indicator.finish()

# TODO get rid of IO using Bio
# TODO figure out pysam; installation requires:
#   $ apt-get install cython zlib1g-dev htslib-dev
#   $ pip install cython pysam
# But then I don't need to do any of the parsing, cf. https://pysam.readthedocs.io/en/latest/api.html

# FIXME this should become a method in SeqDB, or should it?
import sys
import Bio
from uuid import uuid4
def make_sequencing_fixture(genome_file, reads_file, genome_length=1000, **kw):
    """Helper method for tests. Generates a random genome and a random sequence
    of reads from it.

    Args:
        genome_file(str): Path to FASTA file to write genome to.
        reads_file(str): Path to FASTA file to write sequencing reads to.

    Keyword Args:
        genome_length(int): Length of random genome.

    Additionally, all keyword arguments of :func:`randread` are
    allowed and passed as is.
    """
    A = Alphabet('ACGT')
    if 'go_prob' not in kw:
        kw['go_prob'] = 0.1
    if 'go_prob' not in kw:
        kw['ge_prob'] = 0.5
    G = randseq(A, genome_length)
    seqrec = Bio.SeqRecord.SeqRecord(
        Bio.Seq.Seq(str(G)), id='genome', description="(full correct genome)"
    )
    Bio.SeqIO.write([seqrec], genome_file, 'fasta')
    readrecs = []
    for idx, (read, start) in enumerate(G.randread(**kw)):
        seqid = 'R%s_P%d' % (str(uuid4())[:8], start)
        readrecs += [
            Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(str(read)), id=seqid)
        ]
    sys.stderr.write('saving reads to %s.\n' % reads_file)
    Bio.SeqIO.write(readrecs, reads_file, 'fasta')
