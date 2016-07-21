# -*- coding: utf-8 -*-
"""
.. wikisection:: overview
    :title: Random Processes

    The :mod:`biseqt.stochastics` module provides a collection of tools for
    statistical calculations and for simulating random processes.

    >>> from biseqt.sequence import Alphabet
    >>> from biseqt.stochastics import rand_seq, MutationProcess
    >>> A = Alphabet('ACGT')
    >>> seq = rand_seq(A, 10)
    >>> seq
    Sequence(Alphabet(["A","C","G","T"]), \
    contents=(1, 3, 3, 2, 2, 0, 0, 0, 0, 0))
    >>> print(seq)
    CTTGGAAAAA

    >>> P = MutationProcess(A, go_prob=.1, ge_prob=.2, subst_probs=.3)
    >>> mutant, transcript = P.mutate(seq)
    >>> transcript
    EditTranscript("SMIMMMMSSMS")

    >>> from biseqt.io import pw_render_term
    >>> print(pw_render_term(transcript, seq, mutant, colored=False))
    origin[0]: CT-TGGAAAAA
    mutant[0]: ATCTGGATCAT
"""
from math import sqrt, erfc, log
from itertools import product
from scipy.special import erfcinv
import numpy as np

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


def rand_read(seq, len_mean=None, len_sd=1, expected_coverage=None, num=None):
    """Generates a random collection of lossless reads from the given sequence.
    Each read is generated by taking a substring with Gaussian length and
    uniformly chosen starting point. The number of reads is either determined
    as an argument or calculated based on the expected coverage:

    .. math::
        num = \\frac{\\mathbb{E}[coverage]\\times |seq|}{\\mathbb{E}[len]}

    Args:
        seq (sequence.Sequence): The original sequence to be sampled.

    Keyword Args:
        num (int): Number of reads to generate; default is None in which case
            ``expected_coverage`` takes precedence. If both are None a single
            read is returned.
        expected_coverage (float): The expected number of times each letter in
            the sequence ought to appear in the entire read collection; default
            is None in which case ``num`` takes precedence. If both are None
            a single read is returned.
        len_mean (float): The mean of the normal distribution of
            read lengths. Default is 100.
        len_sd (float): The standard deviation of the normal
            distribution or read lengths. Default is 1.

    Yields:
        tuple:
            The read :class:`Sequence` and its starting position ``int``.
    """
    assert len_mean < len(seq), \
        'Expected read length must be smaller than the sequence length'
    assert num is None or expected_coverage is None, \
        'At most one of expected_coverage or num can be specified'
    if num is None:
        if expected_coverage is None:
            num = 1
        else:
            num = int(1. * len(seq) * expected_coverage / len_mean)

    for length in np.random.normal(loc=len_mean, scale=len_sd, size=num):
        # force a minimum read length of 1 and max |S| - 1
        length = max(1, min(len(seq) - 1, int(length)))
        start = np.random.choice(len(seq) - length)
        read = seq[start:start + length]
        yield read, start


def binomial_to_normal(n, p):
    """Given the parameters of a binomial distribution, returns the parameters
    of its normal approximation.

    .. math::
        B(n, p) \\simeq \\mathcal{N}(np, \\sqrt{np(1-p)})

    The approximation approaches identity as :math:`n` grows to infinity for
    any fixed :math:`p`.

    Args:
        n (int): First parameter of binomial distribution (i.e number of
            Bernoulli trials).
        p (float): Second parameter of binomial distribution (i.e Bernoulli
            success probability).

    Returns:
        tuple:
            Mean and standard deviation of the approximating normal
            distribution.
    """
    assert p >= 0 and p <= 1 and n > 0
    mu = n * p
    sd = sqrt(n * p * (1 - p))
    return mu, sd


def normal_neg_log_pvalue(mu, sd, x):
    """Gives the negative log p-value of an observation under the null
    hypothesis of normal distribution; that is, given an observation :math:`x`
    from a random variable:

    .. math::
        X \\sim \\mathcal{N}(\\mu, \\sigma)

    this function calculates :math:`-\\log(\\Pr[X\\ge x])` which is:

    .. math::
        -\log\left(
            \\frac{1}{2} \left[1 - \mathrm{erf} \left(
                                \\frac{x-\mu}{\sigma\sqrt{2}} \\right)
                         \\right]
        \\right)

    Args:
        mu (float): Mean of normal distribution.
        sd (float): Standard deviation of normal distribution.
        x (float): Observation.

    Returns:
        float:
            A positive real number or infinity.

    .. wikisection:: dev
        :title: Log-probability numerics

        It is customary to capture the statistical significance of an
        observation :math:`x` corresponding to the random variable :math:`X`
        by a *score*

        .. math::
            S=-\log(p)

        where :math:`p` is the p-value of the observation, i.e :math:`\Pr(X\ge
        x)` as per the null hypothesis.  If the null hypothesis is :math:`X
        \sim \mathcal{N}(\mu, \sigma)` the score is given in closed form by:

        .. math::
            S = -\log\\left(
                    \\frac{1}{2} \left[1 - \mathrm{erf} \left(
                                    \\frac{x-\mu}{\sigma\sqrt{2}} \\right)
                                 \\right]
                \\right)
              = -\log\left(1 - \mathrm{erf} \left(
                                    \\frac{z}{\sqrt{2}} \\right)
                \\right) + \log(2)

        Calculating this formula numerically runs into precision issues because
        the error function rapidly approaches its limit such that, for
        instance, on 64-bit system, the :math:`\mathrm{erf}` term evaluates to
        ``0.0`` for any z-score larger than 9 leading to infinities where the
        real value of :math:`S` is a small number, e.g. at :math:`z=9` we
        have :math:`S\simeq42.9`.

        One improvement is to use the ``erfc`` `function <erfc_>`_
        (or corresponding wrappers in ``numpy`` or ``scipy.special``) which
        numerically finds :math:`1-\mathrm{erf}(\cdot)`. This postpones
        infinities to z-scores larger than 39 at which point the score jumps
        from roughly 726 to infinity.

        These blowups can have serious consequences in any classification
        algorithm where the cumulative distribution of a population is of
        interest: with too many infinities in a set of scores, the cumulative
        distribution of scores never gets close to 1 (e.g. if a fifth of scores
        are infinities the maximum value of the cumulative distribution is
        0.8). For this reason, it is advisable to use z-scores instead of
        negative log of p-values for classification (note that the ROC curve
        corresponding to the two is identical, cf. below).


        .. _erfc: https://docs.python.org/2/library/math.html#math.erfc

    .. wikisection:: dev
        :title: Effect of filtering on ROC curves
        :parent: Log-probability numerics

        Consider a binary classification problem with two sample sets
        :math:`X,Y` with positive and negative labels respectively. The
        corresponding ROC curve is a paremetric curve in the unit square given
        by:

        .. math::
            \left(F_X(t), F_Y(t)\\right) \in [0,1]\\times[0,1]

        where :math:`F_X,F_Y` are the cumulative distributions of the positive
        and negative samples and :math:`t` varies over the extended real line
        :math:`\mathbb{R}\cup\{\pm\infty\}`.

        Now let :math:`f:\mathbb{R}\\to\mathbb{R}` be a *monotonically
        increasing* real function and for any set :math:`A` let :math:`f(A)`
        denote the set :math:`\{f(a); a\in A\}`. If we consider the sets
        :math:`f(X), f(Y)` as the positive and negative samples, we have:

        .. math::
            F_{f(X)}(a) = \Pr[f(X)\le a] = \Pr[X\le f^{-1}(a)]

        which implies:

        .. math::
            F_{f(X)} = F_X \circ f^{-1}

            F_{f(Y)} = F_Y \circ f^{-1}

        It follows that although the cumulative distributions :math:`F_{f(X)},
        F_{f(Y)}` differ from the original :math:`F_X,F_Y`, the ROC curve
        is merely reparameterized and thus remains unchanged.

        As a special case, this implies that the ROC curve of a collection of
        positive and negative z-score sets is identical to the ROC curve of the
        corresponding sets of p-values or negative-log p-values.
    """
    z_score = (x - mu) / float(sd)
    try:
        return - log(erfc(z_score/sqrt(2)) / 2.)
    except ValueError:
        # NOTE This can only happen if the argument to log is
        # 0, which theoretically means z_score >> 1. In practice, this happens
        # for z scores exceeding 39.
        return float('+inf')


def band_radius(len0, len1, diag, gap_prob=None, sensitivity=None):
    """Calculates the smallest band radius in the dynamic programming table
    such that an overlap alignment, with the given gap probability, stays
    entirely within the diagonal band centered at the given diagonal. This is
    given by:

    .. math::
        r = 2\\sqrt{g(1-g)K}
            \\mathrm{erf}^{-1}\\left(1-\\frac{2\\epsilon}{3}\\right)

    where :math:`g` is the gap probability, :math:`1-\\epsilon` is the desired
    sensitivity, and :math:`K` is the "expected" length of the alignment given
    by:

    .. math::
        K = \\left(\\frac{2}{2 - g}\\right) L

    where :math:`L` is the maximum possible length of the alignment:

    .. math::
        L = \\min(l_0 - d, l_1) + \\min(d, 0)

    with :math:`l_0,l_1` being the length of the sequences (i.e ``len0`` and
    ``len1`` arguments) and :math:`d` the starting diagonal (i.e ``diag``
    argument).

    Diagonals are numbered as follows in the dynamic programming table::

        0 -1 -2 -3 ... -len1
        1
        2
        3
        .
        .
        .
        +len0


    Args:
        len0 (int): Length of the first sequence (the "vertical" sequence in
            the table).
        len1 (int): Length of the second sequence (the "horizontal" sequence in
            the table).
        diag (int): Starting diagonal of alignments to consider.
        gap_prob (float): Probability of indels occuring at any position of an
            alignment.
        sensitivity (float): The probability that an alignment with given gap
            probability remains entirely within the band.
    Returns:
        int: The smallest band radius guaranteeing the required sensitivity.

    """
    assert sensitivity > 0 and sensitivity < 1
    assert gap_prob > 0 and gap_prob < 1

    epsilon = 1. - sensitivity
    adjusted_epsilon = epsilon * 2 / 3

    max_len = min(len0 - diag, len1) + min(diag, 0)
    expected_len = (2. / (2 - gap_prob)) * max_len
    assert expected_len >= 0
    radius = 2 * erfcinv(adjusted_epsilon) * sqrt(
        gap_prob * (1 - gap_prob) * expected_len
    )
    return max(1, int(radius))


def band_radius_calculator(gap_prob=None, sensitivity=None):
    """Creates a function that behaves like :func:`band_radius` for a fixed
    set of parameters (intended for bulk usage).

    Keyword Args:
        gap_prob (float): as in :func:`band_radius`.
        sensitivity (float): as in :func:`band_radius`.

    Returns:
        function: A function with signature ``f(len0, len1, diag)``.
    """
    assert sensitivity > 0 and sensitivity < 1
    assert gap_prob > 0 and gap_prob < 1

    epsilon = 1. - sensitivity
    adjusted_epsilon = epsilon * 2 / 3

    C = 2 * erfcinv(adjusted_epsilon) * sqrt(gap_prob * (1 - gap_prob))

    def calculator(len0, len1, diag):
        max_len = min(len0 - diag, len1) + min(diag, 0)
        expected_len = (2. / (2 - gap_prob)) * max_len
        assert expected_len >= 0
        radius = C * sqrt(expected_len)
        return max(1, int(radius))

    return calculator


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

        go_prob (float): probability of a single indel following a subtitution,
            a match, or an indel of a different kind (i.e from insertion to
            deletion); default is 0 (cf. :attr:`ge_prob`).
        ge_prob (float): The probability of an open gap to be extended by an
            indel of the same kind; default is 0. For consistency, it is always
            required that the gap-extend probability to be at least as large as
            the gap-open probability with equality implying a linear model
            which translates to a linear gap penalty scheme in
            :func:`log_odds_scores`.
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
          to the distribution dictated by :attr:`subst_probs`. This could be a
          match or a substitution.
        - or, a gap is opened (with probability :attr:`go_prob`) the length of
          which follows a geometric distribution with parameter
          :attr:`ge_prob`.

        Accordingly, a transcript is generated which corresponds to the
        performed edit sequence.

        Args:
            seq (Sequence): The original sequence.
        Returns:
            tuple:
                The mutant :class:`Sequence <biseqt.sequence.Sequence>` and the
                corresponding :class:`EditTranscript
                <biseqt.sequence.EditTranscript>`.
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

    def noisy_read(self, seq, **kw):
        """Wraps :func:`rand_read` to generates a collection of lossy reads
        from the given sequence. First, lossless reads are generated by
        :func:`rand_read` (all keyword arguments are passed as is), then each
        read is modified by this mutation process, as in :func:`mutate`.

        Args:
            seq (sequence.Sequence): The original sequence to be sampled.

        Keyword Args:
            num (int): Number of reads to generate; if not given
                ``expected_coverage`` must be specified.
            expected_coverage (float): The expected number of times each letter
                in the sequence ought to appear in the entire read collection;
                if not given ``num`` must be specified.
            len_mean (float): The mean of the normal distribution of read
                lengths.
            len_sd (float): The standard deviation of the normal
                distribution or read lengths; default is 1.

        Yields:
            tuple:
                The noisy read :class:`Sequence <biseqt.sequence.Sequence>`,
                the starting position ``int`` of the original lossless read,
                and the corresponding :class:`EditTranscript
                <biseqt.sequence.EditTranscript>`.
        """
        for read, start in rand_read(seq, **kw):
            read, tx = self.mutate(read)
            yield read, start, tx

    def log_odds_scores(self, null_hypothesis=None):
        """Converts the mutation probabilities (substitution and gap) to
        log-odds scores usable for sequence alignment (cf.
        :class:`AlignScores`).

        Keyword Args:
            null_hypothesis(list): The probability distribution of letters to
                use as the null hypothesis; default is None which is taken to
                mean uniform distribution.

        Returns:
            tuple:
                The substitution score matrix (a list of length equal to the
                alphabet containing lists each of the same length containing
                the substitution scores) and the gap scores as a tuple of
                gap-open and gap-extend scores.

        .. rubric:: Substitution scores

        The scores are natural logs of odds ratios, namely the substitution
        score of letter :math:`a_i` to letter :math:`a_j` is:

        .. math::

            S(a_i\\rightarrow a_j) = \log[(1-g)\Pr(a_j|a_i)] - \log[\Pr(a_j)]

        where :math:`g` is the gap-extend
        probability (see below), :math:`\Pr(a_j|a_i)` is the :attr:`substution
        probability <subst_probs>`, and :math:`\Pr(a_j)` is the null hypothesis
        distribution.

        .. rubric:: Gap scores

        Gap :attr:`open <go_prob>` and :attr:`extend <ge_prob>` probabilities
        are converted to an affine gap penalty scheme (but reported values are
        "scores" not penalties, i.e they are negative). The probabilitstic
        model is as described in :func:`MutationProcess` (also cf.
        :func:`mutate`).  Accordingly, given the gap-open and gap-extend
        probabilities :math:`g_o, g_e`, the score of a gap of length :math:`n
        \\ge 1` is

        .. math::
            \\log g_o + (n-1)\\log g_e =
            \\log\\left( \\frac{g_o}{g_e} \\right) + n\\log g_e

        where the first term in the RHS is the reported gap-score and the
        coefficient of :math:`n` in the second term is the reported gap-extend
        score.

        Note:
            * In the calculation of substitution scores the gap-open
              probability is ignored and a linear model, with its sole
              parameter equal to the gap-extend probability, is assumed. This
              is because under an affine model where the probability of opening
              a gap differs from that of extending a gap, the substitution
              probabilities also become context-dependent (see where :math:`g`
              appears in the formula above) which is not supported by
              ``pwlib``.
            * The above formula for gap scores is the technical reason why we
              always demand that :math:`g_e\\ge g_o` so that the gap score is
              not positive, with equality implying linear gap penalties (i.e
              the gap-open *score* is 0 when the gap-open and gap-extend
              *probabilities* are identical).
        """
        if null_hypothesis is None:
            null_hypothesis = [1./len(self.alphabet)] * len(self.alphabet)
        err = 'Zero probabilities are not allowed for score calculation'
        assert all(x > 0 and x < 1 for x in null_hypothesis), err
        assert all(x > 0 and x < 1 for y in self.subst_probs for x in y), err
        assert self.ge_prob > 0 and self.go_prob > 0, err

        subst_scores = [[0] * len(self.alphabet) for _ in self.alphabet]
        for i, j in product(range(len(self.alphabet)), repeat=2):
            subst_scores[i][j] = log(1 - self.ge_prob) + \
                log(self.subst_probs[i][j]) - log(null_hypothesis[j])
        gap_scores = log(self.go_prob) - log(self.ge_prob), log(self.ge_prob)
        return subst_scores, gap_scores
