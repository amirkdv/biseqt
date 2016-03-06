# Documentation

## Random processes

The following can be generated randomly by the `seq` module:

* Sequence of given length given a letter probabilities.
* Mutants of given sequence given substitution and gap probabilities.
* Random reads from a given sequence given required coverage, read length
  distribution, and substitution/gap probabilities.

## Pairwise alignments

The core dynamic programming algorithm is implemented [in C](https://github.com/amirkdv/oval.py/blob/master/align/liboval.c)
and interfaced using [cffi](https://cffi.readthedocs.org/en/latest/). The
following are supported by `libalign` (the C implementaion) and exposed to
Python by the `pw` module:

* *Global* alignments, i.e Needleman-Wunsch.
* *Local* alignments, i.e Smith-Waterman.
* *Overlap* alignments where only suffix-prefix alignments are considered.
* *Anchored* alignments where only alignments starting/ending at the
  beginning/end of the two sequences are considered,
* *Anchored Overlap* alignments where we demand both the suffix-prefix
  constraint and either of the start or end anchored constraints. This mode of
  alignment is specifically designed for seed expansion in [assembly](#genome-assembly).
* *Affine gap* penalties are supported.
* *Banded* alignments are supported.
* Alignment can be limited to a *frame* for each sequence, i.e using a pair of
  `min_idx,max_idx` provided for both sequences.

The result of an alignment performed by `pw.AlignProblem` is a score and
potentially a traced-back `pw.Transcript` which is represented as a string of
the form:
```
(100,12),12.50:MMMISSSMMMMDDD
```
The first two integers are starting positions of the alignment, the following
number is the score of the alignment and what follows is the `opseq` (edit
transcript) for the alignment.

### Scoring

Gap and substitution probabilities can be translated by `pw.AlignParams` to
corresponding log likelihood scores given the following pieces of information:

* Substitution probabilities,
* Gap probabilities (only a linear gap model is supported for score calculation).
* The random (null hypothesis) distribution of letters.

The calculated score for substitution of letter $a_i$ by letter $a_j$ is given
by:
$$S(a_i,a_j) = \log[1-g] + \log[\Pr(a_j|a_i)] - \log[\Pr(a_j)]$$

Similarly, gap probabilities $g_o$ and $g_e$ are translated to gap scores
where:

* The gap open probability $g_o$ is the probability of a
  single indel following a substitution/match or an indel of a
  different kind.
* The gap extend probability $g_e$ is the probability of a
  single indel following an indel of the same kind.

Note that in the above sense the score (log likelihood) of a gap of
length $n \ge 1$ is:
$$\log g_o + (n-1)\log g_e$$
This differs by the 1 offset from textbook definitions of the
affine gap penalty (and from what `libalign` expects). The two are
equivalent since the above gap penalty function can be rewritten as:
$$\log {g_o \over g_e} + n \log g_e$$

### Complexity

**Time** complexity is $O(mn)$ where $m$ and $n$ are the lengths of the
sequences except for banded global
alignment. In the banded global alignment time complexity is $O(B\max(m,n))$
where $B$ is the band width. As an
example, finding the optimal local alignment for two related sequences (one is
an artificial mutant of the other) with length 5 Kbp takes on
average 3 seconds on a 3.67 GHz quad-core Intel CPU.

**Space** complexity is quadratic too (linear space optimization not
implemented) Currently an average of roughly 100 bytes
is required per cell of the dynamic programming table. For example, it takes
more than 6GB of space to align two sequences of length 8000.

## Tuples

Some tuple methods (aka *k*-mer analysis) are provided by `tuples.TupleDB`
which is backed by [SQLite](https://docs.python.org/2/library/sqlite3.html) to
store, index, and query tuples and by
[Bio.SeqIO](http://biopython.org/wiki/SeqIO) to read FASTA files.
The index provided by `tuples.Index` can be used to index all *k*-mers of a set
of sequences (for given *k*) and to find maximal exactly-matching "seeds".

## Genome assembly

Overlap and layout graphs (i.e OLC minus consensus) can be calculated by methods
provided by `assembly.OverlapBuilder`. All graph algorithms are delegated to
[igraph](http://igraph.org/python/).
Overlap graphs are represented by `assembly.OverlapGraph` (which wraps an
igraph directed graph).
The weighted overlap graph is built as follows:

* For any pair of potentially overlapping reads, find the *shift*
  distribution of all seeds using a rolling sum window. The *shift* of a
  seed with coordinates $(i_S,i_T)$ is the integer $i_S-i_T$.
* Find the ratio of the mode frequency of shifts over the uniform frequency
  (which is 1 over the range of possible shifts). This ratio is taken as
  a measure of "peakedness" of the shift distribution.
    - If the ratio is large enough, the pair of reads are considered
      overlapping with score equal to the overlap length that mode shift
      implies.
    - If the ratio is small enough, the pair of reads are considered
      non-overlapping.
    - If the ratio is neither small or large enough, proceed to seed extension.
* Only considering those seeds with shifts close to the shift mode, try to
  find a seed that extends to a full overlap alignment by consecutive
  start/end-anchored overlap alignments in a moving window along the two
  reads.
    - If any such seed is found, the seeds are considered overlapping with
      score equal to the alignment score of the extended segment.
    - If no such seeds are found, the seeds are considered non-overlapping.

### Cycle breaking

The resulting overlap graph may not be a DAG due to two main reasons:

* wrong weak edges that should not exist.
* strong edges with the wrong direction.
* strong edges that should not exist.

The second case is typically caused by highly overlapping sequences (i.e the
start or end index of end points are too close). Currently such edges are
ignored altogether. The first and third case are delegated to the cycle breaking
algorithm, the latter being the hardest to get rid of.

Regardless, cycle breaking is delegated to `igraph.Graph.feedback_arc_set` which
finds a set of edges the removal of which gives an acyclic graph.
It supports (see [docs](http://igraph.org/python/doc/igraph.GraphBase-class.html#feedback_arc_set))
an optimal, but slow (exponential complexity), integer programming algorithm
(presumably something similar to what is dicussed [here](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.31.5137))
and a suboptimal, but fast, algorithm relying on the [Eades heuristic](http://www.sciencedirect.com/science/article/pii/002001909390079O).
