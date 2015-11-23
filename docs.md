# Documentation

## Random processes

The following can be generated randomly by the `seq` module:

* Sequence of given length given a letter probabilities.
* Mutants of given sequence given substitution and gap probabilities.
* Random reads from a given sequence given required coverage, read length
  distribution, and substitution/gap probabilities.

## Pairwise alignments

The core dynamic programming algorithm is implemented [in C](https://github.com/amirkdv/align.py/blob/master/align/libalign.c)
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

## Alphabet translation

To deal the high rates of homopolymeric errors in certain sequencing methods,
`homopolymeric.HpCondenser` can be used to translate (aka "condense") sequences
into a new alphabet where each homopolymeric substring is translated to
a single letter(e.g `AACCC` becomes `A2C3`).

### Alignment in condensed alphabet

The machinery in `libalign` is capable of aligning sequences in alphabets with
letters longer than a single character. The only requirement is that all letters
in an alphabet have the same length (to avoid a book-keeping mess). This leads
to a [caveat](#a-caveat) discussed below. Aside from the caveat, the following operations
are supported:

* "Condensing" transcripts into transcripts for condensed sequences and
  "expanding" transcripts for condensed sequences back into one for the original
  sequences.
* "Condensing" alignment parameters into corresponding parameters in the
  condensed alphabet,
* "Condensing" a seed (see the section on [assembly](#genome-assembly) below)
  into a seed in the condensed alphabet.

#### Score translation

To get a set of parameters for alignment in the condensed alphabet there are
two options:

* Convert *scores* in the source alphabet directly into scores in the condensed
  alphabet. Additionally scores for homopolymeric indels must be provided.
* Convert *probability parameters* for substitution and indels into
  corresponding probabilities in the condensed alphabet and converting those
  into scores as usual. Additionally homopolymeric indel probabilities must
  be provided.
  The translation formula is the following when the length of letters are identical:
  $$\Pr(x_i \rightarrow y_i) = \Pr(x \rightarrow y)^i(1-g_h)^{i-1}$$
  where $g_h$ is the homopolymeric gap probability (only a linear model is
  supported). When the length of letters differ, say $i<j$, we have the
  following where $\pi(\cdot)$ is the integer partition function:
  $$\Pr(x_i \rightarrow y_j) = \pi(i)\Pr(x \rightarrow y)^i(1-g_h)^{i-1}g_h^{|i-j|}$$


*Note*: These calculations here may have serious errors. In fact, the
calculated probabilities as described above don't necessarily add
up to 1! Returned probability matrix is normalized in each row
to make sure the output is not terribly wrong.

### A Caveat

Condensing a sequence requires specifying a whole number `maxlen`: homopolymeric
substrings longer than `maxlen` are considered to have only `maxlen` characters.
This is needed to ensure constant letter length across the condensed alphabet
(which is required by `libalign`).

Due to this, if source alphabet sequences contain homopolymeric substrings that
are longer than the specified `maxlen`, the condensing process is lossy
(expanding a condensed sequence does not necessarily give its original
sequence). However, if the original sequence is available, expanding an
alignment transcript can be done losslessly to match the original sequence.

## Genome assembly

Overlap and layout graphs (i.e OLC minus consensus) can be calculated by methods
provided by `assembly.OverlapBuilder`. All graph algorithms are delegated to
[igraph](http://igraph.org/python/).
Overlap graphs are represented by `assembly.OverlapGraph` (which wraps an
igraph directed graph).
The weighted overlap DAG is built by seed expansion (see
[Tuples Methods](#tuples)) on all pairs of sequences and the longest
path is reported as the layout. Expansion is done by a rolling window of small
global alignments (see tuning parameters in [Simulations](#simulations)) to find
*overlap* alignments of sequences in the database.

### Cycle breaking

The resulting overlap graph may not be a DAG due to two main reasons:

* wrong weak edges that should not exist.
* strong edges with the wrong direction.

The second case is typically caused by highly overlapping sequences (i.e the
start or end index of end points are too close). Currently such edges are
ignored altogether.

Regardless, cycle breaking is delegated to `igraph.Graph.feedback_arc_set` which
finds a set of edges the removal of which gives an acyclic graph.
It supports (see [docs](http://igraph.org/python/doc/igraph.GraphBase-class.html#feedback_arc_set))
an optimal, but slow (exponential complexity), integer programming algorithm
(presumably something similar to what is dicussed [here](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.31.5137))
and a suboptimal, but fast, algorithm relying on the [Eades heuristic](http://www.sciencedirect.com/science/article/pii/002001909390079O).

### Assembly in condensed alphabet

Assembly can be modified in two places to use condensed alphabets:

* *Indexing*: The sequence of reads can be indexed in the condensed alphabet.
  For example, if we are indexing 5-mers the read ``AAACCGTG`` gives only
  one tuple ``A3C2G1T1G1`` (which is 5 letters in the condensed alphabet).
  Typically, we may want to set the ``maxlen`` of the translator used for
  indexing to a very low number such that we do not miss seeds due to indels
  in long homopolymeric stretches. For example, if `maxlen` is set to 1 then
  the above example yields the tuple `A1C1G1T1G1`.
* *Seed extension*: This phase too can be performed in the condensed alphabet
  (and the translator may be a different one than the one for indexing, i.e
  have a different `maxlen`). This has the added benefit of allowing us to lower
  the penalty of homopolymeric indels.

### Simulations

For the simulated case where the true genome is known a *difference graph*
(which looks like a `diff`, with matching edges in black, missing edges in red,
and added edges in green) can be generated between the true layout path and
the assembled layout path. The key parameters for overlap discovery are:

* Window size for successive alignment frames,
* What constitutes a bad score in a single window,
* Number of consecutive bad scores which disqualifies a seed.
* Number of successful seeds (extending to boundaries) which is enough to call
  two reads overlapping (this is mainly an performance-optimization trick and
  does not seem to introduce errors).

Input generation parameters are:

* Length of the original genome,
* Parameters for the normal distribution of read lengths,
* Expected coverage.
* Substitution and gap probabilities used to mutate reads from true genome.


Usage:
```shell
# creates genome.assembly.fa, reads.assembly.fa, genome.assembly.db
make -f assembly.mk clean genome.assembly.db

# builds overlap.assembly.layout.gml, and compares against the true version.
make -f assembly.mk layout.diff.assembly.pdf
```

To perform assembly in condensed alphabet:
```shell
# creates genome.hp_assembly.fa, reads.hp_assembly.fa, genome.hp_assembly.db
make -f assembly.mk clean genome.hp_assembly.db MODE=hp_assembly

# builds overlap.hp_assembly.layout.gml, and compares against the true version.
make -f assembly.mk clean layout.diff.hp_assembly.pdf MODE=hp_assembly
```

### Behavior

**Good**

i. When compared to the true graph, the assembled overlap graph typically has
  some missing edges (e.g %15 of edges missing) but very few wrong edges are
  added (often none).
i. Generated overlap graphs are (close to) acyclic.
i. As a consequence of the (1), the assembled layout path is consistent
  with the true layout in the sense that the sequence of reads it announces
  as layout (its heaviest path) is a subsequence (i.e in correct order) of the
  correct layout path.

**Bad**

i. When two reads are both mostly overlapping the direction may come out wrong
  and this can cause cycles in the overlap graph.
i. There are occasional insertions too which do not seem to be problematic since
  they are weak (i.e low scoring alignments).
