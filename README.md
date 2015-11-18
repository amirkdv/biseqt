[![Documentation Status](https://readthedocs.org/projects/alignpy/badge/?version=latest)](http://alignpy.readthedocs.org/en/latest/?badge=latest)

`align.py` is a library of sequence alignment and fragment assembly algorithms
in python and C. To get started:

```shell
python setup.py develop
make clean tests
```

## Random processes

The following can be generated randomly by the `seq` module:

* Sequence of given length given a letter probabilities.
* Mutants of given sequence given substitution and gap probabilities.
* Random reads from a given sequence given required coverage, read length
  distribution, and substitution/gap probabilities.

## Pairwise alignments

The core dynamic programming algorithm is implemented [in C](https://github.com/amirkdv/align.py/blob/master/align/libalign.c)
and interfaced using [cffi](https://cffi.readthedocs.org/en/latest/). The
following are supported by `libalign`:

* *Global* alignments, i.e Needleman-Wunsch.
* *Local* alignments, i.e Smith-Waterman.
* *Overlap* alignments where only suffix-prefix alignments are considered.
* *Anchored* alignments where only alignments starting/ending at the
  beginning/end of the two sequences are considered,
* *Affine gap* penalties are supported.
* *Banded* alignments are supported.
* Alignment can be limited to a *frame* for each sequence.

### Scoring

Gap and substitution probabilities can be translated by `align.AlignParams` to
corresponding log likelihood scores. The random (null hypothesis) distribution
of letters must be provided.

Scoring an alignment is performed on the corresponding `align.Transcript` object
and the original `seq.Sequence`'s. Each `align.Transcript` looks like this:
```
(100,12),12.50:MMMISSSMMMMDDD
```
The first two integers are starting positions of the alignment, the following
number is the score of the alignment and what follows is the `opseq` (edit
transcript) for the alignment.

### Complexity

**Time** complexity is quadratic in sequence lengths except for banded global
alignment which is linear with a constant proportional to band width. As an
example, finding the optimal local alignment for two related sequences (one is
an artificial mutant of the other) with length 5 Kbp takes on
average 3 seconds on a 3.67 GHz quad-core Intel CPU.

**Space** complexity is quadratic too. Currently an average of roughly 100 bytes
is required per cell of the dynamic programming table. For example, it takes
more than 6GB of space to align two sequences of length 8000 (See [To
Do](#to-do)).

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

Substitution score matrices for the original alphabet can be translated into a
substitution score matrix for the (larger) condensed alphabet given the
following pieces of information:

i. Gap probability in homopolymeric stretches (only a linear gap model is well
  specified),
i. Gap probability in other regions (here an affine gap model is fine),
i. Substitution probabilities of letters in original alphabet.

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
provided in `assembly` which delegates all graph algorithms to [igraph](http://igraph.org/python/).
A weighted, DAG is built by seed expansion (see
[Tuples Methods](#tuples)) on all pairs of sequences and the longest
path is reported as the layout. Expansion is done by `assembly.OverlapBuilder`
which uses a rolling window of small global alignments (see tuning parameters
in [Simulations](#simulations)) to find *overlap* alignments of sequences in the
database.

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
  For example, if we are indexing all 5-mers the read ``AAACCGTG`` gives only
  one tuple ``A3C2G1T1G1`` (which is 5 letters in the condensed alphabet).
  Typically, we may want to set the ``maxlen`` of the translator used for
  indexing to a very low number such that we do not miss seeds due to indels
  in long homopolymeric stretches.
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
# creates genome.fa, reads.fa, genome.db
make -f assembly.mk genome.db

# builds overlap.dag.gml, overlap.layout.gml, and compares against true versions
make -f assembly.mk layout.diff.assembly.pdf
```

To perform assembly in condensed alphabet:
```shell
make clean
make -f assembly.mk layout.diff.hp_assembly.pdf MODE=hp_assembly
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

## To Do

* Move seed expansion from Python to C.
* Simulations:

    * Test on larger data sets (requires speedup).
    * Separate sanity tests from simulations; write sanity tests for individual
      parts of assembly.
    * Support hompolymeric-specific indel parameters in random generation of genome
      sequencing reads.
    * *Real* data: test against Leishmania dataset.

* Code:

    * Make `align.Transcript` a `namedtuple` as well (unless it's becoming a
      `CffiObject`).

* Improvements:

    * An overlap graph must satisfy two consistency criteria: it is a DAG,
      and for any vertex *u* in it, any pair of outgoing (incoming) neighbors of
      *u* are adjacent.  Assembly overlap graphs are DAG (or close to it) but
      they rarely satisfy the second. The second criteria can be used to find
      missing edges by brute force overlap alignment (this matches the typical
      case of left-out-vertices in simulations). The difficulty is to find a way
      to recover necessary edges for a full layout path without trying to
      recover *all* missing edges.
    * Stop ignoring sequence pairs that are
      mostly overlapping. These are currently ignored since we may get the
      direction wrong on a heavy edge.

* Low priority:

    * Figure out how to pull in docstrings from C code into sphinx (e.g look
      at [Breathe](https://github.com/michaeljones/breathe)).
    * Add an ungapped seed expansion phase.
    * Adapt Karlin-Altschul statistics (references:
      [[1]](http://www.pnas.org/content/87/6/2264.full.pdf),
      [[2]](https://publications.mpi-cbg.de/Altschul_1990_5424.pdf),
      [[3]](http://www.jstor.org/stable/1427732?seq=1#page_scan_tab_contents), and
      chap. 7-9 [[4]](https://books.google.ca/books?id=uZvlBwAAQBAJ)) to the
      problem of finding overlaps.
    * Support [Hirschberg](https://en.wikipedia.org/wiki/Hirschberg\'s_algorithm) -style
      linear space optimization in `libalign`.
    * Make it work with Python 3.
