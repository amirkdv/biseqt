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
number is the score of the aligment and what follows is the `opseq` (edit
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
to a caveat discussed below.

Substitution score matrices for the original alphabet can be translated into a
substitution score matrix for the (larger) condensed alphabet provided gap
parameters for homopolymeric indels are given.

Furthermore, an alignment `opseq` for condensed versions of two sequences can
be translated back to an alignment for the original sequences.

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
provided in `assembly`. All graph processing is delegated to [igraph](http://igraph.org/python/).
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
supports an optimal, but slow, integer programming algorithm and a suboptimal,
but fast, [heuristic](http://www.sciencedirect.com/science/article/pii/002001909390079O)
algorithm.

### Simulations

For the simulated case where the true genome is known a difference graph can be
generated between the true overlap path and the assembled overlap path.
The key parameters for overlap discovery are:

1. Window size for successive alignment frames,
2. What constitutes a bad score in a single window,
3. How many consecutive bad scores disqualifies a seed.

Input generation parameters are:

1. Length of the original genome,
2. Parameters for the normal distribution of read lengths,
3. Expected coverage.


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

#### Good

1. When compared to the true graph, the assembled overlap graph typically has
  some missing edges (e.g %15 of edges missing) but very few wrong edges are
  added (often none).
1. Generated overlap graphs are (close to) acyclic.
1. As a consequence of the (i), the assembled layout path is consistent
  with the true layout in the sense that its sequence of reads is a
  subsequence (i.e in correct order) of the correct layout path.

#### Bad

1. When two reads are both mostly overlapping the direction may come out wrong
  and this can cause cycles in the overlap graph.
1. There are occasional insertions too which do not seem to be problematic since
  they are weak (i.e low scoring alignments).

## To Do

* Perform assembly on condensed sequences.
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
