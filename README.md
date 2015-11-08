# align.py
This is a library of some sequence alignment algorithms in python and C. To get
started:

```sh
python setup.py develop
make clean tests
```

## Random processes

The following can be generated randomly by the `seq` module:

* Sequence of given length given a letter probabilities.
* Mutants of given sequence given substition and gap probabilities.
* Random reads from a given sequence given required coverage, read length
  distribution, and substitution/gap probabilties.

## Pairwise alignments

The core dynamic programming algorithm is implemented [in C](/align/libalign.c)
and interfaced (via `align.align`) using
[`cffi`](https://cffi.readthedocs.org/en/latest/). The following are supported
by `libalign`:

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

Some tuple methods (aka *l*-mer analysis) are provided by `tuples.TupleDB`
which is backed by [SQLite](https://docs.python.org/2/library/sqlite3.html) to
store, index, and query tuples and by
[`Bio.SeqIO`](http://biopython.org/wiki/SeqIO) to read FASTA files.

### Alignment

For any two given sequences (both already indexed) `tuples.OverlapFinder` can be
used to:
* Find maximal exactly matching "seeds".
* Find out if a seed can be extended to a suffix-prefix alignment by repeated
  short global alignments of a fixed window size.

Note that seed extension is specifically geared towards the genome assembly
problem and, unlike BLAST, it does not try to find *all* significant local
alignments, but only those that would correspond to a suffix-prefix alignment.

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
sequence).

## Genome assembly

Overlap and layout graphs (i.e OLC minus consensus) can be calculated by methods
provided in `assembly`. A weighted, DAG is built by seed expansion on all pairs
of sequences and the longest path is reported as the layout.

### Simulations

For the simulated case where the true genome is known a difference graph can be
generated between the true overlap path and the assembled overlap path.
The key paramters for overlap discovery are:

1. Window size for successive alignment frames,
2. What constitutes a bad score in a single window,
3. How many consecutive bad scores disqualifies a seed.

Input generation parameters are:

1. Length of the original genome,
2. Parameters for the normal distribution of read lengths,
3. Expected coverage.


Usage:
```
$ make clean genome.db                # creates genome.fa, reads.fa, genome.db
$ make true_overlap.gml true_overlap.layout.pdf # find the true overlap and layout
$ make overlap.gml overlap.layout.pdf # find the overlap and layout by seed extension
$ make overlap.layout.diff.pdf        # diff against the true overlap graph
```

### Behavior

* Good:
  * When compared to the true graph, the assembled overlap graph typically has
    some missing edges (e.g %15 of edges missing) but very few wrong edges are
    added (often none).
  * Generated overlap graphs are (close to) acyclic.
* Bad:
  * When two reads are both mostly overlapping the direction may come out wrong and
    this can cause cycles in the overlap graph.
  * There are occassional insertions too which do not seem to be problematic since
    they are weak (i.e low scoring alignments).

## To Do

#### Missing features
* Perform assembly on condensed sequences.
* Deal with potential cycles. There are two issues:
  * weak edges: easy to deal with (remove weak edges until the overlap graph is
    acyclic).
  * wrong direction strong edges: see below.

#### Improvements
* Code docs: all of [`OverlapFinder`](/align/tuples.py) and [Assembler](/align/assembly.py).
* For any two reads, do we need to pursue all segments that satisfy the
  score criteria or should we drop out once we find one segment? Note that most
  of the time for overlapping sequences many seeds come from the same correct
  suffix-prefix alignment.
* Deal with the case where two sequences mostly overlap (with close start
  and ends). In such cases a small error in alignment (which would have a small
  score cost) can reverse the direction of the edge. This is problematic because
  such edges are typically really heavy (due to the strong overlap). Currently,
  if the starting or ending indices of an alignment are too close (less than 5)
  the alignment is ignored. Alternatively, we could perform a full global
  suffix-prefix alignment to keep the edge but typically such edges are not
  really informative.
* Make `align.Transcript` a `namedtuple` as well (unless it's becoming a
  `CffiObject`).
* An overlap graph must satisfy two consistency criterions:
  * it is a DAG,
  * For any vertex *u* in it, any pair of outgoing (incoming) neighbors of
    *u* are adjacent.

  Assembly overlap graphs are DAG (or close to it) but they rarely satisfy the
  second. The second criteria can be used to find missing edges by brute force
  overlap alignment (this matches the typical case of left-out-vertices in
  simulations). The difficulty is to find a way to recover necessary edges
  for a full layout path without trying to recover *all* missing edges.
* Move seed expansion from Python to C.

#### Simulations
* Make sure random reads cover the entire genome.
* Test on larger data sets (requires speedup).
* Separate sanity tests from simulations; write sanity tests for individual
  parts of assembly.
* Support hompolymeric-specific indel parameters in random generation of genome
  sequencing reads.
* *Real* data: test against Leishmania dataset.

#### Low priority
* Add an ungapped seed expansion phase.
* Adapt Karlin-Altschul statistics (references:
  [[1]](http://www.pnas.org/content/87/6/2264.full.pdf),
  [[2]](https://publications.mpi-cbg.de/Altschul_1990_5424.pdf),
  [[3]](http://www.jstor.org/stable/1427732?seq=1#page_scan_tab_contents), and
  chap. 7-9 [[4]](https://books.google.ca/books?id=uZvlBwAAQBAJ)) to the
  problem of finding overlaps.
* Support [Hirschberg](https://en.wikipedia.org/wiki/Hirschberg's_algorithm)-style
  linear space optimization (cf. [`libalign::solve()` and `libalign::tracback()`](/align/libalign.c)).
* Make it work with Python 3.
