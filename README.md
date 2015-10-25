# align.py
This is a library of some sequence alignment algorithms in python and C. To get
started:

```sh
python setup.py develop
make clean tests
```

## Pairwise alignments

The core dynamic programming algorithm is [implemented](/align/libalign.c) in C
and [interfaced](/align/align.py) using
[`cffi`](https://cffi.readthedocs.org/en/latest/) The following are implemented;
see the [tests](/align/tests/align.py) for example usage:

* **Global** alignments, i.e Needleman-Wunsch.
* **Local** alignments, i.e Smith-Waterman.
* **Overlap** alignments where only suffix-prefix alignments are considered.
* **Anchored** alignments where only alignments starting/ending at the
  beginning/end of the two sequences are considered,
* **Affine** gap penalties and **Banded** alignments are supported as well.

*Notes*:

3. Alignments are represented (and transported from C to Python) in a special
  format which looks like this:

      (100,12),12.50:BMMMISSSMMMMDDD

  The first two integers are starting positions of the alignment, the following
  number is the score of the aligment and what follows is an edit transcript
  (look for [`opseq`](/align/align.py)) for the alignment.
2. Time complexity is quadratic in sequence lengths except for banded global
  alignment which is linear with a constant proportional to
  band width. As an example, finding the optimal local alignment for two related
  sequences (one is an [artificial mutant](/align/seq.py) of the other)
  with length 5 Kbp takes on average 3 seconds on a 3.67 GHz quad-core Intel
  CPU.
1. Space complexity is quadratic. Currently an average of roughly 100 bytes is
  required per DP table cell. For example, it takes more than 6GB of space to
  align two related sequences of length 8000. See [todo list](#missing) below.


## Tuple methods

[Tuple methods](/align/tuples.py) rely on
[python-sqlite](https://docs.python.org/2/library/sqlite3.html) to store, index,
and query tuples through SQL (also parts of the algorithm to find seeds is
delegated, via SQL, to SQLite). The following operations are implemented; see the
[tests](/align/tests/tuples.py) for example usage:

* Loading a FASTA database of sequences (parsing is done via
  [`Bio.SeqIO`](http://biopython.org/wiki/SeqIO)) into the SQLite database.
* Indexing sequences with an arbitrary word (tuple) length.
* Querying the index for sequences in order of number of hits they share with a
  query string.
* Finding seeds for a given query string by collapsing overlaping hits.
* (**incomplete**) Extending seeds by repeated anchored alignments.

## Alphabet translation

A set of experimental methods are [implemented](/align/distillery.py) to
condense sequences into an alternative alphabet where each homopolymeric
substring is translated to a single letter(e.g `AACCC` becomes `A2C3`). The
following are implemented; see the [tests](/align/tests/distillery.py) for
example usage:

* Translating a sequence from an arbitrary alphabet to a distilled alphabet.
  This requires specifying a whole number `maxlen`:
  homopolymeric substrings longer than `maxlen` are considered to have only
  `maxlen` characters. This is needed to ensure constant letter length across
  the distilled alphabet.
* Translating a substitution score matrix for the original alphabet into a
  substitution score matrix for the (larger) distilled alphabet. This requires
  specifying an affine gap penalty scheme for homopolymeric indels.
* Expanding a distilled sequence back to the original sequence.
* Translating an alignment `opseq` for distilled versions of two sequences back
  to an alignment for the original sequences.

*Notes*:

1. If source alphabet sequences contain homopolymeric substrings that are
  longer than the specified `maxlen`, the distilling process is lossy (expanding
  a distilled sequence does not necessarily give its original sequence).
  Translating edit transcripts (i.e `opseq`s) does not have an issue with this.

## Missing
* complete [seed expansion](/align/tuples.py)
* make tuple methods aware of distilled sequences.
* support hompolymeric-specific indel parameters in [random generation](/align/seq.py)
  of whole genome reads.
* better [tests](/tests).
* support [Hirschberg](https://en.wikipedia.org/wiki/Hirschberg's_algorithm)-style
[linear space](/align/libalign.c) optimization.
* make it work with Python 3.
