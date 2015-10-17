# align.py
This is a library of some sequence alignment algorithms in python and C. To get
started:

```sh
make clean tests
```

## Pairwise alignments

The core dynamic programming algorithm is [implemented](/align/libalign.c) in C
and [interfaced](/align/align.py) using
[`cffi`](https://cffi.readthedocs.org/en/latest/) The following are supported;
see the [tests](/align/tests/align.py) for example usage:

* *Needleman-Wunsch* for global alignments,
* *Smith-Waterman* for local alignments,
* *Overlap* alignments: only suffix-prefix alignments are considered,
* *Anchored* alignments: only alignments starting/ending at the beginning/end
of the two sequences are considered,
* *Banded* alignments: any of the above alignments can optionally be performed
with a specified band width,
* *Affine* gap penalty: any of the above alignments can use an affine gap penalty.

Alignments are represented (and transported from C to Python) in a special
format which looks like this: `(100,12),12.50:BMMMISSSMMMMDDD`. The first two
integers are starting positions of the alignment, the following number is the
score of the aligment and what follows is an [`opseq`](/align/align.py) for the
alignment.

## Tuple methods

Tuple methods rely on [python-sqlite](https://docs.python.org/2/library/sqlite3.html)
to store, index, and query tuples through SQL. The following are supported; see
the [tests](/align/tests/tuples.py) for example usage:

* Loading a FASTA database of sequences (parsing is done via
  [`Bio.SeqIO`](http://biopython.org/wiki/SeqIO)) into the SQLite database.
* Indexing the database of sequences for any word (tuple) length.
* Finding and extending alignment seeds (**incomplete**) for a query sequence.

## Alphabet translation

A set of experimental methods are [avialable](/align/distillery.py) to condense
sequences into an alternative alphabet where each homopolymeric substring is
translated to a single character (e.g `AACCC` becomes `A2C3`). The alignment
algorithms support alphabets with arbitrary letter sizes (but letter sizes must
be fixed across the alphabet). The following are supported; see the
[tests](/align/tests/distillery.py) for example usage:

* Translating a sequence from an arbitrary alphabet to a distilled alphabet.
  *Note*: distilling a sequence requires specifying a whole number `maxlen` for
  homopolymeric substrings where longer substrings are considered as if they had
  `maxlen` characters. This is needed to ensure constant letter length across
  the distilled alphabet.
* Translating a substitution score matrix for the original alphabet into a
  substitution score matrix for the (larger) distilled alphabet. The new matrix
  will have dimensions `maxlen` × A where A is the original alphabet size. The
  new score matrix is derived by assuming a separate affine gap penalty scheme
  for homopolymeric indels (i.e fixes probability of gap opening and geometric
  distribution of gap extension).
* Expanding a distilled sequence back to the original sequence. *Note*: if
  the original sequence contains homopolymeric substrings that are longer than
  the `maxlen` this process is not the exact inverse of the distilling process.
* Translating an alignment opseq generated for the distilled versions of two
  sequences to an alignment for the original sequences.

## Missing
* a `setup.py` and documentation for `python-cffi` backend installation.
* better [tests](/tests).
* complete [seed expansion](/align/tuples.py)
* support [random generation](/align/randseq.py) of assembly reads.
* support [Hirschberg](https://en.wikipedia.org/wiki/Hirschberg's_algorithm)-style
[linear space](/align/libalign.c) optimization. Currently an average of 100
bytes is required per DP table cell. Consequently, with the current quadratic
space implementation it takes more than 6GB of space to align two sequences of
length 8000.
* make it work with Python 3.
