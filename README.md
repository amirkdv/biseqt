[![Documentation Status](https://readthedocs.org/projects/alignpy/badge/?version=latest)](http://alignpy.readthedocs.org/en/latest/?badge=latest)

[`align.py`](https://alignpy.readthedocs.org/) is a library of sequence alignment and
fragment assembly algorithms in python and C. To get started:

```shell
make env
source env/bin/activate
make clean tests
```

## To Do

* Improvements:
    * Move seed expansion from Python to C.
    * The alignment problem in the condensed alphabet seems ill-defined as it
      currently stands. A clear example of this is the fact that we currently
      don't have a way to properly align `AAAAAA` and `AAACCC`: our best option
      is `A6--` and `A3C3`. A possible complicated formulation is to allow
      single letters to be matched to multiple letters in an alignment:
      this requires allowing nonstandard choices in the DP table. There is a
      way of doing this while maintaining polynomial time complexity (at most
      cubic) but the implementation is not trivial.
    * When performing alignments in the condensed alphabet we currently have
      no choice but to score indels (i.e non-homopolymeric indels) identically.
      For example, a deletion of `A9T7` is scored exactly the same as a
      deletion of `A1T1`. This can clearly lead to wrong alignments, but a
      clear effect is not yet observed in simulations. Fixing this would require
      that `libalign` accept content-dependent scores for indels which is
      feasible but not trivial.
    * An overlap graph must satisfy two consistency criteria:
      * it is a DAG, and
      * for any vertex *u* in it, any pair of outgoing (incoming) neighbors of
        *u* are adjacent.

      Assembly overlap graphs are DAG (or close to it) but
      they rarely satisfy the second. The second criteria can be used to find
      missing edges by brute force overlap alignment (this matches the typical
      case of left-out-vertices in simulations). The difficulty is to find a way
      to recover necessary edges for a full layout path without trying to
      recover *all* missing edges.
    * Stop ignoring sequence pairs that are
      mostly overlapping. These are currently ignored since we may get the
      direction wrong on a heavy edge.

* Low priority:
    * Support [Hirschberg](https://en.wikipedia.org/wiki/Hirschberg\'s_algorithm) -style
      linear space optimization in `libalign`.
    * Add an ungapped seed expansion phase.
    * Adapt Karlin-Altschul statistics (references:
      [[1]](http://www.pnas.org/content/87/6/2264.full.pdf),
      [[2]](https://publications.mpi-cbg.de/Altschul_1990_5424.pdf),
      [[3]](http://www.jstor.org/stable/1427732?seq=1#page_scan_tab_contents), and
      chap. 7-9 [[4]](https://books.google.ca/books?id=uZvlBwAAQBAJ)) to the
      problem of finding overlaps.
    * Make it work with Python 3.
