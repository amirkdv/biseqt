[![Documentation Status](https://readthedocs.org/projects/alignpy/badge/?version=latest)](http://alignpy.readthedocs.org/en/latest/?badge=latest)

[`align.py`](https://alignpy.readthedocs.org/) is a library of sequence alignment and
fragment assembly algorithms in python and C. To get started:

```shell
make env
source env/bin/activate
make clean tests
```

For genome assembly simulations:
```shell
# creates genome.assembly.fa, reads.assembly.fa, genome.assembly.db
make -f assembly.mk clean genome.assembly.db

# builds overlap.assembly.layout.gml, and compares against the true version.
make -f assembly.mk layout.diff.assembly.pdf
```

To simulate assembly in condensed alphabet:
```shell
# creates genome.hp_assembly.fa, reads.hp_assembly.fa, genome.hp_assembly.db
make -f assembly.mk clean genome.hp_assembly.db MODE=hp_assembly

# builds overlap.hp_assembly.layout.gml, and compares against the true version.
make -f assembly.mk clean layout.diff.hp_assembly.pdf MODE=hp_assembly
```

## To Do

* Move seed expansion from Python to C.
* Simulations:

    * Test on larger data sets (requires speedup).
    * Separate sanity tests from simulations; write sanity tests for individual
      parts of assembly.
    * Support hompolymeric-specific indel parameters in random generation of genome
      sequencing reads.
    * *Real* data: test against Leishmania dataset.

* Improvements:

    * When performing alignments in the condensed alphabet we currently have
      no choice but to score indels (i.e non-homopolymeric indels) identically.
      For example, a deletion of `A9T7` is scored exactly the same as a
      deletion of `A1T1`. This can clearly lead to wrong alignments, but a
      clear effect is not yet observed in simulations. Fixing this would require
      that `libalign` accept content-dependent scores for indels which is
      feasible but not trivial.
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
