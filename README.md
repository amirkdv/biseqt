[![Documentation Status](https://readthedocs.org/projects/alignpy/badge/?version=latest)](http://alignpy.readthedocs.org/en/latest/?badge=latest)

[`align.py`](https://alignpy.readthedocs.org/) is a library of sequence alignment and
fragment assembly algorithms in python and C. To get started:

```shell
make env
source env/bin/activate
make clean tests
```

## To Do

* Merge the functionality of ad hoc scripts: `num_seeds.py`, `spectra.py`,
  `prepare.py` and `rw.py`.
* Missing docs:
  * p-value calculation for words which discards potentially repetitive
    structures.
  * p-value calculation for shift distributions based on 2d representation.
  * improved algorithm for condensing seeds.
* Allow the same algorithm to be used for aligning noisy long reads against a
  reference genome; the current scheme (using `bwa`) does not seem accurate
  enough. Additionally:
  * We are currently setting the end position of a read as the sum of its length
    and its start position. This is clearly inaccurate due to high indel rates.
  * The idea is this: the starting segments for extension can be null (e.g. a
    global alignment problem can be solved by starting with a null segment
    starting at `(0,0)` and allowing arbitrarily many falls in the score. This
    can be applied to the overlap alignment problem by introducing many null
    segments and trying to extend all of them with a reasonable `max_new_mins`.
    This can be further sped up using the hint given by shift distribution which
    is very strong and accurate for aligning reads against a reference genome.
* Stop cheating with reverse complements.
* Once we have a robust way of aligning reads against a reference genome we can
  judge whether *all* seeds are necessary for seed extension or those very close
  to the shift mode.
* Make it work with Python 3.
* The alignment problem in the condensed alphabet seems ill-defined as it
  currently stands. A clear example of this is the fact that we currently
  don't have a way to properly align `AAAAAA` and `AAACCC`: our best option
  is `A6--` and `A3C3`. A possible complicated formulation is to allow
  single letters to be matched to multiple letters in an alignment:
  this requires allowing nonstandard choices in the DP table. There is a
  way of doing this while maintaining polynomial time complexity (at most
  cubic) but the implementation is not trivial.
