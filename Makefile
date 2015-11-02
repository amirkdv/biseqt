GCC = gcc -shared -Wall -fPIC -std=c99 -g

# To get coredumps:
# 	$ ulimit -c unlimited
# To debug coredump:
# 	$ gdb align/libalign.so core
align/libalign.so: align/libalign.c align/libalign.h
	$(GCC) align/libalign.c -o $@

clean:
	find . -regex .\*.pyc | while read f; do rm -f $$f; done
	find . -regex .\*.swp | while read f; do rm -f $$f; done
	rm -f align/libalign.so
	rm -f core
	rm -f genome.fa reads.fa genome.db
	rm -f overlap.svg

tests: overlap_align.svg overlap_tuple.svg overlap_true.svg
	python -m align.tests.homopolymeric
	python -m align.tests.align

overlap_align.svg: genome.db
	python -c 'import align.tests.tuples as T; T.overlap_by_alignment("genome.db", "$@")'

overlap_tuple.svg: genome.db
	python -c 'import align.tests.tuples as T; T.overlap_by_tuple_extension("genome.db", "$@")'

overlap_true.svg: genome.db
	python -c 'import align.tests.tuples as T; T.overlap_by_known_order("genome.db", "$@")'

genome.db: align/libalign.so
	python -c 'import align.tests.tuples as T; T.create_example("$@")'

.PHONY: clean tests *.svg
