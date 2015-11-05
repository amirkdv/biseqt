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

tests: align/libalign.so genome.db overlap_align.svg overlap_tuple.svg overlap_true.svg
	python -m align.tests.homopolymeric
	python -m align.tests.align

overlap_align.svg:
	python -c 'import align.tests.assembly as T; T.overlap_by_alignment("genome.db", "$@")'

overlap_tuple.svg:
	python -c 'import align.tests.assembly as T; T.overlap_by_tuple_extension("genome.db", "$@")'

genome.db: align/libalign.so
	python -c 'import align.tests.assembly as T; T.create_example("$@")'
	python -c 'import align.tests.assembly as T; T.overlap_by_known_order("genome.db", "overlap_true.svg")'

.PHONY: clean tests *.svg
