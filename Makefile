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

tests: align/libalign.so genome.db
	python -m align.tests.align
	python -m align.tests.homopolymeric
	python -m align.tests.tuples

genome.db: genome.fa
	python -c 'from align.tuples import TuplesDB as DB; \
		from align.seq import Alphabet as A; \
		B = DB("genome.db", wordlen=10, alphabet=A("ACGT")); \
		B.initdb(); B.populate("reads.fa", lim=-1); B.index()'

genome.fa:
	python -c "from align.seq import make_sequencing_fixture as msf; msf('genome.fa', 'reads.fa', 'query.fa');"

.PHONY: clean tests
