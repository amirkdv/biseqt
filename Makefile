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

tests: align/libalign.so genome.db overlap_true.svg
	python -m align.tests.align
	python -m align.tests.homopolymeric
	python -m align.tests.tuples

overlap_true.svg: genome.db
	python -c '\
		from align.tuples import TuplesDB as DB; \
		from align.seq import Alphabet as A; \
		from align.assembly import save_overlap_graph, overlap_graph_by_known_order; \
		D = DB(db="genome.db", wordlen=10, alphabet=A("ACGT")); \
		G = overlap_graph_by_known_order(D); \
		save_overlap_graph(G, "overlap_true.svg", figsize=(100,100))'

genome.db: align/libalign.so genome.fa
	python -c '\
		from align.tuples import TuplesDB as DB; \
		from align.seq import Alphabet as A; \
		B = DB("genome.db", wordlen=10, alphabet=A("ACGT")); \
		B.initdb(); B.populate("reads.fa", lim=-1); B.index()'

genome.fa: align/libalign.so
	python -c '\
		from align.seq import make_sequencing_fixture as msf; \
		msf("genome.fa", "reads.fa", genome_length=1000, coverage=5, len_mean=100, len_var=25);'

.PHONY: clean tests
