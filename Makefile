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
	rm -f *.db
	rm -f test.fa query.fa

tests: align/libalign.so test.fa
	python -m align.tests.align
	python -m align.tests.distillery
	python -m align.tests.tuples

test.fa:
	python -c "import sys; from align.seq import Alphabet; A = Alphabet('ACGT'); x = A.randseq(160);p=[[0.94 if k==i else 0.02 for k in range(4)] for i in range(4)]; print '> orig\n' + str(x); print '\n'.join(['>seq_{}\n'.format(i) + str(x.mutate(go_prob=0.1, ge_prob=0.5, subst_probs=p)[0]) for i in range(10)]);" > /tmp/$@
	head -n +2 /tmp/$@ | tail -n1 > query.fa
	tail -n +3 /tmp/$@ > $@
	rm -f /tmp/$@

.PHONY: clean tests
