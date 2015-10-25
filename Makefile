GCC = gcc -shared -Wall -fPIC -std=c99 -g

align/libalign.so: align/libalign.c align/libalign.h
	$(GCC) align/libalign.c -o $@

# To get coredumps: ulimit -c unlimited
# To debug coredump: gdb align.so core
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
	python -c "from align.seq import randgen, mutate; x = randgen(60); print '> orig\n' + x; print '\n'.join(['> seq_%d\n' % i + mutate(x, gap_open=0.1, rates={i:{k:0.94 if k==i else 0.02 for k in 'ACGT'} for i in 'ACGTG'})[0] for i in range(10)]);" > /tmp/$@
	head -n +2 /tmp/$@ | tail -n1 > query.fa
	tail -n +3 /tmp/$@ > $@
	rm -f /tmp/$@

.PHONY: clean tests
