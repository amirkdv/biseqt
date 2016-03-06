# FIXME actually use the virtualenv
GCC = gcc -shared -Wall -fPIC -std=c99 -g
# To get coredumps:
# 	$ ulimit -c unlimited
# To debug coredump:
# 	$ gdb soda/libsoda.so core
soda/libsoda.so: soda/lib/pw.c soda/lib/overlap.c soda/lib/libsoda.h
	$(GCC) soda/pw.c soda/overlap.c -o $@

clean:
	@find . -regextype posix-extended -regex '.*.(pyc|swp|egg|egg-info)' | grep -v '^./.git' | tee /dev/stderr  | while read f; do rm -rf $$f; done
	rm -rf build dist env
	rm -f soda/libsoda.so core
	rm -rf docs/$(DOCS_OUT)

tests: soda/libsoda.so
	python -m soda.tests.homopolymeric
	python -m soda.tests.soda
	make -f assembly.mk clean
	make -f assembly.mk genome.assembly.db    overlap.assembly.layout.pdf    layout.diff.assembly.pdf
	make -f assembly.mk genome.hp_assembly.db overlap.hp_assembly.layout.pdf layout.diff.hp_assembly.pdf MODE=hp_assembly

loc:
	find soda -type f -regex '.*\(\.py\|\.c\|\.h\)' | xargs wc -l

CAIRO=$(shell python -c 'import site, os.path; print filter(lambda e:os.path.isdir(e + "/cairo"), site.getsitepackages())[0] + "/cairo"')
env:
	virtualenv --system-site-packages $@
	. env/bin/activate && python setup.py install
	ln -s "$(CAIRO)" env/lib/python2.7/cairo

DOCS_EXCLUDE = soda/tests
DOCS_OUT = _build
docs: soda/libsoda.so docs/docs.rst docs/doxygen
	sphinx-apidoc -e -o $@ soda/ $(DOCS_EXCLUDE)
	cd $@ && sphinx-build -b html . $(DOCS_OUT)
	cd $@ && sphinx-build -b latex . $(DOCS_OUT)
	@echo "Find the docs at file://`readlink -f $@/$(DOCS_OUT)/index.html`"

docs/doxygen:
	doxygen docs/doxygen.conf

docs/docs.rst:
	pandoc -f markdown -t rst -o $@ -i docs.md

.PHONY: clean tests *.pdf loc docs docs/docs.rst docs/doxygen
