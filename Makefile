GCC = gcc -shared -Wall -fPIC -std=c99 -g

# To get coredumps:
# 	$ ulimit -c unlimited
# To debug coredump:
# 	$ gdb align/libalign.so core
align/libalign.so: align/libalign.c align/libalign.h
	$(GCC) align/libalign.c -o $@

clean:
	@find . -regextype posix-extended -regex '.*.(pyc|swp|egg|egg-info)' | grep -v '^./.git' | tee /dev/stderr  | while read f; do rm -rf $$f; done
	rm -rf build dist env
	rm -f align/libalign.so core
	rm -rf docs/$(DOCS_OUT)
	make -f assembly.mk clean

tests: align/libalign.so
	python -m align.tests.homopolymeric
	python -m align.tests.align
	make -f assembly.mk clean
	make -f assembly.mk genome.assembly.db    overlap.assembly.layout.pdf    layout.diff.assembly.pdf
	make -f assembly.mk genome.hp_assembly.db overlap.hp_assembly.layout.pdf layout.diff.hp_assembly.pdf MODE=hp_assembly

loc:
	find align -type f -regex '.*\(\.py\|\.c\|\.h\)' | xargs wc -l

CAIRO=$(shell python -c 'import site, os.path; print filter(lambda e:os.path.isdir(e + "/cairo"), site.getsitepackages())[0] + "/cairo"')
env:
	virtualenv --system-site-packages $@
	. env/bin/activate && python setup.py install
	ln -s "$(CAIRO)" env/lib/python2.7/cairo

DOCS_EXCLUDE = align/tests
DOCS_OUT = _build
docs: align/libalign.so docs/docs.rst
	sphinx-apidoc -f -e -o $@ align/ $(DOCS_EXCLUDE)
	cd $@ && sphinx-build -b html . $(DOCS_OUT)
	cd $@ && sphinx-build -b latex . $(DOCS_OUT)
	@echo "Find the docs at file://`readlink -f $@/$(DOCS_OUT)/index.html`"

docs/docs.rst:
	pandoc -f markdown -t rst -o $@ -i docs.md

.PHONY: clean tests *.pdf loc docs docs/docs.rst
