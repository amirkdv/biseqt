# FIXME actually use the virtualenv
GCC = gcc -shared -Wall -fPIC -std=c99 -g
# To get coredumps:
# 	$ ulimit -c unlimited
# To debug coredump:
# 	$ gdb biseqt/libbiseqt.so core
LIBDIR = biseqt/pwlib
CFILES = $(LIBDIR)/_pw_internals.c $(LIBDIR)/pw.c $(LIBDIR)/seedext.c
$(LIBDIR)/pwlib.so: $(CFILES) $(LIBDIR)/pwlib.h
	$(GCC) $(CFILES) -o $@

clean:
	@find biseqt -regextype posix-extended -regex '.*.(pyc|swp|egg|egg-info)' | \
		grep -v '^./.git' | tee /dev/stderr  | while read f; do rm -rf $$f; done
	rm -rf build dist env
	rm -f core biseqt/pwlib/pwlib.so
	rm -rf docs/biseqt*.rst docs/_build docs/doxygen

tests: $(LIBDIR)/pwlib.so
	python -m biseqt.tests.pw

loc:
	find biseqt -type f -regex '.*\(\.py\|\.c\|\.h\)' | xargs wc -l

todo:
	find biseqt Makefile *.mk *.py -type f -regex '.*\(\.py\|\.c\|\.h\|Makefile\|\.mk\)' | xargs grep -A2 -nP --color 'FIXME|TODO|BUG'

env:
	# install numpy separately before, cf. setup.py
	virtualenv $@
	. env/bin/activate && \
		pip install numpy && \
		pip install -e . && \
		pip install -e .[docs]

DOCS_EXCLUDE = biseqt/tests
DOCS_OUT = _build
# in RTD docs/doxygen is built by docs/conf.py
docs: $(LIBDIR)/pwlib.so docs/docs.rst docs/doxygen
	sphinx-apidoc -e -o $@ biseqt/ $(DOCS_EXCLUDE)
	python docs/cdocs.py $(LIBDIR)/pwlib.h > docs/biseqt.pwlib.rst
	cd $@ && sphinx-build -b html . $(DOCS_OUT)
	cd $@ && sphinx-build -b latex . $(DOCS_OUT)
	@echo "Find the docs at file://`readlink -f $@/$(DOCS_OUT)/index.html`"

docs/doxygen:
	doxygen docs/doxygen.conf

docs/docs.rst:
	pandoc -f markdown -t rst -o $@ -i docs.md

DOCKER_IMG=amirkdv/biseqt-base
docker_build:
	cat Dockerfile | docker build -t $(DOCKER_IMG) -


.PHONY: clean tests *.pdf loc docs todo docs/docs.rst docs/doxygen
