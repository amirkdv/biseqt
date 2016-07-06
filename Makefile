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

FLAKE8_EXCLUDE = biseqt/overlap,biseqt/mapping.py
flake8:
	flake8 biseqt tests --exclude=$(FLAKE8_EXCLUDE)

tests: $(LIBDIR)/pwlib.so flake8
	cd tests && PYTHONPATH=.. py.test --cov biseqt --capture=no --verbose

loc:
	find biseqt tests -type f -regex '.*\(\.py\|\.c\|\.h\)' | xargs wc -l

todo:
	find biseqt Makefile *.mk *.py -type f -regex '.*\(\.py\|\.c\|\.h\|Makefile\|\.mk\)' | xargs grep -A2 -nP --color 'FIXME|TODO|BUG'

deps:
	# install numpy separately before, cf. setup.py
	pip install numpy && \
	pip install -e .  && \
	pip install -e .[docs] .[tests]

env:
	virtualenv $@

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
	rm -rf $@
	doxygen docs/doxygen.conf

docs/docs.rst:
	pandoc -f markdown -t rst -o $@ -i docs.md

DOCKER_IMG=amirkdv/biseqt-base
docker_build:
	cat Dockerfile | docker build -t $(DOCKER_IMG) -


.PHONY: clean tests *.pdf loc docs todo docs/docs.rst docs/doxygen flake8
