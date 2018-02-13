GCC = gcc -shared -Wall -fPIC -std=c99 -g
# To get coredumps:
# 	$ ulimit -c unlimited
# To debug coredump:
# 	$ gdb biseqt/libbiseqt.so core
LIBDIR = biseqt/pwlib
CFILES = $(LIBDIR)/_pw_internals.c $(LIBDIR)/pw.c
$(LIBDIR)/pwlib.so: $(CFILES) $(LIBDIR)/pwlib.h
	$(GCC) $(CFILES) -o $@

clean:
	@find biseqt tests experiments -regextype posix-extended -regex '.*.(pyc|swp|swo|egg|egg-info)' | \
		grep -v '^./.git' | tee /dev/stderr  | while read f; do rm -rf $$f; done
	rm -f core biseqt/pwlib/pwlib.so
	rm -rf docs/_build docs/doxygen

FLAKE8_INCLUDE = biseqt tests # experiments
FLAKE8_EXCLUDE =
flake8:
	flake8 $(FLAKE8_INCLUDE) --exclude=$(FLAKE8_EXCLUDE)

PYTEST_OPTS=--cov biseqt --capture=no --cov-report=term-missing --verbose
COVERAGE_BADGE=coverage.svg  # relative to tests/
tests: $(LIBDIR)/pwlib.so flake8
	rm -f tests/$(COVERAGE_BADGE)
	cd tests && PYTHONPATH=.. py.test $(PYTEST_OPTS) && coverage-badge -o $(COVERAGE_BADGE)

experiments:
	cd $@ && PYTHONPATH=.. python blot_stats.py

loc:
	find biseqt tests -type f -regex '.*\(\.py\|\.c\|\.h\)' | grep -v __pycache__ | xargs wc -l

todo:
	find biseqt Makefile *.mk *.py -type f -regex '.*\(\.py\|\.c\|\.h\|Makefile\|\.mk\)' | xargs grep -A2 -nP --color 'FIXME|TODO|BUG'

env:
	virtualenv $@

DOCS_OUT = _build
# in RTD docs/doxygen is built by docs/conf.py
docs: $(LIBDIR)/pwlib.so docs/doxygen
	rm -rf docs/$(DOCS_OUT)
	# sphinx-apidoc -e -o $@ biseqt/
	docs/doxygen-glue.py $(LIBDIR)/pwlib.h > docs/biseqt.pwlib.rst
	cd $@ && sphinx-build -b html . $(DOCS_OUT)
	cd $@ && sphinx-build -b latex . $(DOCS_OUT)
	@echo "Find the docs at file://`readlink -f $@/$(DOCS_OUT)/index.html`"

gh-coverage-badge:
	@echo "Updating github's cache of test coverage badge"
	curl https://github.com/amirkdv/biseqt/blob/master/README.md 2>/dev/null | \
		grep coverage.svg | \
		sed 's/.*img src="\([^"]*\).*$$/\1/' | \
		xargs curl -X PURGE

docs/doxygen:
	rm -rf $@
	mkdir $@
	doxygen docs/doxygen.conf

DOCKER_IMG=amirkdv/biseqt-base
docker_build:
	cat Dockerfile | docker build -t $(DOCKER_IMG) -


.PHONY: clean tests *.pdf loc docs todo experiments docs/doxygen flake8
