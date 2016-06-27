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
	@find . -regextype posix-extended -regex '.*.(pyc|swp|egg|egg-info)' | grep -v '^./.git' | tee /dev/stderr  | while read f; do rm -rf $$f; done
	rm -rf build dist env
	rm -f biseqt/pwlib/pwlib.so core
	find docs/$(DOCS_OUT) -maxdepth 1 -mindepth 1 | grep -v _static | xargs rm -rf
	rm -rf docs/doxygen

tests: $(LIBDIR)/pwlib.so
	python -m biseqt.tests.pw

loc:
	find biseqt -type f -regex '.*\(\.py\|\.c\|\.h\)' | xargs wc -l

todo:
	find biseqt Makefile *.mk *.py -type f -regex '.*\(\.py\|\.c\|\.h\|Makefile\|\.mk\)' | xargs grep -A2 -nP --color 'FIXME|TODO|BUG'

CAIRO=$(shell python -c 'import site, os.path; print filter(lambda e:os.path.isdir(e + "/cairo"), site.getsitepackages())[0] + "/cairo"')
env:
	virtualenv --system-site-packages $@
	. env/bin/activate && python setup.py install
	ln -s "$(CAIRO)" env/lib/python2.7/cairo

DOCS_EXCLUDE = biseqt/tests
DOCS_OUT = _build
# in RTD docs/doxygen is built by docs/conf.py
# FIXME depends on doxygen
# FIXME use sphinx entrypoints from virtualenv
# FIXME who is creating the latex dir?
# FIXME need to explicitly announce sphinx_rtd_theme as dependency, cf. https://github.com/EUDAT-B2SAFE/B2HANDLE/pull/63/commits
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

.PHONY: clean tests *.pdf loc docs todo docs/docs.rst docs/doxygen
