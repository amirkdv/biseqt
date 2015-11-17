GCC = gcc -shared -Wall -fPIC -std=c99 -g
TRUE_GRAPH = true_overlap
ASSEMBLED_GRAPH = overlap
DB = genome.db
ASSEMBLY_TEST = align.tests.assembly

# To get coredumps:
# 	$ ulimit -c unlimited
# To debug coredump:
# 	$ gdb align/libalign.so core
align/libalign.so: align/libalign.c align/libalign.h
	$(GCC) align/libalign.c -o $@

clean:
	find . -regex .\*.pyc | while read f; do rm -rf $$f; done
	find . -regex .\*.swp | while read f; do rm -rf $$f; done
	find . -regex .\*.egg | while read f; do rm -rf $$f; done
	find . -regex .\*.egg-info | while read f; do rm -rf $$f; done
	rm -rf env build dist
	rm -f align/libalign.so core
	rm -f genome.fa reads.fa $(DB)
	rm -f *.gml *.pdf
	rm -rf docs/$(DOCS_OUT)

tests: align/libalign.so $(DB) $(TRUE_GRAPH).gml $(ASSEMBLED_GRAPH).gml layout_diff.pdf
	python -m align.tests.homopolymeric
	python -m align.tests.align

$(ASSEMBLED_GRAPH).gml: $(TRUE_GRAPH).gml
	python -c 'import $(ASSEMBLY_TEST) as T; T.overlap_by_seed_extension("$(DB)", "$@")'
	python -c 'import align.assembly as A, igraph as ig, sys; A.OverlapGraph(ig.read("$(TRUE_GRAPH).gml")).diff_text(A.OverlapGraph(ig.read("$(ASSEMBLED_GRAPH).gml")), sys.stdout)'

$(ASSEMBLED_GRAPH).dag.gml: $(ASSEMBLED_GRAPH).gml
	python -c 'import align.assembly as A, igraph as ig; G = A.OverlapGraph(ig.read("$(ASSEMBLED_GRAPH).gml")); G.break_cycles(); G.save("$@")'

$(TRUE_GRAPH).gml:
	python -c 'import $(ASSEMBLY_TEST) as A, align.seq as S, align.tuples as T; A.overlap_graph_by_known_order("$(DB)").save("$@")'

$(ASSEMBLED_GRAPH).pdf:
	python -c 'import align.assembly as A, igraph as ig; A.OverlapGraph(ig.read("$(ASSEMBLED_GRAPH).gml")).draw("$@");'

$(TRUE_GRAPH).pdf: $(TRUE_GRAPH).gml
	python -c 'import align.assembly as A, igraph as ig; A.OverlapGraph(ig.read("$(TRUE_GRAPH).gml")).draw("$@");'

# The true graph doesn't have any cycles, no need for a .dag.gml.
$(TRUE_GRAPH).layout.gml: $(TRUE_GRAPH).gml
	python -c 'import align.assembly as A, igraph as ig; A.OverlapGraph(ig.read("$(TRUE_GRAPH).gml")).layout().save("$@")'

$(ASSEMBLED_GRAPH).layout.gml: $(ASSEMBLED_GRAPH).dag.gml
	python -c 'import align.assembly as A, igraph as ig; A.OverlapGraph(ig.read("$(ASSEMBLED_GRAPH).dag.gml")).layout().save("$@")'


$(TRUE_GRAPH).layout.pdf: $(TRUE_GRAPH).gml
	python -c 'import align.assembly as A, igraph as ig; A.OverlapGraph(ig.read("$(TRUE_GRAPH).gml")).draw("$@", longest_path=True);'

# When drawing the layout show all paths; the .gml file contains the longest path only.
$(ASSEMBLED_GRAPH).layout.pdf: $(ASSEMBLED_GRAPH).layout.gml
	python -c 'import align.assembly as A, igraph as ig; g = A.OverlapGraph(ig.read("$(ASSEMBLED_GRAPH).dag.gml")); g.draw("$@", highlight_paths=g.all_longest_paths());'

layout_diff.pdf: $(ASSEMBLED_GRAPH).layout.gml $(TRUE_GRAPH).layout.gml
	python -c 'import align.assembly as A, igraph as ig; A.OverlapGraph(ig.read("$(TRUE_GRAPH).layout.gml")).diff_draw(A.OverlapGraph(ig.read("$(ASSEMBLED_GRAPH).layout.gml")), "$@")'

$(DB): align/libalign.so
	python -c 'import $(ASSEMBLY_TEST) as T; T.create_example("$@")'

loc:
	find align -type f -regex '.*\(\.py\|\.c\|\.h\)' | xargs wc -l

CAIRO=$(shell python -c 'import site, os.path; print filter(lambda e:os.path.isdir(e + "/cairo"), site.getsitepackages())[0] + "/cairo"')
env:
	virtualenv --system-site-packages $@
	. env/bin/activate && python setup.py install
	ln -s "$(CAIRO)" env/lib/python2.7/cairo

DOCS_EXCLUDE = align/tests
DOCS_OUT = _build
docs: align/libalign.so docs/README.rst
	sphinx-apidoc -f -e -o $@ align/ $(DOCS_EXCLUDE)
	cd $@ && sphinx-build -b html . $(DOCS_OUT)
	@echo "Find the docs at file://`readlink -f $@/$(DOCS_OUT)/index.html`"

docs/README.rst:
	pandoc -f markdown -t rst -o $@ -i README.md

.PHONY: clean tests *.pdf loc docs docs/README.rst
