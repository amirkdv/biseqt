GCC = gcc -shared -Wall -fPIC -std=c99 -g
TRUE_GRAPH = true_overlap
ASSEMBLED_GRAPH = overlap
DB = genome.db

# To get coredumps:
# 	$ ulimit -c unlimited
# To debug coredump:
# 	$ gdb align/libalign.so core
align/libalign.so: align/libalign.c align/libalign.h
	$(GCC) align/libalign.c -o $@

clean:
	find . -regex .\*.pyc | while read f; do rm -f $$f; done
	find . -regex .\*.swp | while read f; do rm -f $$f; done
	rm -f align/libalign.so core
	rm -f genome.fa reads.fa $(DB)
	rm -f *.gml *.pdf

tests: align/libalign.so $(DB) $(TRUE_GRAPH).gml $(ASSEMBLED_GRAPH).gml $(ASSEMBLED_GRAPH).layout.diff.pdf
	python -m align.tests.homopolymeric
	python -m align.tests.align

$(ASSEMBLED_GRAPH).gml:
	python -c 'import align.tests.assembly as T; T.overlap_by_tuple_extension("$(DB)", "$@")'
	python -c 'import align.tests.assembly as T; T.compare_results("$(TRUE_GRAPH).gml", "$(ASSEMBLED_GRAPH).gml")'

$(ASSEMBLED_GRAPH).pdf:
	python -c 'import align.assembly as A, networkx as nx; A.draw_digraph(nx.read_gml("$(ASSEMBLED_GRAPH).gml"), "$@");'

$(ASSEMBLED_GRAPH).layout.gml:
	python -c 'import align.assembly as A, networkx as nx; A.save_graph(A.layout_graph(nx.read_gml("$(ASSEMBLED_GRAPH).gml")), "$@")'

$(ASSEMBLED_GRAPH).layout.pdf:
	python -c 'import align.assembly as A, networkx as nx; A.draw_digraph(nx.read_gml("$(ASSEMBLED_GRAPH).gml"), "$@", longest_path=True);'

$(ASSEMBLED_GRAPH).layout.diff.pdf: $(ASSEMBLED_GRAPH).layout.gml $(TRUE_GRAPH).layout.gml
	python -c 'import align.assembly as A, networkx as nx; A.diff_graph(nx.read_gml("$(TRUE_GRAPH).layout.gml"), nx.read_gml("$(ASSEMBLED_GRAPH).layout.gml"), "$@")'

$(TRUE_GRAPH).gml:
	python -c 'import align.tests.assembly as T; T.overlap_by_known_order("$(DB)", "$@")'

$(TRUE_GRAPH).pdf:
	python -c 'import align.assembly as A, networkx as nx; A.draw_digraph(nx.read_gml("$(TRUE_GRAPH).gml"), "$@");'

$(TRUE_GRAPH).layout.gml:
	python -c 'import align.assembly as A, networkx as nx; A.save_graph(A.layout_graph(nx.read_gml("$(TRUE_GRAPH).gml")), "$@")'

$(TRUE_GRAPH).layout.pdf:
	python -c 'import align.assembly as A, networkx as nx; A.draw_digraph(nx.read_gml("$(TRUE_GRAPH).gml"), "$@", longest_path=True);'

$(DB): align/libalign.so
	python -c 'import align.tests.assembly as T; T.create_example("$@")'

loc:
	find . -type f -regex '.*\(\.py\|\.c\|\.h\)' | xargs wc -l

docs:
	sphinx-apidoc -F -o $@ align/
	cd $@ && sphinx-build -b html . _build
	@echo "Find the docs at file://`readlink -f $@/_build/index.html`"

.PHONY: clean tests *.gml *.pdf loc docs
