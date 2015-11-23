# one of assembly or hp_assembly
MODE = assembly
ASSEMBLY_TEST = align.tests.$(MODE)
TRUE_GRAPH = true_overlap.$(MODE)
ASSEMBLED_GRAPH = overlap.$(MODE)
DB = genome.$(MODE).db
READS = reads.$(MODE).fa

clean:
	find . -regex .\*.fa  | grep -v '^./.git' | while read f; do rm -rf $$f; done
	find . -regex .\*.db  | grep -v '^./.git' | while read f; do rm -rf $$f; done
	find . -regex .\*.gml | grep -v '^./.git' | while read f; do rm -rf $$f; done
	find . -regex .\*.svg | grep -v '^./.git' | while read f; do rm -rf $$f; done

$(READS):
	python -c 'import $(ASSEMBLY_TEST) as T; T.create_example("$@", "$(READS)");'
$(DB): $(READS)
	python -c 'import $(ASSEMBLY_TEST) as T; T.create_db("$@", "$(READS)")'

$(ASSEMBLED_GRAPH).gml: $(DB) $(TRUE_GRAPH).gml
	python -c 'import $(ASSEMBLY_TEST) as T; T.overlap_by_seed_extension("$(DB)", "$@")'
	python -c 'import align.assembly as A, igraph as ig, sys; \
		g = A.OverlapGraph(ig.read("$(TRUE_GRAPH).gml")); \
		g.diff_text(A.OverlapGraph(ig.read("$(ASSEMBLED_GRAPH).gml")))'

$(ASSEMBLED_GRAPH).dag.gml: $(ASSEMBLED_GRAPH).gml
	python -c 'import align.assembly as A, igraph as ig; \
		G = A.OverlapGraph(ig.read("$(ASSEMBLED_GRAPH).gml")); \
		G.break_cycles(); G.save("$@")'

$(TRUE_GRAPH).gml: $(DB)
	python -c 'import $(ASSEMBLY_TEST) as A; \
		A.overlap_graph_by_known_order("$(DB)").save("$@")'

$(ASSEMBLED_GRAPH).pdf:
	python -c 'import align.assembly as A, igraph as ig; \
		A.OverlapGraph(ig.read("$(ASSEMBLED_GRAPH).gml")).draw("$@");'

$(TRUE_GRAPH).pdf: $(TRUE_GRAPH).gml
	python -c 'import align.assembly as A, igraph as ig; \
		A.OverlapGraph(ig.read("$(TRUE_GRAPH).gml")).draw("$@");'

# The true graph doesn't have any cycles, no need for a .dag.gml.
$(TRUE_GRAPH).layout.gml: $(TRUE_GRAPH).gml
	python -c 'import align.assembly as A, igraph as ig; \
		A.OverlapGraph(ig.read("$(TRUE_GRAPH).gml")).layout().save("$@")'

$(ASSEMBLED_GRAPH).layout.gml: $(ASSEMBLED_GRAPH).dag.gml
	python -c 'import align.assembly as A, igraph as ig; \
		A.OverlapGraph(ig.read("$(ASSEMBLED_GRAPH).dag.gml")).layout(full=True).save("$@")'

$(TRUE_GRAPH).layout.svg: $(TRUE_GRAPH).gml
	python -c 'import align.assembly as A, igraph as ig; \
		g = A.OverlapGraph(ig.read("$(TRUE_GRAPH).gml")); \
		g.draw("$@", highlight_paths=g.all_longest_paths());'

# When drawing the layout show all paths; the .gml file contains the longest path only.
$(ASSEMBLED_GRAPH).layout.svg: $(ASSEMBLED_GRAPH).layout.gml
	python -c 'import align.assembly as A, igraph as ig; \
		g = A.OverlapGraph(ig.read("$(ASSEMBLED_GRAPH).dag.gml")); \
		g.draw("$@", highlight_paths=g.all_longest_paths());'

layout.diff.$(MODE).svg: $(ASSEMBLED_GRAPH).layout.gml $(TRUE_GRAPH).layout.gml
	python -c 'import align.assembly as A, igraph as ig; \
		g = A.OverlapGraph(ig.read("$(TRUE_GRAPH).layout.gml")); \
		g.diff_draw(A.OverlapGraph(ig.read("$(ASSEMBLED_GRAPH).layout.gml")), "$@")'
