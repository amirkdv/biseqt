TRUE_GRAPH = true_overlap
ASSEMBLED_GRAPH = overlap
DB = genome.db
READS = reads.fa

clean:
	rm -f $(ASSEMBLED_GRAPH).gml $(ASSEMBLED_GRAPH).dag.gml $(ASSEMBLED_GRAPH).layout.gml
	rm -f $(DB) $(READS)
	rm -f $(TRUE_GRAPH).gml $(TRUE_GRAPH).layout.gml layout.diff.$(ASSEMBLED_GRAPH).svg

$(READS):
	python -c 'import align.tests.hp_assembly as T; T.create_example("$@", "$(READS)");'
$(DB): $(READS)
	python -c 'import align.tests.hp_assembly as T; T.create_db("$@", "$(READS)")'

WORDLEN = $(shell echo $$WORDLEN)
word_pvalues.$(WORDLEN).png:
	python -c 'from align import tuples, seq; \
		A = seq.Alphabet("ACGT"); \
		B = tuples.TuplesDB("$(DB)", alphabet=A); \
		I = tuples.Index(tuplesdb=B, wordlen=$(WORDLEN)); \
		I.plot_pvalues("$@");'

$(ASSEMBLED_GRAPH).gml:
	python -c 'import align.tests.hp_assembly as T; T.overlap_by_seed_extension("$(DB)", "$@")'
	$(MAKE) -f assembly.mk diff SUMMARY_ONLY=True

SUMMARY_ONLY = False
diff:
	python -c 'import align.assembly as A, igraph as ig, sys; \
		g = A.OverlapGraph(ig.read("$(TRUE_GRAPH).gml")); \
		g.diff_text(A.OverlapGraph(ig.read("$(ASSEMBLED_GRAPH).gml")), summary_only=$(SUMMARY_ONLY))'

FAS_METHOD = ip
$(ASSEMBLED_GRAPH).dag.gml: $(ASSEMBLED_GRAPH).gml
	python -c 'import align.assembly as A, igraph as ig; \
		G = A.OverlapGraph(ig.read("$(ASSEMBLED_GRAPH).gml")); \
		G.break_cycles(method="$(FAS_METHOD)"); G.save("$@")'

$(TRUE_GRAPH).gml:
	python -c 'import align.tests.hp_assembly as A; \
		A.overlap_graph_by_known_order("$(DB)").save("$@")'

$(ASSEMBLED_GRAPH).svg:
	python -c 'import align.assembly as A, igraph as ig; \
		A.OverlapGraph(ig.read("$(ASSEMBLED_GRAPH).gml")).draw("$@");'

$(TRUE_GRAPH).svg: $(TRUE_GRAPH).gml
	python -c 'import align.assembly as A, igraph as ig; \
		A.OverlapGraph(ig.read("$(TRUE_GRAPH).gml")).draw("$@");'

# The true graph doesn't have any cycles, no need for a .dag.gml.
$(TRUE_GRAPH).layout.gml: $(TRUE_GRAPH).gml
	python -c 'import align.assembly as A, igraph as ig; \
		A.OverlapGraph(ig.read("$(TRUE_GRAPH).gml")).layout(equal_weights=True).save("$@")'

$(ASSEMBLED_GRAPH).layout.gml: $(ASSEMBLED_GRAPH).dag.gml
	python -c 'import align.assembly as A, igraph as ig; \
		A.OverlapGraph(ig.read("$(ASSEMBLED_GRAPH).dag.gml")).layout(full=True).save("$@")'

$(TRUE_GRAPH).layout.svg: $(TRUE_GRAPH).gml
	python -c 'import align.assembly as A, igraph as ig; \
		g = A.OverlapGraph(ig.read("$(TRUE_GRAPH).gml")); \
		g.draw("$@", highlight_paths=g.all_longest_paths(equal_weights=True));'

# When drawing the layout show all paths; the .gml file contains the longest path only.
$(ASSEMBLED_GRAPH).layout.svg: $(ASSEMBLED_GRAPH).dag.gml
	python -c 'import align.assembly as A, igraph as ig; \
		g = A.OverlapGraph(ig.read("$(ASSEMBLED_GRAPH).dag.gml")); \
		g.draw("$@", highlight_paths=g.all_longest_paths());'

layout.diff.$(ASSEMBLED_GRAPH).svg: $(ASSEMBLED_GRAPH).layout.gml $(TRUE_GRAPH).layout.gml
	python -c 'import align.assembly as A, igraph as ig; \
		g = A.OverlapGraph(ig.read("$(TRUE_GRAPH).layout.gml")); \
		g.diff_draw(A.OverlapGraph(ig.read("$(ASSEMBLED_GRAPH).layout.gml")), "$@")'
