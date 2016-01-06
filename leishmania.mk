READS = leishmania/reads.fa
ANNOTATED_READS = leishmania/reads.annotated.fa
GENOME = leishmania/genome.fa
REFERENCE = leishmania/reference.fa
ASSEMBLY_DB = genome.leishmania.db
ASSEMBLED_GRAPH = leishmania
TRUE_GRAPH = leishmania_true
ASSEMBLY_OPTS = DB=$(ASSEMBLY_DB) READS=$(ANNOTATED_READS) ASSEMBLED_GRAPH=$(ASSEMBLED_GRAPH) TRUE_GRAPH=$(TRUE_GRAPH)

$(GENOME):
	wget http://www.cs.mcgill.ca/~blanchem/561/pacbio/GCA_000227135.2_ASM22713v2_genomic.fna -O $@

$(READS):
	wget http://www.cs.mcgill.ca/~blanchem/561/pacbio/readsMappingToChr1.fa -O $@

# read the very first chromosome and write it to the reference file
$(REFERENCE): $(GENOME)
	python -c 'from Bio import SeqIO as I, SeqRecord as R; \
		rec = I.parse("$(GENOME)", "fasta").next(); \
		I.write([R.SeqRecord(rec.seq, id=rec.id)], "$@", "fasta")'

# annotate all reads with their actual position in the genome
$(ANNOTATED_READS): $(REFERENCE) $(READS)
	bwa index $(REFERENCE) 2>&1 | sed 's/^/|  /' >&2
	READS=$(READS) DB=$(REFERENCE) python prepare.py $@

ASSEMBLY_TARGET = leishmania.gml
assembly:
	$(MAKE) -f assembly.mk $(ASSEMBLY_TARGET) $(ASSEMBLY_OPTS)

clean:
	rm -f $(REFERENCE)*
	rm -f $(ANNOTATED_READS)
	rm -f $(ASSEMBLY_DB)
	make -f assembly.mk clean $(ASSEMBLY_OPTS)

.PHONY: clean
