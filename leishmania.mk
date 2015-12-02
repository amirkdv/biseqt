READS = leishmania/reads.fa
ANNOTATED_READS = leishmania/reads.annotated.fa
GENOME = leishmania/genome.fa
REFERENCE = leishmania/reference.fa
BLAST_DB = leishmania/blast.db
ASSEMBLY_DB = genome.leishmania.db
MODE = hp_assembly
ASSEMBLED_GRAPH = leishmania
TRUE_GRAPH = leishmania_true
ASSEMBLY_OPTS = MODE=$(MODE) DB=$(ASSEMBLY_DB) READS=$(ANNOTATED_READS) ASSEMBLED_GRAPH=$(ASSEMBLED_GRAPH) TRUE_GRAPH=$(TRUE_GRAPH)
NUM_READS = -1

$(GENOME):
	wget http://www.cs.mcgill.ca/~blanchem/561/pacbio/GCA_000227135.2_ASM22713v2_genomic.fna -O $@

$(READS):
	wget http://www.cs.mcgill.ca/~blanchem/561/pacbio/readsMappingToChr1.fa -O $@

# read the very first chromosome and write it to the reference file
$(REFERENCE): $(GENOME)
	python -c 'from Bio import SeqIO as I, SeqRecord as R; rec = I.parse("$(GENOME)", "fasta").next(); I.write([R.SeqRecord(rec.seq, id=rec.id)], "$@",
	"fasta")'

$(BLAST_DB): $(REFERENCE)
	makeblastdb -dbtype nucl -in $(REFERENCE) -out $@ 2>&1 | sed 's/^/|  /'

# annotate all reads with their actual position in the genome
$(ANNOTATED_READS): $(READS) $(BLAST_DB)
	READS=$(READS) DB=$(BLAST_DB) NUM_READS=$(NUM_READS) python prepare.py $@

$(ASSEMBLY_DB): $(ANNOTATED_READS)
	python -c 'import align.tests.$(MODE) as A; A.create_db("$@", "$(ANNOTATED_READS)")'

ASSEMBLY_TARGET = leishmania.gml
assembly:
	@echo $(ASSEMBLY_TARGET)
	$(MAKE) -f assembly.mk $(ASSEMBLY_TARGET) $(ASSEMBLY_OPTS)

clean:
	rm -f $(BLAST_DB)*
	rm -f $(ANNOTATED_READS)
	rm -f $(ASSEMBLY_DB)
	make -f assembly.mk clean $(ASSEMBLY_OPTS)

.PHONY: clean
