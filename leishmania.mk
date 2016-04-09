READS = leishmania/reads.fa
GENOME = leishmania/genome.fa
REFERENCE = leishmania/reference.fa
ASSEMBLY_DB = genome.leishmania.db
ASSEMBLED_GRAPH = leishmania
TRUE_GRAPH = leishmania_true
MAPPINGS = leishmania/bwa.mappings.txt
ASSEMBLY_OPTS = DB=$(ASSEMBLY_DB) READS=$(READS) ASSEMBLED_GRAPH=$(ASSEMBLED_GRAPH) TRUE_GRAPH=$(TRUE_GRAPH) MAPPINGS=$(MAPPINGS)

$(GENOME):
	wget http://www.cs.mcgill.ca/~blanchem/561/pacbio/GCA_000227135.2_ASM22713v2_genomic.fna -O $@

$(READS):
	wget http://www.cs.mcgill.ca/~blanchem/561/pacbio/readsMappingToChr1.fa -O $@

# read the very first chromosome and write it to the reference file
$(REFERENCE): $(GENOME)
	python -c 'from Bio import SeqIO as I, SeqRecord as R; \
		rec = I.parse("$(GENOME)", "fasta").next(); \
		I.write([R.SeqRecord(rec.seq, id=rec.id)], "$@", "fasta")'

leishmania/bwa.mappings.txt: $(REFERENCE) $(READS)
	python -c 'from biseqt.mapping import BwaReadMapper as M, save_mappings; \
		save_mappings("$@", M(read_src="$(READS)", ref_src="$(REFERENCE)").mappings());'

LASTZ_PATH=/home/amir/lastz-distrib/bin/lastz
leishmania/lastz.mappings.txt: $(REFERENCE) $(READS)
	python -c 'from biseqt.mapping import LastzReadMapper as M, save_mappings; \
		save_mappings("$@", M(read_src="$(READS)", ref_src="$(REFERENCE)", lastz_path="$(LASTZ_PATH)").mappings());'

ASSEMBLY_TARGET = leishmania.gml
assembly:
	$(MAKE) -f assembly.mk $(ASSEMBLY_TARGET) $(ASSEMBLY_OPTS)

clean:
	rm -f $(REFERENCE)*
	rm -f $(ASSEMBLY_DB)
	make -f assembly.mk clean $(ASSEMBLY_OPTS)

.PHONY: clean
