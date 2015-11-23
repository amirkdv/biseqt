READS = leishmania/reads.fa
ANNOTATED_READS = leishmania/reads.annotated.fa
GENOME = leishmania/genome.fa
REFERENCE = leishmania/reference.fa
BLAST_DB = leishmania/blast.db
ASSEMBLY_DB = genome.leishmania.db

$(GENOME):
	wget http://www.cs.mcgill.ca/~blanchem/561/pacbio/GCA_000227135.2_ASM22713v2_genomic.fna -O $@

$(READS):
	wget http://www.cs.mcgill.ca/~blanchem/561/pacbio/readsMappingToChr1.fa -O $@

# read the very first chromosome and write it to the reference file
$(REFERENCE): $(GENOME)
	python -c 'from Bio import SeqIO as I, SeqRecord as R; rec = I.parse("$(GENOME)", "fasta").next(); I.write([R.SeqRecord(rec.seq, id=rec.id)], "$@", "fasta")'

# annotate all reads with their actual position in the genome
$(ANNOTATED_READS): $(REFERENCE) $(READS)
	makeblastdb -dbtype nucl -in $(REFERENCE) -out $(BLAST_DB)
	READS=$(READS) DB=$(BLAST_DB) python prepare.py $@

$(ASSEMBLY_DB): $(ANNOTATED_READS)
	python -c 'import align.tests.assembly as T; T.create_db("$@", "$(ANNOTATED_READS)")'
