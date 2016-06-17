READS = leishmania/reads.fa
GENOME = leishmania/genome.fa
REFERENCE = leishmania/reference.fa
ASSEMBLY_DB = leishmania/genome.db
ASSEMBLED_GRAPH = leishmania/graph
TRUE_GRAPH = leishmania/true
MAPPINGS = leishmania/blasr.mappings.txt
ASSEMBLY_OPTS = DB=$(ASSEMBLY_DB) READS=$(READS) ASSEMBLED_GRAPH=$(ASSEMBLED_GRAPH) TRUE_GRAPH=$(TRUE_GRAPH) MAPPINGS=$(MAPPINGS)

DOCKER_OPTS=-it -v "$$PWD:/mnt" -v /tmp:/tmp -w /mnt

$(GENOME):
	wget http://www.cs.mcgill.ca/~blanchem/561/pacbio/GCA_000227135.2_ASM22713v2_genomic.fna -O $@

$(READS):
	wget http://www.cs.mcgill.ca/~blanchem/561/pacbio/readsMappingToChr1.fa -O $@

# read the very first chromosome and write it to the reference file
$(REFERENCE): $(GENOME)
	python -c 'from Bio import SeqIO as I, SeqRecord as R; \
		rec = I.parse("$(GENOME)", "fasta").next(); \
		I.write([R.SeqRecord(rec.seq, id=rec.id)], "$@", "fasta")'

bwa:
	cat Dockerfile.bwa | docker build -t $(BWA_IMG) -

BWA_SAM=leishmania/bwa.mappings.sam
BWA_IMG=bwa
BWA_CMD=bwa index $(REFERENCE) && bwa mem -x pacbio $(REFERENCE) $(READS) > $(BWA_SAM)
BWA=docker run $(DOCKER_OPTS) $(BWA_IMG) sh -c "$(BWA_CMD)"
leishmania/bwa.mappings.txt: $(REFERENCE) $(READS)
	$(BWA)
	python -c 'from biseqt.mapping import parse_mappings, save_mappings; \
		save_mappings("$@", parse_mappings("$(REFERENCE)", "$(READS)", "$(BWA_SAM)"));'

#FIXME document installation steps, this works nowhere else
LASTZ_PATH=/home/amir/lastz-distrib/bin/lastz
LASTZ_OPTS=--ambiguous=n --format=softsam --chain # --gap=31 --match=12
LASTZ=$(LASTZ_PATH) $(LASTZ_OPTS) $(REFERENCE) $(READS)
LASTZ_SAM=leishmania/lastz.mappings.sam
leishmania/lastz.mappings.txt: $(REFERENCE) $(READS)
	$(LASTZ) > $(LASTZ_SAM)
	python -c 'from biseqt.mapping import parse_mappings, save_mappings; \
		save_mappings("$@", parse_mappings("$(REFERENCE)", "$(READS)", "$(LASTZ_SAM)"));'

# name of docker image
blasr:
	cat Dockerfile.blasr | docker build -t $(BLASR_IMG) -

BLASR_IMG=blasr
BLASR_SAM=leishmania/blasr.mappings.sam
BLASR=docker run $(DOCKER_OPTS) $(BLASR_IMG) /var/blasr/blasr -bestn 1 $(READS) $(REFERENCE) -sam -out $(BLASR_SAM)
leishmania/blasr.mappings.txt: $(REFERENCE) $(READS)
	$(BLASR)
	python -c 'from biseqt.mapping import parse_mappings, save_mappings; \
		save_mappings("$@", parse_mappings("$(REFERENCE)", "$(READS)", "$(BLASR_SAM)", blasr=1));'


MAPPING_DB = leishmania/genome.mapping.db
$(MAPPING_DB):
	python -c 'import biseqt.tests.assembly as T; T.create_mapping_db("$@", "$(READS)", "$(REFERENCE)")'

MAPPING_TO_REFS=leishmania/our.mappings.txt
$(MAPPING_TO_REFS):
	python -c 'import biseqt.tests.assembly as T; T.map_to_refs("$@", "$(MAPPING_DB)")'

ASSEMBLY_TARGET = leishmania.gml
assembly:
	$(MAKE) -f assembly.mk $(ASSEMBLY_TARGET) $(ASSEMBLY_OPTS)

clean:
	rm -f $(REFERENCE)*
	rm -f $(ASSEMBLY_DB)
	make -f assembly.mk clean $(ASSEMBLY_OPTS)

.PHONY: clean blasr bwa
