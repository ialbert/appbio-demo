#
# Download genomic data from NCBI, one directory per type.
# Index FASTA (samtools faidx) and GFF (bgzip + tabix) for IGV.
#

BASE = https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/889/155/GCF_000889155.1_ViralProj51245
ACC = GCF_000889155.1_ViralProj51245
NAME = Bundibugyo

FA_GZ = fasta/${NAME}.fna.gz
FA = fasta/${NAME}.fa
GFF = gff/${NAME}.gff.gz

JAVA_HOME ?= /Library/Java/JavaVirtualMachines/jdk1.8.0_144.jdk/Contents/Home
JAVAC = ${JAVA_HOME}/bin/javac
JAVA = ${JAVA_HOME}/bin/java

all: ${FA}.fai ${GFF}.tbi

run: all src/browser.class
	${JAVA} -cp src browser ${FA} ${GFF}

src/browser.class: src/browser.java
	${JAVAC} src/browser.java

fasta gff:
	mkdir -p $@

${FA_GZ}: | fasta
	curl -L ${BASE}/${ACC}_genomic.fna.gz -o fasta/${ACC}_genomic.fna.gz
	mv fasta/${ACC}_genomic.fna.gz $@

${GFF}: | gff
	curl -L ${BASE}/${ACC}_genomic.gff.gz -o gff/${ACC}_genomic.gff.gz
	mv gff/${ACC}_genomic.gff.gz $@

${FA}: ${FA_GZ}
	gunzip -c $< > $@

${FA}.fai: ${FA}
	samtools faidx $<

${GFF}.tbi: ${GFF}
	gunzip -c ${GFF} > gff/${NAME}.gff
	(grep '^#' gff/${NAME}.gff; grep -v '^#' gff/${NAME}.gff | sort -k1,1 -k4,4n) | bgzip > ${GFF}.tmp
	mv ${GFF}.tmp ${GFF}
	rm -f gff/${NAME}.gff
	tabix -p gff ${GFF}

clean:
	rm -rf fasta gff src/*.class

.PHONY: all clean run
