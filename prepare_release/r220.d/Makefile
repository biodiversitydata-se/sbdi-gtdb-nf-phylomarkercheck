clean:
	rm -f MANIFEST*.txt
	rm -f README.txt
	rm -f *.fna
	rm -f *.alnfna
	rm -f *.treefile

sbdi-gtdb-sativa.r09rs220.1genome.addSpecies.fna.gz: $(wildcard ../../r220/correct/*n1.addSpecies.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r09rs220.1genome.assignTaxonomy.fna.gz: $(wildcard ../../r220/correct/*n1.assignTaxonomy.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r09rs220.5genomes.addSpecies.fna.gz: $(wildcard ../../r220/correct/*n5.addSpecies.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r09rs220.5genomes.assignTaxonomy.fna.gz: $(wildcard ../../r220/correct/*n5.assignTaxonomy.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r09rs220.20genomes.addSpecies.fna.gz: $(wildcard ../../r220/correct/*n20.addSpecies.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r09rs220.20genomes.assignTaxonomy.fna.gz: $(wildcard ../../r220/correct/*n20.assignTaxonomy.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r09rs220.20genomes.fna.gz: $(wildcard ../../r220/correct/*n20.correct.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r09rs220.ar53.sprepr.alnfna: ../../r220/iqtree/arc-ssu-r220-sprep.alnfna
	ln -s $< $@

sbdi-gtdb-sativa.r09rs220.ar53.sprepr.brlenopt.newick: ../../r220/iqtree/arc-ssu-r220-sprep.brlenopt.treefile
	ln -s $< $@

sbdi-gtdb-sativa.r09rs220.ar53.sprep.taxonomy.tsv: ../../r220/iqtree/arc-ssu-r220-sprep.taxonomy.tsv
	ln -s $< $@

sbdi-gtdb-sativa.r09rs220.bac120.sprepr.alnfna: ../../r220/iqtree/bac-ssu-r220-sprep.alnfna
	ln -s $< $@

sbdi-gtdb-sativa.r09rs220.bac120.sprepr.brlenopt.newick: ../../r220/iqtree/bac-ssu-r220-sprep.brlenopt.treefile
	ln -s $< $@

sbdi-gtdb-sativa.r09rs220.bac120.sprep.taxonomy.tsv: ../../r220/iqtree/bac-ssu-r220-sprep.taxonomy.tsv
	ln -s $< $@

README.txt: ../../README.txt
	ln -s $< $@

FILES = sbdi-gtdb-sativa.r09rs220.1genome.addSpecies.fna.gz sbdi-gtdb-sativa.r09rs220.1genome.assignTaxonomy.fna.gz sbdi-gtdb-sativa.r09rs220.5genomes.addSpecies.fna.gz sbdi-gtdb-sativa.r09rs220.5genomes.assignTaxonomy.fna.gz sbdi-gtdb-sativa.r09rs220.20genomes.addSpecies.fna.gz sbdi-gtdb-sativa.r09rs220.20genomes.assignTaxonomy.fna.gz sbdi-gtdb-sativa.r09rs220.20genomes.fna.gz sbdi-gtdb-sativa.r09rs220.ar53.sprepr.alnfna sbdi-gtdb-sativa.r09rs220.ar53.sprepr.brlenopt.newick sbdi-gtdb-sativa.r09rs220.ar53.sprep.taxonomy.tsv sbdi-gtdb-sativa.r09rs220.bac120.sprepr.alnfna sbdi-gtdb-sativa.r09rs220.bac120.sprepr.brlenopt.newick sbdi-gtdb-sativa.r09rs220.bac120.sprep.taxonomy.tsv

MANIFEST.txt: README.txt $(FILES)
	ls -lLh $^ | sed 's/,/./g' | sed 's/.* \([0-9.]\+\)\([KM]\).* \(.*\)/\3 [\1 \2iB]/' > $@

sbdi-gtdb.r09rs220.figshare.tar.gz: MANIFEST.txt README.txt $(FILES)
	tar cvfzh $@ $^
