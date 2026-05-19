FILES = sbdi-gtdb-sativa.r11rs232-1.1genome.addSpecies.fna.gz sbdi-gtdb-sativa.r11rs232-1.1genome.assignTaxonomy.fna.gz sbdi-gtdb-sativa.r11rs232-1.5genomes.addSpecies.fna.gz sbdi-gtdb-sativa.r11rs232-1.5genomes.assignTaxonomy.fna.gz sbdi-gtdb-sativa.r11rs232-1.20genomes.addSpecies.fna.gz sbdi-gtdb-sativa.r11rs232-1.20genomes.assignTaxonomy.fna.gz sbdi-gtdb-sativa.r11rs232-1.20genomes.fna.gz sbdi-gtdb-sativa.r11rs232-1.ar53.sprepr.alnfna sbdi-gtdb-sativa.r11rs232-1.ar53.sprepr.brlenopt.newick sbdi-gtdb-sativa.r11rs232-1.ar53.sprepr.taxonomy.tsv sbdi-gtdb-sativa.r11rs232-1.bac120.sprepr.alnfna sbdi-gtdb-sativa.r11rs232-1.bac120.sprepr.brlenopt.newick sbdi-gtdb-sativa.r11rs232-1.bac120.sprepr.taxonomy.tsv

clean:
	rm -f MANIFEST*.txt
	rm -f README.txt
	rm -f *.fna
	rm -f *.alnfna
	rm -f *.treefile

sbdi-gtdb-sativa.r11rs232-1.1genome.addSpecies.fna.gz: $(wildcard ../../r232-1/correct/*n1.addSpecies.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r11rs232-1.1genome.assignTaxonomy.fna.gz: $(wildcard ../../r232-1/correct/*n1.assignTaxonomy.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r11rs232-1.5genomes.addSpecies.fna.gz: $(wildcard ../../r232-1/correct/*n5.addSpecies.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r11rs232-1.5genomes.assignTaxonomy.fna.gz: $(wildcard ../../r232-1/correct/*n5.assignTaxonomy.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r11rs232-1.20genomes.addSpecies.fna.gz: $(wildcard ../../r232-1/correct/*n20.addSpecies.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r11rs232-1.20genomes.assignTaxonomy.fna.gz: $(wildcard ../../r232-1/correct/*n20.assignTaxonomy.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r11rs232-1.20genomes.fna.gz: $(wildcard ../../r232-1/correct/*n20.correct.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r11rs232-1.ar53.sprepr.alnfna: ../../r232-1/iqtree/arc-ssu-r232-sprep.alnfna
	ln -s $< $@

sbdi-gtdb-sativa.r11rs232-1.ar53.sprepr.brlenopt.newick: ../../r232-1/iqtree/arc-ssu-r232-sprep.brlenopt.treefile
	ln -s $< $@

sbdi-gtdb-sativa.r11rs232-1.ar53.sprepr.taxonomy.tsv: ../../r232-1/iqtree/arc-ssu-r232-sprep.taxonomy.tsv
	ln -s $< $@

sbdi-gtdb-sativa.r11rs232-1.bac120.sprepr.alnfna: ../../r232-1/iqtree/bac-ssu-r232-sprep.alnfna
	ln -s $< $@

sbdi-gtdb-sativa.r11rs232-1.bac120.sprepr.brlenopt.newick: ../../r232-1/iqtree/bac-ssu-r232-sprep.brlenopt.treefile
	ln -s $< $@

sbdi-gtdb-sativa.r11rs232-1.bac120.sprepr.taxonomy.tsv: ../../r232-1/iqtree/bac-ssu-r232-sprep.taxonomy.tsv
	ln -s $< $@

README.txt: ../../README.txt
	ln -s $< $@

MANIFEST.txt: README.txt $(FILES)
	ls -lLh $^ | sed 's/,/./g' | sed 's/.* \([0-9.]\+\)\([KM]\).* \(.*\)/\3 [\1 \2iB]/' > $@

sbdi-gtdb.r11rs232-1.figshare.tar.gz: MANIFEST.txt README.txt $(FILES)
	tar cvfzh $@ $^
