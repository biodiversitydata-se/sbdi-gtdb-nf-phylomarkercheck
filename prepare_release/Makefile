clean:
	rm -f MANIFEST*.txt
	rm -f README.txt
	rm -f *.fna
	rm -f *.alnfna
	rm -f *.treefile

sbdi-gtdb-sativa.r08rs214.1genome.addSpecies.fna.gz: $(wildcard ../r214/correct/*n1.addSpecies.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r08rs214.1genome.assignTaxonomy.fna.gz: $(wildcard ../r214/correct/*n1.assignTaxonomy.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r08rs214.5genomes.addSpecies.fna.gz: $(wildcard ../r214/correct/*n5.addSpecies.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r08rs214.5genomes.assignTaxonomy.fna.gz: $(wildcard ../r214/correct/*n5.assignTaxonomy.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r08rs214.20genomes.addSpecies.fna.gz: $(wildcard ../r214/correct/*n20.addSpecies.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r08rs214.20genomes.assignTaxonomy.fna.gz: $(wildcard ../r214/correct/*n20.assignTaxonomy.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r08rs214.20genomes.fna.gz: $(wildcard ../r214/correct/*n20.correct.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r08rs214.ar53.sprepr.alnfna: ../r214/iqtree/arc-ssu-r214-sprep.alnfna
	ln -s $< $@

sbdi-gtdb-sativa.r08rs214.ar53.sprepr.brlenopt.newick: ../r214/iqtree/arc-ssu-r214-sprep.brlenopt.treefile
	ln -s $< $@

sbdi-gtdb-sativa.r08rs214.ar53.sprep.taxonomy.tsv: ../r214/iqtree/arc-ssu-r214-sprep.taxonomy.tsv
	ln -s $< $@

sbdi-gtdb-sativa.r08rs214.bac120.sprepr.alnfna: ../r214/iqtree/bac-ssu-r214-sprep.alnfna
	ln -s $< $@

sbdi-gtdb-sativa.r08rs214.bac120.sprepr.brlenopt.newick: ../r214/iqtree/bac-ssu-r214-sprep.brlenopt.treefile
	ln -s $< $@

sbdi-gtdb-sativa.r08rs214.bac120.sprep.taxonomy.tsv: ../r214/iqtree/bac-ssu-r214-sprep.taxonomy.tsv
	ln -s $< $@

README.txt: ../README.txt
	ln -s $< $@

FILES = sbdi-gtdb-sativa.r08rs214.1genome.addSpecies.fna.gz sbdi-gtdb-sativa.r08rs214.1genome.assignTaxonomy.fna.gz sbdi-gtdb-sativa.r08rs214.5genomes.addSpecies.fna.gz sbdi-gtdb-sativa.r08rs214.5genomes.assignTaxonomy.fna.gz sbdi-gtdb-sativa.r08rs214.20genomes.addSpecies.fna.gz sbdi-gtdb-sativa.r08rs214.20genomes.assignTaxonomy.fna.gz sbdi-gtdb-sativa.r08rs214.20genomes.fna.gz sbdi-gtdb-sativa.r08rs214.ar53.sprepr.alnfna sbdi-gtdb-sativa.r08rs214.ar53.sprepr.brlenopt.newick sbdi-gtdb-sativa.r08rs214.ar53.sprep.taxonomy.tsv sbdi-gtdb-sativa.r08rs214.bac120.sprepr.alnfna sbdi-gtdb-sativa.r08rs214.bac120.sprepr.brlenopt.newick sbdi-gtdb-sativa.r08rs214.bac120.sprep.taxonomy.tsv

MANIFEST.txt: README.txt $(FILES)
	ls -lLh $^ | sed 's/,/./g' | sed 's/.* \([0-9.]\+\)\([KM]\).* \(.*\)/\3 [\1 \2iB]/' > $@

sbdi-gtdb.r08rs214.figshare.tar.gz: MANIFEST.txt README.txt $(FILES)
	tar cvfzh $@ $^
