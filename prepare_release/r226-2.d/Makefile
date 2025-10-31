clean:
	rm -f MANIFEST*.txt
	rm -f README.txt
	rm -f *.fna
	rm -f *.alnfna
	rm -f *.treefile

sbdi-gtdb-sativa.r10rs226-2.1genome.addSpecies.fna.gz: $(wildcard ../../r226-2.d/correct/*n1.addSpecies.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r10rs226-2.1genome.assignTaxonomy.fna.gz: $(wildcard ../../r226-2.d/correct/*n1.assignTaxonomy.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r10rs226-2.5genomes.addSpecies.fna.gz: $(wildcard ../../r226-2.d/correct/*n5.addSpecies.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r10rs226-2.5genomes.assignTaxonomy.fna.gz: $(wildcard ../../r226-2.d/correct/*n5.assignTaxonomy.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r10rs226-2.20genomes.addSpecies.fna.gz: $(wildcard ../../r226-2.d/correct/*n20.addSpecies.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r10rs226-2.20genomes.assignTaxonomy.fna.gz: $(wildcard ../../r226-2.d/correct/*n20.assignTaxonomy.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r10rs226-2.20genomes.fna.gz: $(wildcard ../../r226-2.d/correct/*n20.correct.fna.gz)
	gunzip -c $^ | gzip -c > $@

sbdi-gtdb-sativa.r10rs226-2.ar53.sprepr.alnfna: ../../r226-2.d/iqtree/arc-ssu-r226-sprep.alnfna
	ln -s $< $@

sbdi-gtdb-sativa.r10rs226-2.ar53.sprepr.brlenopt.newick: ../../r226-2.d/iqtree/arc-ssu-r226-sprep.brlenopt.treefile
	ln -s $< $@

sbdi-gtdb-sativa.r10rs226-2.ar53.sprepr.taxonomy.tsv: ../../r226-2.d/iqtree/arc-ssu-r226-sprep.taxonomy.tsv
	ln -s $< $@

sbdi-gtdb-sativa.r10rs226-2.bac120.sprepr.alnfna: ../../r226-2.d/iqtree/bac-ssu-r226-sprep.alnfna
	ln -s $< $@

sbdi-gtdb-sativa.r10rs226-2.bac120.sprepr.brlenopt.newick: ../../r226-2.d/iqtree/bac-ssu-r226-sprep.brlenopt.treefile
	ln -s $< $@

sbdi-gtdb-sativa.r10rs226-2.bac120.sprepr.taxonomy.tsv: ../../r226-2.d/iqtree/bac-ssu-r226-sprep.taxonomy.tsv
	ln -s $< $@

README.txt: ../../README.txt
	ln -s $< $@

FILES = sbdi-gtdb-sativa.r10rs226-2.1genome.addSpecies.fna.gz sbdi-gtdb-sativa.r10rs226-2.1genome.assignTaxonomy.fna.gz sbdi-gtdb-sativa.r10rs226-2.5genomes.addSpecies.fna.gz sbdi-gtdb-sativa.r10rs226-2.5genomes.assignTaxonomy.fna.gz sbdi-gtdb-sativa.r10rs226-2.20genomes.addSpecies.fna.gz sbdi-gtdb-sativa.r10rs226-2.20genomes.assignTaxonomy.fna.gz sbdi-gtdb-sativa.r10rs226-2.20genomes.fna.gz sbdi-gtdb-sativa.r10rs226-2.ar53.sprepr.alnfna sbdi-gtdb-sativa.r10rs226-2.ar53.sprepr.brlenopt.newick sbdi-gtdb-sativa.r10rs226-2.ar53.sprepr.taxonomy.tsv sbdi-gtdb-sativa.r10rs226-2.bac120.sprepr.alnfna sbdi-gtdb-sativa.r10rs226-2.bac120.sprepr.brlenopt.newick sbdi-gtdb-sativa.r10rs226-2.bac120.sprepr.taxonomy.tsv

MANIFEST.txt: README.txt $(FILES)
	ls -lLh $^ | sed 's/,/./g' | sed 's/.* \([0-9.]\+\)\([KM]\).* \(.*\)/\3 [\1 \2iB]/' > $@

sbdi-gtdb.r10rs226-2.figshare.tar.gz: MANIFEST.txt README.txt $(FILES)
	tar cvfzh $@ $^
