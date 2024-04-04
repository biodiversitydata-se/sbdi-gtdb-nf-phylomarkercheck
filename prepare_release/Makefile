clean:
	rm MANIFEST*.txt
	rm README.txt
	rm *.fna.gz

MANIFEST_r214.txt: $(wildcard ../r214/correct/*n*.correct.fna.gz) $(wildcard ../r214/correct/*n*.assignTaxonomy.fna.gz) $(wildcard ../r214/correct/*n*.addSpecies.fna.gz)
	ln -s ../README.txt .
	ln -s $^ .
	ls -lLh *.fna.gz README.txt | sed 's/.*\([0-9]\+\)\([KM]\).* \(.*\)/\3 [\1 \2iB]/' > $@
	ln -s $@ MANIFEST.txt

sbdi-gtdb.r214.figshare.tar.gz: MANIFEST.txt README.txt $(wildcard *.fna.gz)
