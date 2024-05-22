all: arc-ssu-r220.fna bac-ssu-r220.fna arc-ssu-r214.fna bac-ssu-r214.fna

ssu_all_r220.fna.gz:
	wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/genomic_files_all/$@

arc-ssu-r220.fna: ssu_all_r220.fna.gz
	zgrep -A 1 'd__Arc' $< | grep -v '^--' > $@

bac-ssu-r220.fna: ssu_all_r220.fna.gz
	zgrep -A 1 'd__Bac' $< | grep -v '^--' > $@

arc-ssu-r214.fna: ssu_all_r214.fna.gz
	zgrep -A 1 'd__Arc' $< | grep -v '^--' > $@

bac-ssu-r214.fna: ssu_all_r214.fna.gz
	zgrep -A 1 'd__Bac' $< | grep -v '^--' > $@
