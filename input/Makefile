all: arc-ssu-r232.fna bac-ssu-r232.fna arc-ssu-r226.fna bac-ssu-r226.fna arc-ssu-r220.fna bac-ssu-r220.fna

ssu_all_r232.fna.gz:
	wget https://data.gtdb.ecogenomic.org/releases/release232/232.0/genomic_files_all/$@

arc-ssu-r232.fna: ssu_all_r232.fna.gz
	zgrep -A 1 'd__Arc' $< | grep -v '^--' > $@

bac-ssu-r232.fna: ssu_all_r232.fna.gz
	zgrep -A 1 'd__Bac' $< | grep -v '^--' > $@

ssu_all_r226.fna.gz:
	wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/genomic_files_all/$@

arc-ssu-r226.fna: ssu_all_r226.fna.gz
	zgrep -A 1 'd__Arc' $< | grep -v '^--' > $@

bac-ssu-r226.fna: ssu_all_r226.fna.gz
	zgrep -A 1 'd__Bac' $< | grep -v '^--' > $@

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
