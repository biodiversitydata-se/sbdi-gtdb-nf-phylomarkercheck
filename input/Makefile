all: arc-ssu-r214.fna.gz bac-ssu-r214.fna.gz

arc-ssu-r214.fna.gz: ssu_all_r214.fna.gz
	zgrep -A 1 'd__Arc' $< | grep -v '^--' | gzip -c > $@

bac-ssu-r214.fna.gz: ssu_all_r214.fna.gz
	zgrep -A 1 'd__Bac' $< | grep -v '^--' | gzip -c > $@
