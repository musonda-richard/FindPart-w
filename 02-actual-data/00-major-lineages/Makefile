all:freq_by_day.csv S2_Fig.tif major_variants.csv full_variant_counts.csv count_variants.csv 

full_variant_counts.csv:
	cp ../../full_variant_counts.csv .

freq_by_day.csv:full_variant_counts.csv
	R -f freq_by_day.R

major_variants.csv S2_Fig.tif:freq_by_day.csv
	R -f select_major_variants_normalz.R

count_variants.csv:major_variants.csv full_variant_counts.csv
	R -f count_variants.R

clean:
	rm -f freq_by_day.csv
	rm -f major_variants.csv
	rm -f count_variants.csv
	rm -f full_variant_counts.csv 
	rm -f S2_Fig.tif
