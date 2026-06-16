
all: lineage_tree.tiff pango_distance_matrix.csv

lineage_tree.tiff pango_distance_matrix.csv:../01-USA-data-analysis/rgs-AIC.csv 
	R -f lineage_tree.R

clean:
	rm -f lineage_tree.tiff
	rm -f pango_distance_matrix.csv
