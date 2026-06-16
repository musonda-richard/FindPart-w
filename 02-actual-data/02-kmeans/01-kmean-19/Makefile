THREADS ?=16
FTOL_REL=1e-9
ALPHA=4.109
THETA=0.73
BASELINE=KP.3.1.1
NUM_TRIALS=19

all: copy_files kmeans-rgs-AIC.csv estimates.csv loglikelihood.csv consensus_matrix.tiff

copy_files:
	cp ../../00-major-lineages/count_variants.csv .
	cp ../../../RelRe.jl .
	cp ../../../lib_setpartition.jl .
	cp ../../../FindPart-w.jl .
	cp ../../../kmeans.jl .

estimates.csv: count_variants.csv
	julia --threads=$(THREADS) RelRe.jl -b $(BASELINE) -a $(ALPHA) -t $(THETA) -i $< -d 1 --ftol_rel $(FTOL_REL) -u

loglikelihood.csv: estimates.csv

kmeans-rgs-AIC.csv: count_variants.csv estimates.csv loglikelihood.csv
	julia --threads=$(THREADS) kmeans.jl count_variants.csv $(BASELINE) estimates.csv loglikelihood.csv $(FTOL_REL) $(ALPHA) $(THETA) $(NUM_TRIALS)

consensus_matrix.tiff: kmeans-rgs-AIC.csv count_variants.csv
	R -f plot-consensus.R

clean:
	find . -maxdepth 1 -name "*loglikelihood.csv" -delete
	find . -maxdepth 1 -name "*estimates.csv" -delete
	find . -maxdepth 1 -name "*blocks.csv" -delete
	rm -f kmeans-rgs-AIC.csv
	rm -f variant_order.csv
	rm -f nohup.out
	rm -f count_variants.csv
	rm -f RelRe.jl
	rm -f lib_setpartition.jl
	rm -f FindPart-w.jl
	rm -f *sus_k*
	rm -f consensus_matrix.tiff
