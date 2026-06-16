THREADS ?= 16
FTOL_REL=1e-9
ALPHA=4.109
THETA=0.73
BASELINE=KP.3.1.1

all: copy_files hierarchical-rgs-AIC.csv estimates.csv loglikelihood.csv

copy_files:
	cp ../../RelRe.jl .
	cp ../../lib_setpartition.jl .
	cp ../../FindPart-w.jl .
	cp ../00-major-lineages/count_variants.csv .

estimates.csv: count_variants.csv
	julia --threads=$(THREADS) RelRe.jl -b $(BASELINE) -a $(ALPHA) -t $(THETA) -i $< -d 1 --ftol_rel $(FTOL_REL) -u

loglikelihood.csv: estimates.csv

hierarchical-rgs-AIC.csv: count_variants.csv estimates.csv loglikelihood.csv
	julia --threads=$(THREADS) hclust.jl count_variants.csv $(BASELINE) estimates.csv $(FTOL_REL) $(ALPHA) $(THETA)

clean:
	find . -maxdepth 1 -name "*loglikelihood.csv" -delete
	find . -maxdepth 1 -name "*estimates.csv" -delete
	find . -maxdepth 1 -name "*blocks.csv" -delete
	rm -f hierarchical-rgs-AIC.csv
	rm -f variant_order.csv
	rm -f nohup.out
	rm -f count_variants.csv
	rm -f RelRe.jl
	rm -f lib_setpartition.jl
	rm -f FindPart-w.jl
