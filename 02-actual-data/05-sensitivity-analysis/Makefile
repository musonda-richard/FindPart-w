THREADS=250
FTOL_REL=1e-9
ALPHA=4.109
THETA=0.73
BASELINE=KP.3.1.1
WIDTH=1

all: copy_files estimates.csv rgs-AIC.csv

copy_files:
	cp ../01-USA-data-analysis/count_variants.csv .
	cp ../01-USA-data-analysis/RelRe.jl .
	cp ../01-USA-data-analysis/lib_setpartition.jl .
	cp ../01-USA-data-analysis/FindPart-w.jl .

estimates.csv: count_variants.csv
	time julia --threads=$(THREADS) RelRe.jl -b $(BASELINE) -a $(ALPHA) -t $(THETA) -i count_variants.csv --ftol_rel $(FTOL_REL) -u

loglikelihood.csv: estimates.csv

rgs-AIC.csv: count_variants.csv estimates.csv loglikelihood.csv
	julia --threads=$(THREADS) FindPart-w.jl count_variants.csv $(BASELINE) estimates.csv loglikelihood.csv $(FTOL_REL) $(ALPHA) $(THETA) $(WIDTH)

clean:
	find . -maxdepth 1 -name "*loglikelihood.csv" -delete
	find . -maxdepth 1 -name "*estimates.csv" -delete
	find . -maxdepth 1 -name "*blocks.csv" -delete
	find . -maxdepth 1 -name "*.log" -delete
	rm -f rgs-AIC.csv
	rm -f rgs_no_chance.csv
	rm -f elapsed.csv
	rm -f nohup.out
	rm -f count_variants.csv
	rm -f RelRe.jl
	rm -f lib_setpartition.jl
	rm -f FindPart-w.jl
