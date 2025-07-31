THREADS=16
FTOL_REL = 1e-9
ALPHA = 4.109
THETA = 0.73
BASELINE = KP.3.1.1 
WIDTH=1

all:estimates.csv rgs-AIC.csv best_estimates.csv best_partition.csv original_RelRe_model.csv trajectory_best.pdf trajectory_original.pdf

estimates.csv:count_variants.csv 
	time julia --threads=$(THREADS) RelRe.jl -b $(BASELINE) -a $(ALPHA) -t $(THETA) -i count_variants.csv --ftol_rel $(FTOL_REL) -q -u -c 

original_RelRe_model.csv:estimates.csv 
	R -f original_RelRe_model.R

loglikelihood.csv:estimates.csv

frequencies.csv:estimates.csv

rgs-AIC.csv:count_variants.csv estimates.csv loglikelihood.csv
	time julia --threads=$(THREADS) FindPart-w.jl count_variants.csv $(BASELINE) estimates.csv loglikelihood.csv $(FTOL_REL) $(ALPHA) $(THETA) $(WIDTH)  

best-rgs.txt:rgs-AIC.csv
	head -n 2 rgs-AIC.csv | tail -n 1 | cut -d, -f 1 | tr '-' ' ' > best-rgs.txt

best_estimates.csv:count_variants.csv estimates.csv best-rgs.txt
	time julia --threads=$(THREADS) RelRe.jl -b $(BASELINE) -a $(ALPHA) -t $(THETA) -i count_variants.csv -d 1 --ftol_rel $(FTOL_REL) -v estimates.csv -o best -r $$(cat best-rgs.txt) -u -q -c 

best_frequencies.csv:best_estimates.csv

best_partition.csv:best_estimates.csv
	R -f best_partition.R

trajectory_best.pdf:best_frequencies.csv
	R -f plot_frequencies_best.R 
trajectory_original.pdf:frequencies.csv
	R -f plot_frequencies_original.R 
 
clean:
	find . -maxdepth 1 -name "*loglikelihood.csv" -delete
	find . -maxdepth 1 -name "*estimates.csv" -delete
	find . -maxdepth 1 -name "*blocks.csv" -delete
	find . -maxdepth 1 -name "*.log" -delete
	rm -f frequencies.csv
	rm -f rgs-AIC.csv
	rm -f rgs_no_chance.csv
	rm -f elapsed.csv
	rm -f best_*.csv
	rm -f original_RelRe_model.csv
	rm -f *.pdf 
	rm -f *.txt
	rm -f nohup.out

