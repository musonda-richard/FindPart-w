 FindPart-w: A Program to Identify SARS-CoV-2 Lineages Sharing the Same Relative Reproduction Numbers

## Authors

* Kimihito Ito
* Richard Musonda

---

## Overview

FindPart-w is a program for identifying groups of SARS-CoV-2 lineages that share the same relative effective reproduction number. This repository contains scripts for analyzing both simulated hypothetical datasets and real-world SARS-CoV-2 sequence count data from the United States.

The analyses include:

* Estimation of relative effective reproduction numbers using the original RelRe model
* Identification of lineage partitions using the FindPart-w algorithm
* Comparison with k-means clustering
* Comparison with centroid-linkage hierarchical clustering
* Sensitivity analysis for the assumed generation-time distribution
* Visualization of lineage relationships and model-comparison results

---

## Repository Structure

```text
01-Hypothetical-Dataset/
    Scripts for simulating and analyzing hypothetical SARS-CoV-2 lineage datasets

02-Actual-Dataset/
    00-major-lineages/
        Selection of major lineages and creation of count_variants.csv

    01-USA-data-analysis/
        Main real-world data analysis using the original RelRe model and FindPart-w

    02-kmeans/
        k-means clustering analyses

    03-hierarchical/
        Centroid-linkage hierarchical clustering analysis

    04-lineage_tree/
        Lineage tree visualization

    05-sensitivity-analysis/
        Sensitivity analysis for the generation-time distribution

    06-furtherplots/
        Additional plots and model-comparison summaries
```

---

## Required Packages

### Julia Packages

The following Julia packages are required to run the program:

* CSV
* DataFrames
* Distributions
* NLopt
* ArgParse
* Combinatorics
* SHA
* Dates
* Statistics

Install the required Julia packages using:

```bash
julia install.packages.jl
```

### R Packages

The following R packages are required to run the R scripts:

* tidyverse
* igraph
* ggraph
* ggrepel
* grid

Install the required R packages using:

```bash
R -f install.packages.R
```

Note: The `grid` package is normally included with base R.

---

## Reproducibility

Each analysis directory contains a `Makefile` or script that reproduces the corresponding analysis. The main commands are shown in the sections below.

Users can modify the number of threads using the `THREADS` variable in the Makefiles depending on the capacity of their computer system.

---

# Section 1: Simulating and Analyzing Hypothetical SARS-CoV-2 Datasets

## Purpose

This analysis uses simulated SARS-CoV-2 lineage replacement data for two purposes:

1. To evaluate the accuracy of the FindPart-w algorithm in identifying lineages that share the same relative reproduction number.
2. To evaluate the computational time required by the FindPart-w algorithm.

---

## Step 1: Change to the directory

```bash
cd 01-Hypothetical-Dataset
```

---

## Step 2: Run the full pipeline

```bash
make
```

This command executes the complete pipeline described below through the Makefile.

---

## 1. Generate true parameters

This step creates `parameters.csv`, which contains:

* Relative reproduction numbers (`k`)
* Initial frequencies (`qt`)

Command:

```bash
julia param_gen_rank.jl NUM_LIN RANK
```

Customizable parameters:

* `NUM_LIN`: number of lineages
* `RANK`: rank used in the simulation

---

## 2. Simulate daily lineage counts

This step simulates daily lineage counts that mimic real-world SARS-CoV-2 lineage replacement over time. The output file is `count_variants.csv`, which contains daily lineage counts.

Command:

```bash
julia RelRe_simul.jl parameters.csv NUM_SAMPLES
```

Customizable parameter:

* `NUM_SAMPLES`: total number of lineages sampled each day

---

## 3. Estimate parameters using the original RelRe model

This step creates `estimates.csv`, which contains estimates of:

* Relative effective reproduction numbers (`k`)
* Initial frequencies (`qt`)

Command:

```bash
julia RelRe.jl BASELINE count_variants.csv
```

Customizable parameter:

* `BASELINE`: baseline lineage used as the reference lineage

---

## 4. Identify the best partition using FindPart-w

This step uses the FindPart-w algorithm to find the best partition in which some lineages are constrained to share the same relative effective reproduction number.

Command:

```bash
julia FindPart-w.jl
```

Customizable parameter:

* `WIDTH`: limits the number of partitions evaluated at each rank, where the rank is less than the number of lineages

---

## 5. Final output

The final output of FindPart-w is:

```text
rgs-AIC.csv
```

This file contains the best partition found by FindPart-w, in which some lineages share the same relative effective reproduction number.

---

# Section 2: Analyzing Real-world SARS-CoV-2 Data from the United States

## Purpose

This analysis uses actual SARS-CoV-2 lineage sequence count data from the United States for two purposes:

1. To use FindPart-w to identify the best partition in which some lineages are constrained to share the same relative reproduction numbers.
2. To evaluate whether the partition found by FindPart-w achieves a better AIC than the partition found using the original RelRe model, k-means clustering, and centroid-linkage hierarchical clustering.

---

## Step 1: Selection of major lineages

Change to the following directory:

```bash
cd 02-Actual-Dataset/00-major-lineages/
```

Run the full pipeline:

```bash
make
```

### Purpose of this step

This step selects the major SARS-CoV-2 lineages used in the real-world data analysis and creates the lineage count file used by the downstream analyses.

The pipeline first copies `full_variant_counts.csv` into this directory. Then, `freq_by_day.R` creates `freq_by_day.csv`, which contains daily lineage frequencies. Next, `select_major_variants_normalz.R` calculates the mean frequency of each lineage across all days, ranks the lineages in decreasing order of mean frequency, and calculates the cumulative mean frequency.

Lineages were selected as major lineages if their cumulative mean frequency was within the 90% threshold. Finally, `count_variants.R` creates `count_variants.csv`, which contains the daily sequence counts for the selected major lineages. This file is used as the input file for the downstream real-world analyses.

### Main output files

* `freq_by_day.csv`: contains daily lineage frequencies.
* `major_variants.csv`: contains the selected major lineages, their mean frequencies, and their cumulative mean frequencies.
* `S2_Fig.tif`: shows the cumulative percentage of sequence counts against the number of ranked lineages.
* `count_variants.csv`: contains daily sequence counts for the selected major lineages and is used by the downstream analyses.

---

## Step 2: USA data analysis using the original RelRe model and FindPart-w

Change to the following directory:

```bash
cd 02-Actual-Dataset/01-USA-data-analysis/
```

Run the full pipeline:

```bash
make
```

This command executes the complete USA data analysis pipeline defined in the Makefile.

### Purpose of this step

This step estimates the parameters of the original RelRe model and then applies the FindPart-w algorithm to identify the best partition in which some lineages share the same relative effective reproduction number.

### Main output files

* `loglikelihood.csv`: contains the log-likelihood and AIC of the original RelRe model.
* `estimates.csv`: contains the estimated relative effective reproduction numbers (`k`) and initial frequencies (`qt`).
* `rgs-AIC.csv`: contains the AIC values of partitions evaluated by FindPart-w. The first row corresponds to the partition with the minimum AIC.
* `*_blocks.csv`: contains the lineage groups identified by FindPart-w.
* `*_estimates.csv`: contains the parameter estimates for the selected FindPart-w partition.

---

## Step 3: k-means clustering analysis

The k-means clustering analysis is performed using two settings: 19 trials and 100 trials.

---

### Step 3.1: k-means clustering with 19 trials

Change to the following directory:

```bash
cd 02-Actual-Dataset/02-kmeans/01-kmean-19/
```

Run the full pipeline:

```bash
make
```

### Purpose of this step

This step compares the partition found by FindPart-w with partitions obtained using k-means clustering. In this analysis, k-means clustering is performed with 19 trials.

### Main output files

* `kmeans-rgs-AIC.csv`: contains the AIC values of partitions obtained from k-means clustering. The first row corresponds to the partition with the minimum AIC.
* `consensus_k*.csv`: contains the consensus matrix for the selected number of clusters.
* `*_blocks.csv`: contains the lineage groups identified by k-means clustering.
* `*_estimates.csv`: contains the parameter estimates for the selected k-means partition.

---

### Step 3.2: k-means clustering with 100 trials

Change to the following directory:

```bash
cd 02-Actual-Dataset/02-kmeans/02-kmean-100/
```

Run the full pipeline:

```bash
make
```

### Purpose of this step

This step repeats the k-means clustering analysis using 100 trials to evaluate the stability of the clustering result.

### Main output files

* `kmeans-rgs-AIC.csv`: contains the AIC values of partitions obtained from k-means clustering. The first row corresponds to the partition with the minimum AIC.
* `consensus_k*.csv`: contains the consensus matrix for the selected number of clusters.
* `*_blocks.csv`: contains the lineage groups identified by k-means clustering.
* `*_estimates.csv`: contains the parameter estimates for the selected k-means partition.

---

## Step 4: Centroid-linkage hierarchical clustering analysis

Change to the following directory:

```bash
cd 02-Actual-Dataset/03-hierarchical/
```

Run the full pipeline:

```bash
make
```

This command executes the centroid-linkage hierarchical clustering analysis defined in the Makefile.

### Purpose of this step

This step compares the partition found by FindPart-w with partitions obtained using centroid-linkage hierarchical clustering. The hierarchical clustering procedure starts from individual lineages and progressively merges clusters. At each hierarchical level, the resulting partition is evaluated using AIC.

### Main output files

* `hierarchical-rgs-AIC.csv`: contains the AIC values of partitions obtained from hierarchical clustering. The first row corresponds to the partition with the minimum AIC.
* `*_blocks.csv`: contains the lineage groups identified by hierarchical clustering.
* `*_estimates.csv`: contains the parameter estimates for the selected hierarchical clustering partition.

---

## Step 5: Lineage tree visualization

Change to the following directory:

```bash
cd 02-Actual-Dataset/04-lineage_tree/
```

Run the full pipeline:

```bash
make
```

This command generates the lineage tree visualization.

### Purpose of this step

This step visualizes the relationships among the analyzed Pango lineages and displays the relative effective reproduction number estimates on the lineage tree.

### Main output files

* Lineage tree figure files showing the relationships among lineages.
* Output files used to visualize lineage relationships, recombinant lineages, and relative effective reproduction number estimates.

---

## Step 6: Sensitivity analysis

Change to the following directory:

```bash
cd 02-Actual-Dataset/05-sensitivity-analysis/
```

Run the full sensitivity analysis:

```bash
bash run.all.sh
```

### Purpose of this step

This step evaluates whether the best partition found by FindPart-w is robust to changes in the assumed generation-time distribution. The sensitivity analysis varies the mean generation time and variance, then reruns the FindPart-w analysis for each parameter setting.

The script tests the following mean generation times:

* 2.0 days
* 3.0 days
* 4.0 days
* 5.0 days
* 6.0 days

The script also tests the following variance values:

* 1.0
* 2.19
* 3.0

For each combination of mean generation time and variance, the script calculates the corresponding gamma-distribution parameters:

* `ALPHA = mean^2 / variance`
* `THETA = variance / mean`

The script then runs the FindPart-w analysis using these values and saves the resulting AIC output file for each parameter setting.

### Main output files

* `rgs-AIC_mu*_var*.csv`: contains the best partition and AIC values obtained under each mean-generation-time and variance setting.
* `combined_results.csv`: summarizes the results from all sensitivity-analysis runs, including the mean generation time, variance, `ALPHA`, `THETA`, partition, log-likelihood, number of parameters, and AIC.
* Output files generated during each FindPart-w run, including temporary estimate, log-likelihood, and block files. These temporary files are removed by the clean step after the analysis is completed.

---

## Step 7: Additional plots and model comparison

Change to the following directory:

```bash
cd 02-Actual-Dataset/06-furtherplots/
```

Run the full pipeline:

```bash
make
```

This command generates additional plots for comparing the original RelRe model, FindPart-w, k-means clustering, and centroid-linkage hierarchical clustering.

### Purpose of this step

This step summarizes and visualizes the results from the real-world dataset analysis. It includes comparisons of AIC values, the number of maximum likelihood estimations, computational time.

### Main output file
* Plots comparing the number of maximum likelihood estimations and computational time.

---

# Summary of the Real-world Data Analysis Workflow

The complete workflow is:

```bash
cd 02-Actual-Dataset/00-major-lineages/
make

cd ../01-USA-data-analysis/
make

cd ../02-kmeans/01-kmean-19/
make

cd ../02-kmeans/02-kmean-100/
make

cd ../../03-hierarchical/
make

cd ../04-lineage_tree/
make

cd ../05-sensitivity-analysis/
bash run.all.sh

cd ../06-furtherplots/
make
```

The key final outputs from this section are the selected major lineages, `count_variants.csv`, the AIC comparison results, the best partition found by FindPart-w, the k-means and centroid-linkage hierarchical clustering comparison results, the sensitivity-analysis results, and the figures used to summarize the real-world SARS-CoV-2 lineage analysis.

---

## Notes

* The number of threads can be modified through the `THREADS` variable in the Makefiles.
* The baseline lineage and optimization parameters can be modified in the corresponding Makefiles.
* Temporary files may be removed by running `make clean` in directories where a clean rule is provided.
