# The Nyu Test

## A Bootstrap-based Test for Interspecific Phenotype Shifts in Phylogenetic Data

### Description

Nyu is a computational pipeline designed to analyze quantitative traits across phylogenetic trees. It uses bootstrap-based simulations under stochastic models of evolution (e.g., Brownian Motion or Ornstein-Uhlenbeck) to identify interspecific phenotype shifts and evaluate their statistical significance.

### Features
* Flexible Models: Supports Brownian Motion (BM) and Ornstein-Uhlenbeck (OUrandomRoot) models for evolutionary simulations.
* Z-score Standardization: Observed and simulated trait values are standardized for robust comparison.
* Outlier Detection: Identifies significant pairwise differences based on a 95% confidence interval (nyu95).
* Seamless Integration: Combines Python and R for efficient simulation and analysis.

### Files

**ntest.py**

A Python script that:
* Reads a trait database and phylogenetic tree.
* Prunes the tree to match the species in the dataset.
* Standardizes observed values and generates pairwise comparisons.
* Runs evolutionary simulations in R.
* Calculates pairwise and trait-wide nyu95 values.

**Usage:**

`python ntest.py -d <database_file> -t <tree_file> -v <trait_name> -n <number_of_simulations> -m`

**Arguments:**
```
-d, --database: Path to the dataset file (default: database/nhp.phenomic.dataset.tsv).
-t, --tree: Path to the phylogenetic tree file (default: database/phylogeny.nw).
-v, --trait: Trait to analyze (default: Body_mass).
-n, --nsim: Number of simulations (default: 1000).
-m, --model: Model for simulation, either BM (Brownian Motion) or OUrandomRoot (default: BM).
```

**Output Files**
* <trait_name>.single_values.csv: Observed trait values with z-scores.
* <trait_name>.simulations.csv: Simulated trait values.
* <trait_name>.simulations.zscore.csv: Z-scores of simulated values.
* <trait_name>.nyu.tab: Pairwise nyu values for all species combinations.
* <trait_name>.outliers.tab: Outliers based on nyu95 threshold.
* <trait_name>.nyu95.txt: Overall nyu95 score.

**Installation and Dependencies**
* Python (≥3.7)
* Required libraries: pandas, dendropy
* R (≥4.0)
* Required libraries: ape, geiger, phytools, phylolm
