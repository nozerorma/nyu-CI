#!/bin/bash
# Run nyu test

phenotypes="malignant_prevalence neoplasia_prevalence"
run_modes="BM OU"

for p in ${phenotypes}; do
    for r in ${run_modes}; do

    mkdir -p test

    python ntest.py -d "/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/1.Cancer_data/Neoplasia_species360/cancer_traits-FINAL-DATASET-120424.csv" -t "/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/5.Phylogeny/science.abn7829_data_s4.nex.tree" -n 100000 -v $p -m $r -o test

    done
done
