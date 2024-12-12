#!/usr/bin/env python3

# __   __
# \ \ / /
#  \ V / 
#   \_/  

# A Test for Quantitative Traits Over a Phylogeny
# Authors: Fabio Barteri (Universitat Pompeu Fabra - Barcelona)
#          Omar Cornejo (UCSC - Santa Cruz, CA)
# Date: December 2024

### IMPORTS
'''
USAGE:
python ntest.py $
'''


import pandas as pd
import dendropy
from itertools import combinations as cb
import os
import argparse

### INPUTS

# Configure the argparse
parser = argparse.ArgumentParser(description="A Test for Quantitative Traits Over a Phylogeny")

# Add the arguments
parser.add_argument(
    "-d", "--database", 
    type=str, 
    default="database/nhp.phenomic.dataset.tsv", 
    help="Path to the dataset file (default: database/nhp.phenomic.dataset.tsv)"
)
parser.add_argument(
    "-t", "--tree", 
    type=str, 
    default="database/phylogeny.nw", 
    help="Path to the phylogenetic tree file (default: database/phylogeny.nw)"
)
parser.add_argument(
    "-v", "--trait", 
    type=str, 
    default="Body_mass", 
    help="Trait to analyze (default: Body_mass)"
)
parser.add_argument(
    "-n", "--nsim", 
    type=int, 
    default=1000, 
    help="Number of simulations (default: 1000)"
)

parser.add_argument(
    "-m", "--model", 
    type=str, 
    default="BM", 
    help="Stochastic model for the simulations (BM = Brownian Motion, OU = Osteirn-Uhlenbeck)"
)

parser.add_argument(
    "-o", "--outfolder", 
    type=str, 
    default="results", 
    help="The output folder (default: test)"
)

# Parse the arguments
args = parser.parse_args()

# Assign input values
database_path = args.database
tree_path = args.tree
trait = args.trait
nsim = args.nsim
model = args.model

model = model.replace("OU", "OUrandomRoot")

mainfolder = "test"

outfolder = mainfolder + "/" + trait

try:
    os.mkdir(outfolder)
except:
    pass
# The output tag

output_tag = outfolder + "/" + trait + "." + model

### FUNCTIONS
#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@

#@@ FUNCTION calculate_abs_diff_row()
def calculate_abs_diff(row, trait_db, trait):
    '''
    Calculates the absolute difference between two species and it returns it to
    a new row in the trait dataframe. 
    '''
    value1 = trait_db.loc[row["Species1"], trait]  # Species1
    value2 = trait_db.loc[row["Species2"], trait]  # Species2

    return abs(value1 - value2)

#@@ FUNCTION set_distance()
def set_distance(row, distance_dictionary):
    '''
    Calculates the distance between two species. 
    '''
    s1 = row["Species1"] # Species1
    s2 = row["Species2"]  # Species2

    try:
        distance = distance_dictionary[s1 + "," + s2]
    except:
        distance = distance_dictionary[s2 + "," + s1]

    return distance

#@@ FUNCTION calc_nyu()
def calc_nyu(row, dataset, dictionary):
    '''
    Calculates the nyu value for a pair of species. 
    '''
    
    # Retrieve the species
    s1 = row["Species1"] # Species1
    s2 = row["Species2"]  # Species2

    # Retrieve the simulated values
    try:
        key = s1 + "@" + s2
        simulation_list = dictionary[key]
    except:
        key = s2 + "@" + s1
        simulation_list = dictionary[key]
    
    condition = (
        ((dataset["Species1"] == s1) & (dataset["Species2"] == s2))
    )
    # Retrieve the value from AbsDiff_z
    value = dataset.loc[condition, "AbsDiff_z"].values
    smaller = len(list(x for x in simulation_list if x < value))

    nyu = smaller / len(simulation_list)

    return nyu




#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@


#@ Step 0: the session.
#@#@#@#@#@#@#@#@#@#@#@#@#@

print("Starting session on", trait, "from", database_path, "database with", tree_path, "phylogeny")                      #@@ STDOUT PRINT




#@ Step 1: import primate traits database information (phylgeny and traits table)
#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@

#   [TRAITS TABLE] Import the traits table
db = pd.read_csv(database_path, sep='\t', header=0, index_col="SpeciesBROAD")

#   [PHYLOGENY] Import the tree
tree = dendropy.Tree.get(
    path = tree_path,
    schema='newick')

#   [PHYLOGENY] Extract the species
species_in_tree = []
species_objects = []

for i, t1 in enumerate(tree.taxon_namespace[:-1]):
    species_in_tree.append(str(t1).replace("'", "").replace(" ", "_"))
    species_objects.append((t1))

#   [PHYLOGENY] Load Distance Matrix

print("Loading distance matrix")                      #@@ STDOUT PRINT

distance_dictionary = {}

pdc = tree.phylogenetic_distance_matrix()

tree_species_combinations = list(cb(species_objects, 2))

distance_dictionary = {}

for x in tree_species_combinations:
    species_a = str(x[0]).replace(" ", "_").replace("'", "")
    species_b = str(x[1]).replace(" ", "_").replace("'", "")
    distance_dictionary[species_a + "," + species_b] = pdc(x[0], x[1])


#@ Step 2: Isolate the column realtive to the trait
#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#

print("Extracting trait values")                      #@@ STDOUT PRINT
#   [CHECKPOINT] Check if the trait exists in the table                                                         [!!!CHECKPOINT]
if trait not in db.columns:
    raise ValueError(f"The trait '{trait}' is not found in the dataset. Please verify the column names.")

trait_db = db[[trait]]
trait_db.dropna(subset=[trait], inplace=True)

#   Conversion to z-values of the observed values
obs_mean = trait_db[trait].mean()
obs_std = trait_db[trait].std()

print("Calculated observed z-values")                      #@@ STDOUT PRINT

trait_db["zvalues"] = (trait_db[trait] - obs_mean) / obs_std

#   Export the output

# Esporta il DataFrame come file CSV                                                                            [!!!OUTPUT] 
trait_obs_values = output_tag + ".single_values.csv"
trait_db.to_csv(trait_obs_values, sep='\t', index=True)                                                              

#@ Step 3: Calculate combined values and create the output database
#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@

#   Obtain species with trait values
species_with_values = trait_db.index.tolist()
species_intersection = set(species_in_tree).intersection(set(species_with_values))

#   Extract pairwise combinations
combinations = list(cb(species_intersection, 2))



#@ Step 4: Prune the tree with the species included in the species tree
#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@

# Taxa to we neex to keep
taxa_to_keep = [taxon for taxon in tree.taxon_namespace if str(taxon).replace(" ", "_").replace("'", "") in species_intersection]

# Copy of the tree
pruned_tree = tree.clone()

# Prune the tree
pruned_tree.prune_taxa([taxon for taxon in tree.taxon_namespace if taxon not in taxa_to_keep])

# Export the tree in the session folder                                                                                 [!!!OUTPUT]

trait_pruned_tree = output_tag + ".pruned_tree.nw"
pruned_tree.write(path= trait_pruned_tree, schema="newick")                                                          

#@ Step 5: Run the simulation (Rscript simulate.r ...)
#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@

print("Launching simulation with", model)                      #@@ STDOUT PRINT

command = f"Rscript simulate.r -d {trait_obs_values} -t {tree_path} -o {outfolder} -v {trait} -n {nsim} -m {model}"
os.system(command)

simulations_path = "/".join([output_tag + ".simulations.csv"])
simulations = pd.read_csv(simulations_path, sep=',', header=0, index_col="SpeciesBROAD")

#   Conversion to z-values of the simulated values

print("Calculating simulation z-values")                      #@@ STDOUT PRINT
simulations_z = simulations.apply(
    lambda col: (col - col.mean()) / col.std(), axis=0
)

simulations_z_path = output_tag + ".simulations.zscore.csv"
simulations_z.to_csv(simulations_z_path, sep='\t', index=True)                                                          

#@ Step 6: Calculate the modular difference between the simulated values
#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@

print("Calculating pairwise values for observed trait")                      #@@ STDOUT PRINT

sim_dict_d = simulations_z.to_dict(orient='index') 
sim_dict = {x : list(sim_dict_d[x].values()) for x in sim_dict_d.keys()}
sim_dict_key_combs = cb(sim_dict.keys(), 2)

pair2couple = {x[0] + "@" + x[1] : list(zip(sim_dict[x[0]], sim_dict[x[1]])) for x in sim_dict_key_combs}

pair2sim_dist = {
    x: [abs(y[0] - y[1]) for y in pair2couple[x]]  # Calcola la differenza assoluta per ogni coppia nella lista
    for x in pair2couple.keys()
}

#   The difference combination
results_table = pd.DataFrame(combinations, columns=["Species1", "Species2"])
results_table.insert(0, "Trait", trait)
results_table["PatristicDistance"] = results_table.apply(
    lambda row: set_distance(row, distance_dictionary=distance_dictionary), axis=1)

results_table["AbsDiff"] = results_table.apply(
    lambda row: calculate_abs_diff(row, trait_db=trait_db, trait=trait), axis=1)

results_table["AbsDiff_z"] = results_table.apply(
    lambda row: calculate_abs_diff(row, trait_db=trait_db, trait="zvalues"), axis=1)

#@ Step 7: Determine the pairwise nyu value
#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@

print("Cacluating species-to-species pairwise nyu")                      #@@ STDOUT PRINT

results_table["PairWise_NYU"] = results_table.apply(
    lambda row: calc_nyu(row, dataset=results_table, dictionary = pair2sim_dist), axis=1)

results_table = results_table.sort_values(by=["Species1", "Species2"], ascending=[True, True])
results_table = results_table.reset_index(drop=True)


nyu_table = output_tag + ".nyu.tab"
results_table.to_csv(nyu_table, sep='\t', index=True)                                                          

#@ Step 8: Calculate the trait nyu95
#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#

print("Calculating trait-wide nyu95")                      #@@ STDOUT PRINT

outliers = results_table[
    (results_table["PairWise_NYU"] < 0.025) | (results_table["PairWise_NYU"] > 0.975)
]

tops = results_table[
    (results_table["PairWise_NYU"] > 0.975)
]

bottoms = results_table[
    (results_table["PairWise_NYU"] < 0.025)
]

# Mostra le combinazioni filtrate
ol_table = output_tag + ".outliers.tab"
tops_table = output_tag + ".outliers.tops.tab"
bottoms_table = output_tag + ".outliers.bottoms.tab"

outliers.to_csv(ol_table, sep='\t', index=True)
tops.to_csv(tops_table, sep='\t', index=True)
bottoms.to_csv(bottoms_table, sep='\t', index=True)

nyu95 = float(outliers.shape[0])/results_table.shape[0]

nyu_95_txt = output_tag + ".nyu95.txt"

with open(nyu_95_txt, "w") as nyuout:
    print("\t".join([trait, args.model, str(nyu95), str(results_table.shape[0])]), file= nyuout)


print("\n\nDone, you can find the information in the following files:")                                                     #@@ STDOUT PRINT
print("#################################################################\n")                                                #@@ STDOUT PRINT
print("Single Observed values (with z-score):", "\t",trait_obs_values)                                                      #@@ STDOUT PRINT
print("Simulated values with", nsim, "simulations with", model, "\t", simulations_path)                                     #@@ STDOUT PRINT
print("Simulated values, z-scores", "\t", simulations_z_path)                                                               #@@ STDOUT PRINT
print("Total Pairwise nyu results", "\t", nyu_table)                                                                        #@@ STDOUT PRINT
print("Outliers (95%)", "\t", ol_table)                                                                                     #@@ STDOUT PRINT
print("Outliers Top Values (95%)", "\t", tops_table)                                                                        #@@ STDOUT PRINT
print("Outliers Bottoms Values (95%)", "\t", bottoms_table)                                                                 #@@ STDOUT PRINT