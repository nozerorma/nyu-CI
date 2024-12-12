#!/usr/bin/env Rscript

# Load required libraries
library(optparse)
library(ape)
library(geiger)
library(phytools)
library(phylolm)

# Define command-line options
option_list <- list(
  make_option(c("-d", "--datafile"), type = "character", default = "test/Body_mass.OU.single_values.csv",
              help = "Path to the trait data file [default= %default]", metavar = "character"),
  make_option(c("-t", "--treefile"), type = "character", default = "test/Body_mass.OU.pruned_tree.nw",
              help = "Path to the pruned tree file [default= %default]", metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = "test",
              help = "Directory to save outputs [default= %default]", metavar = "character"),
  make_option(c("-v", "--trait"), type = "character", default = "Body_mass",
              help = "Name of the trait column [default= %default]", metavar = "character"),
  make_option(c("-n", "--nsim"), type = "integer", default = 100,
              help = "Number of simulation cycles [default= %default]", metavar = "integer"),
  make_option(c("-m", "--model"), type = "character", default = "OUrandomRoot",
              help = "Evolutionary model to use (BM or OUrandomRoot) [default= %default]", metavar = "character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Assign variables from parsed arguments
trait_obs_values <- opt$datafile
trait_pruned_tree <- opt$treefile
output_dir <- opt$outdir
trait_value <- opt$trait
nsim <- opt$nsim  # Number of simulations
model_type <- opt$model  # Model type (BM or OUrandomRoot)

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the log file with trait name and model type included
log_file <- file.path(output_dir, paste0(trait_value, ".", model_type, ".simulation.log"))

# Start sinking output to the log file
sink(log_file, split = TRUE)  # Redirect output and print to console as well

# Load the traits data
cat("Loading trait data from file:", trait_obs_values, "\n")
traits_df <- read.table(
  trait_obs_values,
  sep = "\t", header = TRUE,
  stringsAsFactors = FALSE
)

# Check if the required columns exist
if (!all(c(trait_value, "SpeciesBROAD") %in% colnames(traits_df))) {
  stop("Error: Missing required columns in the trait data file.")
}

cat("Successfully loaded trait data with", nrow(traits_df), "rows.\n")

# Associate trait values with species
trait_data <- setNames(traits_df[[trait_value]], traits_df$SpeciesBROAD)

# Load the pruned tree
cat("Loading phylogenetic tree from file:", trait_pruned_tree, "\n")
pruned_tree <- read.tree(trait_pruned_tree)

# Adjust the tree to include only species present in the trait data
cat("Pruning tree to match species in trait data.\n")
pruned_tree <- drop.tip(pruned_tree, setdiff(pruned_tree$tip.label, names(trait_data)))
trait_data <- trait_data[names(trait_data) %in% pruned_tree$tip.label]

cat("Tree pruned. Remaining tips:", length(pruned_tree$tip.label), "\n")

# Fit the selected evolutionary model
if (model_type == "BM") {
  cat("Fitting BM model using geiger.\n")
  model <- fitContinuous(pruned_tree, trait_data, model = "BM")
  params <- model$opt
  cat("Fitted BM model parameters:\n")
  print(params)
  
  # Simulate traits under BM model
  simulations <- fastBM(pruned_tree, sig2 = params$sigsq, nsim = nsim)
  
} else if (model_type == "OUrandomRoot") {
  cat("Fitting OU model using phylolm (OUrandomRoot).\n")
  
  # Prepare data for phylolm
  phylolm_data <- data.frame(trait = trait_data, species = names(trait_data))
  
  # Fit OUrandomRoot model
  ou_model <- phylolm(
    trait ~ 1,            # Model formula
    data = phylolm_data,
    phy = pruned_tree,
    model = "OUrandomRoot"
  )
  
  # Log parameters
  alpha <- ou_model$optpar  # Alpha parameter for the OU process
  sigma <- sqrt(ou_model$sigma2)  # Sigma parameter
  theta <- ou_model$coefficients[1]  # Intercept (mean trait value)
  
  cat("Fitted OUrandomRoot model parameters:\n")
  cat("Alpha:", alpha, "\n")
  cat("Sigma:", sigma, "\n")
  cat("Theta:", theta, "\n")
  
  # Simulate traits under OU model
  simulations <- tryCatch({
    replicate(nsim, rTraitCont(pruned_tree, model = "OU", alpha = alpha, sigma = sigma, theta = theta))
  }, error = function(e) {
    cat("Error during OU simulation:", e$message, "\n")
    NULL
  })
  
  if (is.null(simulations)) {
    stop("Simulations failed. Check model parameters and input data.")
  }
} else {
  stop("Unsupported model type. Please choose 'BM' or 'OUrandomRoot'.")
}

# Verify that simulations are consistent with the tree's tip labels
if (nrow(simulations) != length(pruned_tree$tip.label)) {
  stop("Error: Number of simulated values does not match the number of tree tips.")
}

# Convert to a data frame with species names as rows and simulations as columns
simulations_df <- data.frame(SpeciesBROAD = pruned_tree$tip.label, simulations)

# Save the simulations to a CSV file named after the trait
simulations_file <- file.path(output_dir, paste0(trait_value, ".", model_type, ".simulations.csv"))
write.table(
  simulations_df, file = simulations_file,
  sep = ",",
  row.names = FALSE,  # Do not save row names
  col.names = TRUE    # Include column names
)

cat("Simulations saved to:", simulations_file, "\n")

# Stop sinking and restore default output
sink()
cat("Logging completed. Log saved to:", log_file, "\n")