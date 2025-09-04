# ================================================================================
# FIGURE 4: PARAMETER SPACE EXPLORATION - TRADITION EMERGENCE HEATMAPS
# ================================================================================
#
# This script generates Figure 4 from the publication, which presents systematic
# exploration of the parameter space to identify conditions promoting tradition
# emergence. The figure consists of heatmaps showing the proportion of simulation
# replicates in which traditions were detected across different combinations of:
#
# - n_obs (observation capacity): Number of matings females can observe
# - c (copying fidelity): Probability of following learned preference  
# - mutation_rate (genetic transmission fidelity): Rate of trait mutation
#
# Key Features:
# - Cultural Model (top panel): Traditions require majority exaggeration
# - Gene-Culture Model (bottom panel): Traditions possible without majority exaggeration
# - Striped patterns indicate parameter combinations with majority exaggeration
# - Color intensity shows proportion of replicates with detected traditions
#
# ================================================================================

#### Import Simulation Results from Parameter Explorations ####

# These datasets contain tradition duration measurements from comprehensive
# parameter space explorations. Each row represents one simulation replicate
# with specific parameter combinations (100 replicates per combination).

# CULTURAL MODEL RESULTS
# Social model (c > 0.5)
trad_cultural_social <- read.csv("traditions_cultural_social.csv", sep = ";")
# Remove intermediate mutation rates and maximum n_obs for cleaner visualization
trad_cultural_social <- subset(trad_cultural_social, mutation_rate != 0.5 & mutation_rate != 0.4 & 
                          mutation_rate != 0.3 & n_obs != 224)

# Null model: Random mate choice baseline (c = 0.5) for tradition threshold calculation
trad_cultural_null <- read.csv("traditions_cultural_null.csv", sep = ";")
trad_cultural_null <- subset(trad_cultural_null, mutation_rate != 0.5 & mutation_rate != 0.4 & 
                        mutation_rate != 0.3 & n_obs != 224)

# GENE-CULTURE MODEL RESULTS  
# Social model (c > 0.5)
trad_geneculture_social <- read.csv("traditions_geneculture_social.csv", sep = ";")
trad_geneculture_social <- subset(trad_geneculture_social, mutation_rate != 0.5 & mutation_rate != 0.4 & 
                          mutation_rate != 0.3 & n_obs != 224)

# Null model: Random choice under encounter constraints
trad_geneculture_null <- read.csv("traditions_geneculture_null.csv", sep = ";")
trad_geneculture_null <- subset(trad_geneculture_null, mutation_rate != 0.5 & mutation_rate != 0.4 & 
                        mutation_rate != 0.3 & n_obs != 224)

# Convert parameters to factors for proper statistical analysis and visualization
# This ensures ggplot treats them as discrete categories rather than continuous variables

trad_cultural_social$n_obs <- as.factor(trad_cultural_social$n_obs)
trad_cultural_social$c <- as.factor(trad_cultural_social$c)
trad_cultural_social$mutation_rate <- as.factor(trad_cultural_social$mutation_rate)

trad_cultural_null$n_obs <- as.factor(trad_cultural_null$n_obs)
trad_cultural_null$c <- as.factor(trad_cultural_null$c)
trad_cultural_null$mutation_rate <- as.factor(trad_cultural_null$mutation_rate)

trad_geneculture_social$n_obs <- as.factor(trad_geneculture_social$n_obs)
trad_geneculture_social$c <- as.factor(trad_geneculture_social$c)
trad_geneculture_social$mutation_rate <- as.factor(trad_geneculture_social$mutation_rate)

trad_geneculture_null$n_obs <- as.factor(trad_geneculture_null$n_obs)
trad_geneculture_null$c <- as.factor(trad_geneculture_null$c)
trad_geneculture_null$mutation_rate <- as.factor(trad_geneculture_null$mutation_rate)

#### Parameter Space Definition ####

# Define the complete parameter space used in the simulations

# Base parameters (held constant across all simulations)
p <- list(
  nb_class = c(1,1),                    # Age structure: 1 juvenile + 1 adult class
  K = 1000,                             # Population carrying capacity (individuals)
  female_strategy = "conformity",       # Social learning strategy
  survival = 0.9,                       # Baseline survival probability per time step
  ageing = 0.1,                         # Probability of staying in same age class
  initial_trait_frequency = 0.5,        # Initial frequency of trait 1 (neutral start)
  n_matings = 100,                      # Maximum matings per male per generation
  mutation_rate = 0.1,                  # Default mutation rate (overridden by grid)
  n_obs = 6,                            # Default observation capacity (overridden by grid)
  c = 0.7,                              # Default copying fidelity (overridden by grid)
  n_steps_removed = 100                 # Burn-in period before tradition detection
)

# Parameter variation grids
# These define the systematic exploration of social learning and genetic parameters

# Complete parameter space (original full exploration)
param_variation_total <- expand.grid(
  n_obs = c(seq(1, 30, 2), 224),                       # Observation capacities
  mutation_rate = c(1, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05), # Genetic transmission fidelities
  c = seq(0.5, 1, 0.05)                                # Copying fidelities
)

# Reduced parameter space (used for visualization to avoid overcrowding)
param_variation <- expand.grid(
  n_obs = seq(1, 30, 2),                               # Regular observation capacity intervals
  mutation_rate = c(1, 0.2, 0.1, 0.05),               # Key mutation rate values
  c = seq(0.5, 1, 0.05)                               # Full copying fidelity range
)

# Null model parameter grids (for tradition threshold calculation)
param_variation_null_total <- expand.grid(
  n_obs = c(seq(1, 30, 2), 224),
  mutation_rate = c(1, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05),
  c = 0.5                                              # Fixed at random choice
)

param_variation_null <- expand.grid(
  n_obs = seq(1, 30, 2),
  mutation_rate = c(1, 0.2, 0.1, 0.05),
  c = 0.5
)

# Other parameters
n_replicates <- 100                     # Number of stochastic replicates per parameter combination
n_steps <- 1000                         # Maximum generations per simulation

#### Functions for Tradition Detection and Analysis ####

#' Calculate Average Number of Matings for Parameter Combination
#'
#' Computes the mean number of mating demonstrations available for observation
#' in simulations with given parameters. This reflects the actual demographic
#' conditions under which social learning operates.
#'
#' @param n_obs Observation capacity parameter
#' @param c_index Copying fidelity parameter
#' @param df_social Dataframe containing social model simulation results
#' @return Average number of matings across replicates
average_n_matings <- function(n_obs, c_index, df_social) {
  # Filter data to specific parameter combination
  df_subset <- df_social[which(df_social$n_obs == n_obs & df_social$c == c_index), ]
  
  # Calculate mean number of demonstrations across all replicates
  # This represents the actual observational context for social learning
  average_matings <- round(mean(df_subset$mean_demos, na.rm = TRUE))
  
  return(average_matings)
}

#' Calculate 95th Percentile Tradition Duration from Null Model
#'
#' Determines the tradition persistence threshold by analyzing the null model
#' (random mate choice) distribution. Traditions in social models must exceed
#' this threshold to be considered significant (not just random drift).
#'
#' @param nobs Observation capacity
#' @param mutation Mutation rate  
#' @param df_null Dataframe containing null model simulation results
#' @return 95th percentile of tradition durations under random choice
quantile95 <- function(nobs, mutation, df_null) {
  # Filter null model data to matching parameters
  df_subset <- df_null[which(df_null$n_obs == nobs & df_null$mutation_rate == mutation), ]
  
  # Extract tradition durations from random choice simulations
  durations <- df_subset[, "tradition_duration"]
  
  # Calculate 95th percentile as tradition detection threshold
  quantile_95 <- quantile(durations, probs = 0.95)
  
  return(quantile_95)
}

#' Calculate Proportion of Simulations with Detected Traditions
#'
#' Determines what fraction of simulation replicates showed tradition persistence
#' beyond the threshold established by the corresponding null model. This is the
#' primary outcome measure for tradition emergence analysis.
#'
#' @param n_obs Observation capacity
#' @param c_index Copying fidelity
#' @param mutation_rate Mutation rate
#' @param df_null Null model data for threshold calculation
#' @param df_social Social model data for tradition detection
#' @return Proportion of replicates with detected traditions (0 to 1)
prop_traditions <- function(n_obs, c_index, mutation_rate, df_null, df_social) {
  # Filter social model data to specific parameter combination
  df_subset <- df_social[which(df_social$n_obs == n_obs & df_social$c == c_index & 
                                 df_social$mutation_rate == mutation_rate), ]
  
  # Get tradition threshold from corresponding null model
  quantile_95 <- quantile95(n_obs, mutation_rate, df_null)
  
  # Count replicates exceeding the tradition threshold
  # This represents significant tradition emergence beyond random expectation
  p_traditions <- nrow(subset(df_subset, tradition_duration > quantile_95)) / nrow(df_subset)
  
  return(p_traditions)
}

#### Compute Tradition Proportions Across Parameter Space ####

# Initialize dataframes to store tradition analysis results
# These will contain proportions of traditions for each parameter combination

col_names <- c("simu", "n_obs", "c", "mutation_rate", "n_demos", 
               "quantile_95", "prop_traditions", "majority_exaggeration")

# CULTURAL MODEL ANALYSIS
prop_traditions_cultural <- as.data.frame(matrix(nrow = nrow(param_variation), ncol = length(col_names)))
colnames(prop_traditions_cultural) <- col_names

# GENE-CULTURE MODEL ANALYSIS  
prop_traditions_geneculture <- as.data.frame(matrix(nrow = nrow(param_variation), ncol = length(col_names)))
colnames(prop_traditions_geneculture) <- col_names

# Calculate tradition proportions for Cultural Model
for (i in 1:nrow(param_variation)) {
  prop_traditions_cultural[i, "simu"] <- i
  # Store parameter combination
  prop_traditions_cultural[i, colnames(param_variation)] <- param_variation[i,]
  
  # Calculate average number of demonstrations for this parameter set
  prop_traditions_cultural[i, "n_demos"] <- average_n_matings(
    param_variation[i, "n_obs"], param_variation[i, "c"], trad_cultural_social)
  
  # Determine tradition threshold from null model
  prop_traditions_cultural[i, "quantile_95"] <- quantile95(
    param_variation[i, "n_obs"], param_variation[i, "mutation_rate"], trad_cultural_null)
  
  # Calculate proportion of replicates with detected traditions
  prop_traditions_cultural[i, "prop_traditions"] <- prop_traditions(
    param_variation[i, "n_obs"], param_variation[i, "c"], 
    param_variation[i, "mutation_rate"], trad_cultural_null, trad_cultural_social)
}

# Calculate tradition proportions for Gene-Culture Model (identical structure)
for (i in 1:nrow(param_variation)) {
  prop_traditions_geneculture[i, "simu"] <- i
  prop_traditions_geneculture[i, colnames(param_variation)] <- param_variation[i,]
  
  prop_traditions_geneculture[i, "n_demos"] <- average_n_matings(
    param_variation[i, "n_obs"], param_variation[i, "c"], trad_geneculture_social)
  
  prop_traditions_geneculture[i, "quantile_95"] <- quantile95(
    param_variation[i, "n_obs"], param_variation[i, "mutation_rate"], trad_geneculture_null)
  
  prop_traditions_geneculture[i, "prop_traditions"] <- prop_traditions(
    param_variation[i, "n_obs"], param_variation[i, "c"], 
    param_variation[i, "mutation_rate"], trad_geneculture_null, trad_geneculture_social)
}

#### Functions for Majority Exaggeration Analysis ####

#' Calculate Conditional Choice Probability Function for Given Parameters
#'
#' Computes the probability that females choose type 1 males across all possible
#' population frequencies of choices for type 1 in the previous generation, 
#' implementing the hypergeometric sampling model.
#' This function generates the theoretical conformity curve for given parameters.
#'
#' @param n_obs Number of matings observed by each female
#' @param c Copying fidelity (probability of following learned preference)
#' @param n_matings Total number of matings available for observation
#' @return Vector of choice probabilities for each possible population frequency
conformity_function <- function(n_obs, c, n_matings) {
  
  # Create frequency spectrum from 0 to 1
  observable_ratios <- seq(0, 1, (1/n_matings))
  
  # Convert frequencies to counts of type 1 matings
  matings_type1 = round(observable_ratios*n_matings)
  
  # Adjust observation capacity if population has fewer demonstrations
  n_obs_temp <- min(n_obs, n_matings)
  
  # Majority detection threshold (hypergeometric distribution uses > operator)
  quorum_needed = floor(n_obs_temp/2)
  
  # Calculate probabilities of detecting majorities using hypergeometric distribution
  # P(observe majority for type 1)
  p_type1_quorum = phyper(quorum_needed, matings_type1, n_matings - matings_type1, 
                          n_obs_temp, lower.tail = FALSE)
  
  # P(observe majority for type 0)  
  p_type0_quorum = phyper(quorum_needed, n_matings - matings_type1, matings_type1, 
                          n_obs_temp, lower.tail = FALSE)
  
  # P(observe no clear majority - tie situation)
  null_pref = round(1 - p_type1_quorum - p_type0_quorum, 5)
  
  # Calculate final conditional choice probability incorporating copying fidelity and tie-breaking
  p_choice_type1 = c*p_type1_quorum +           # Choose type 1 when preferring type 1
    (1-c)*p_type0_quorum +       # Choose type 1 when preferring type 0 (copying error)
    (1-p_type1_quorum-p_type0_quorum)*0.5  # Random choice when tied observations
  
  return(p_choice_type1)
}

#' Detect Presence of Majority Exaggeration in Conformity Response
#'
#' Determines whether given parameters produce majority exaggeration, where
#' the probability of choosing the majority type exceeds its frequency in the
#' population. This is crucial for tradition emergence in the Cultural Model.
#'
#' @param n_obs Number of matings observed by each female
#' @param c Copying fidelity parameter
#' @param n_matings Total matings available for observation  
#' @return Character: "1" if majority exaggeration detected, "0" otherwise
exaggeration_existence <- function(n_obs, c, n_matings) {
  
  majority_exaggeration <- c()
  
  # Calculate theoretical choice probabilities across all possible frequencies
  observable_ratios <- seq(0, 1, (1/n_matings))
  p_choice_type1 <- conformity_function(n_obs, c, n_matings)
  
  # Focus analysis on majority frequencies (≥ 0.5) where exaggeration is relevant
  # Minority exaggeration is less important for tradition stability
  ratios <- which(observable_ratios >= 0.5)
  observable_ratios <- round(observable_ratios[ratios], 3)  # Avoid floating-point comparison errors
  p_choice_type1 <- round(p_choice_type1[ratios], 3)
  
  # Test for disproportionate copying: does choice probability exceed frequency?
  count <- 1
  while (count <= length(observable_ratios)) {
    if (p_choice_type1[count] <= observable_ratios[count]) {
      # No exaggeration at this frequency point
      majority_exaggeration <- "0"
      count <- count + 1
    } else {
      # Exaggeration detected - probability exceeds frequency
      # This indicates positive feedback that can sustain traditions
      majority_exaggeration <- "1" 
      count <- length(observable_ratios) + 1  # Exit loop early
    }
  }
  
  return(majority_exaggeration)
}

#### Add Majority Exaggeration Classification to Tradition Data ####

# CULTURAL MODEL: Classify each parameter combination for majority exaggeration
for (i in 1:nrow(prop_traditions_cultural)) {
  prop_traditions_cultural[i, "majority_exaggeration"] <- exaggeration_existence(
    prop_traditions_cultural[i, "n_obs"], 
    prop_traditions_cultural[i, "c"], 
    prop_traditions_cultural[i, "n_demos"]
  )
}

# Convert numeric codes to descriptive labels
prop_traditions_cultural$majority_exaggeration <- replace(prop_traditions_cultural$majority_exaggeration, 
                                                           prop_traditions_cultural$majority_exaggeration == "0", "no")
prop_traditions_cultural$majority_exaggeration <- replace(prop_traditions_cultural$majority_exaggeration, 
                                                           prop_traditions_cultural$majority_exaggeration == "1", "yes")

# GENE-CULTURE MODEL: Apply identical classification
for (i in 1:nrow(prop_traditions_geneculture)) {
  prop_traditions_geneculture[i, "majority_exaggeration"] <- exaggeration_existence(
    prop_traditions_geneculture[i, "n_obs"], 
    prop_traditions_geneculture[i, "c"], 
    prop_traditions_geneculture[i, "n_demos"]
  )
}

# Convert to descriptive labels
prop_traditions_geneculture$majority_exaggeration <- replace(prop_traditions_geneculture$majority_exaggeration, 
                                                           prop_traditions_geneculture$majority_exaggeration == "0", "no")
prop_traditions_geneculture$majority_exaggeration <- replace(prop_traditions_geneculture$majority_exaggeration, 
                                                           prop_traditions_geneculture$majority_exaggeration == "1", "yes")

#### Create Heatmap Visualizations ####

# Load libraries for advanced visualization
library(ggplot2)      # Grammar of graphics plotting system
library(ggpattern)    # Add patterns (stripes) to indicate majority exaggeration
library(viridis)      # Perceptually uniform color palettes
library(patchwork)    # Advanced plot composition and layout

# Create panel labels for mutation rate facets
# These appear as headers in the faceted heatmaps
mutation_rate_labs <- c()
for (i in 1:length(unique(param_variation$mutation_rate))){
  mutation_rate_labs[i] <- paste("μ =", unique(param_variation$mutation_rate)[i])
}
names(mutation_rate_labs) <- unique(param_variation$mutation_rate)

# Convert parameters to factors for proper ggplot handling
prop_traditions_cultural$n_obs <- as.factor(prop_traditions_cultural$n_obs)
prop_traditions_cultural$c <- as.factor(prop_traditions_cultural$c)
prop_traditions_cultural$mutation_rate <- as.factor(prop_traditions_cultural$mutation_rate)

prop_traditions_geneculture$n_obs <- as.factor(prop_traditions_geneculture$n_obs)
prop_traditions_geneculture$c <- as.factor(prop_traditions_geneculture$c)
prop_traditions_geneculture$mutation_rate <- as.factor(prop_traditions_geneculture$mutation_rate)

#### CULTURAL MODEL HEATMAP ####
# Shows tradition emergence patterns when mate choice is unconstrained
heatmap_cultural <- ggplot(prop_traditions_cultural, 
                                    aes(x = n_obs, y = c, pattern = majority_exaggeration, fill = prop_traditions)) +
  
  # Base heatmap tiles colored by proportion of traditions
  geom_tile() +
  
  # Apply perceptually uniform color scale optimized for scientific visualization
  # Higher intensity = higher proportion of simulations with detected traditions
  scale_fill_viridis(option = "F", direction = -1, begin = 0.2, end = 1) +
  
  # Configure x-axis to show observation capacity with mathematical notation
  scale_x_discrete(expression(n[obs])) +
  
  # Configure labels and title
  labs(fill = "Proportion of traditions    ",  # Extra spaces for legend alignment
       pattern = "Majority exaggeration", 
       title = "Cultural Model") +
  
  # Create separate panels for each mutation rate
  facet_wrap(vars(mutation_rate), labeller = labeller(mutation_rate = mutation_rate_labs), 
             ncol = length(mutation_rate_labs)) +
  
  # Add stripe patterns to indicate majority exaggeration
  # Striped tiles show parameter combinations that produce majority exaggeration
  # Solid tiles show combinations without majority exaggeration
  geom_tile_pattern(pattern_color = NA,           # No pattern border color
                    pattern_fill = "darkgrey",    # Gray stripes for visibility
                    pattern_angle = 45,           # Diagonal stripes
                    pattern_density = 0.05,       # Stripe spacing
                    pattern_spacing = 0.05,       # Pattern density
                    pattern_key_scale_factor = 1) +
  
  # Define pattern mapping: stripes only for majority exaggeration cases
  scale_pattern_manual(values = c("no" = "none", "yes" = "stripe")) +
  
  # Configure legend appearance  
  guides(pattern = guide_legend(override.aes = list(fill = "#F9DDC9FF"))) +
  
  # Apply publication-ready theme
  theme_minimal() +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 15)), 
        axis.title.y = element_text(size = 12, margin = margin(r = 15)),
        legend.title = element_text(size = 12),  
        strip.background = element_blank(),      # Clean facet appearance
        legend.position = "none",                # Hide legend (shown in Gene-Culture plot)
        plot.title = element_text(size = 17))

#### GENE-CULTURE MODEL HEATMAP ####  
# Shows tradition emergence under encounter-constrained mate choice
heatmap_geneculture <- ggplot(prop_traditions_geneculture, 
                                    aes(x = n_obs, y = c, pattern = majority_exaggeration, fill = prop_traditions)) +
  
  # Identical visual structure to Cultural Model for direct comparison
  geom_tile() +
  scale_fill_viridis(option = "F", direction = -1, begin = 0.2, end = 1) +
  scale_x_discrete(expression(n[obs])) +
  
  labs(fill = "Proportion of traditions    ", 
       pattern = "Majority exaggeration", 
       title = "Gene-Culture Model") +
  
  facet_wrap(vars(mutation_rate), labeller = labeller(mutation_rate = mutation_rate_labs), 
             ncol = length(mutation_rate_labs)) +
  
  # Apply identical pattern system for consistency
  geom_tile_pattern(pattern_color = NA, pattern_fill = "darkgrey", pattern_angle = 45,                 
                    pattern_density = 0.05, pattern_spacing = 0.05, pattern_key_scale_factor = 1) +
  scale_pattern_manual(values = c("no" = "none", "yes" = "stripe")) +
  
  # Show legend in bottom panel with sample fill color
  guides(pattern = guide_legend(override.aes = list(fill = "#751F58FF"))) +
  
  theme_minimal() +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 15)), 
        axis.title.y = element_text(size = 12, margin = margin(r = 15)),
        legend.title = element_text(size = 12),  
        strip.background = element_blank(),
        legend.position = "bottom",              # Show shared legend
        plot.title = element_text(size = 17))

#### Combine Heatmaps into Final Publication Figure ####

# Stack heatmaps vertically: Cultural Model above, Gene-Culture Model below
# This layout facilitates direct comparison of tradition emergence patterns
heatmaps <- heatmap_cultural / heatmap_geneculture

# Export to PDF format for publication
# Dimensions optimized for two-column academic journal layout
pdf(file = "figure4.pdf", width = 13, height = 8.5)
heatmaps
dev.off()

# Alternative SVG export (commented out)  
# Provides scalable vector graphics suitable for online publication
# svg(file = "figure4.svg", width = 13, height = 8.5)
# heatmaps
# dev.off()