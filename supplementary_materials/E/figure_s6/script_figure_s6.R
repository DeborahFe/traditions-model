# ==============================================================================
# SCRIPT FOR GENERATING TRADITION HEATMAPS (FIGURE 6)
# ==============================================================================
# This script generates Figure 10 (supplementary material) showing tradition 
# existence across Cultural and Gene-Culture models, for different values of 
# the survival rate s.
# ==============================================================================

#### opening dataframes from simulations ####

# ------------------------------------------------------------------------------
# DATA LOADING AND PREPROCESSING
# ------------------------------------------------------------------------------
# Load simulation results comparing social learning scenarios with null models
# Each dataset contains tradition duration data from multiple simulation runs

# Cultural Model datasets:
trad_social_cultural <- read.csv("traditions_cultural_social.csv", sep = ";")
trad_social_cultural <- subset(trad_social_cultural, mutation_rate != 0.5 & mutation_rate != 0.4 & mutation_rate != 0.3 & n_obs != 224)
trad_null_cultural <- read.csv("traditions_cultural_null.csv", sep = ";")
trad_null_cultural <- subset(trad_null_cultural, mutation_rate != 0.5 & mutation_rate != 0.4 & mutation_rate != 0.3 & n_obs != 224)

# Gene-Culture Model datasets:
trad_social_geneculture <- read.csv("traditions_geneculture_social.csv", sep = ";")
trad_social_geneculture <- subset(trad_social_geneculture, mutation_rate != 0.5 & mutation_rate != 0.4 & mutation_rate != 0.3 & n_obs != 224)
trad_null_geneculture <- read.csv("traditions_geneculture_null.csv", sep = ";")
trad_null_geneculture <- subset(trad_null_geneculture, mutation_rate != 0.5 & mutation_rate != 0.4 & mutation_rate != 0.3 & n_obs != 224)


# Convert key experimental parameters to factors for proper categorical analysis

# Cultural Model factor conversions
trad_social_cultural$n_obs <- as.factor(trad_social_cultural$n_obs)          
trad_social_cultural$c <- as.factor(trad_social_cultural$c)                   
trad_social_cultural$mutation_rate <- as.factor(trad_social_cultural$mutation_rate)  

trad_null_cultural$n_obs <- as.factor(trad_null_cultural$n_obs)
trad_null_cultural$c <- as.factor(trad_null_cultural$c)
trad_null_cultural$mutation_rate <- as.factor(trad_null_cultural$mutation_rate)

# Gene-Culture Model factor conversions
trad_social_geneculture$n_obs <- as.factor(trad_social_geneculture$n_obs)
trad_social_geneculture$c <- as.factor(trad_social_geneculture$c)
trad_social_geneculture$mutation_rate <- as.factor(trad_social_geneculture$mutation_rate)

trad_null_geneculture$n_obs <- as.factor(trad_null_geneculture$n_obs)
trad_null_geneculture$c <- as.factor(trad_null_geneculture$c)
trad_null_geneculture$mutation_rate <- as.factor(trad_null_geneculture$mutation_rate)

# ------------------------------------------------------------------------------
# SIMULATION PARAMETERS DOCUMENTATION
# ------------------------------------------------------------------------------
#### the parameters used in those simulations ####

## Global simulation parameters used across all runs
# These parameters define the basic population dynamics and biological constraints
p <- list(
  nb_class = c(1,1),              # Age structure: 1 juvenile + 1 adult class each for males/females
  K = 1000,                       # Population carrying capacity (total individuals)
  female_strategy = "conformity",  # Female mate choice strategy (vs "koinophilia")
  survival = 0.9,                 # Per-timestep survival probability (must be < 1)
  ageing = 0.1,                   # Probability to advance to next age class per timestep
  initial_trait_frequency = 0.5,  # Starting frequency of trait 1 in male population
  n_matings = 100,                # Maximum number of matings per male per timestep
  mutation_rate = 0.1,            # Rate of genetic mutation (heritability = 1 - mutation_rate)
  n_obs = 6,                      # Baseline number of couples observed per female
  c = 0.7,                        # Baseline conformity strength (0.5 = random, 1 = full conformity)
  n_steps_removed = 100)          # Number of initial timesteps excluded from tradition analysis

## Parameter grid for social learning simulations
# This grid defines all parameter combinations tested in the main experiment
param_variation <- expand.grid(n_obs = seq(1, 30, 2),                    
                               mutation_rate = 0.1,                       
                               c = seq(0.5, 1, 0.05),                    
                               survival = c(0.6, 0.7, 0.8, 0.9)) # variation of the survival rate

## Parameter grid for null model simulations  
# Null model uses random choice (c = 0.5) to establish baseline tradition durations
param_variation_null <- expand.grid(n_obs = seq(1, 30, 2),               
                                    mutation_rate = 0.1,                 
                                    c = 0.5,                             
                                    survival = c(0.6, 0.7, 0.8, 0.9))

## Other parameters
n_replicates <- 100                                        # Number of independent simulation runs per parameter set
n_steps <- 1000                                           # Total timesteps per simulation run

# ------------------------------------------------------------------------------
# ANALYTICAL FUNCTIONS FOR TRADITION ANALYSIS
# ------------------------------------------------------------------------------
#### functions to compute the proportion of traditions ####

# Function to calculate average number of observable matings for given parameters
average_n_matings <- function(n_obs, c_index, df_social) {
  # Filter data for specific parameter combination
  df_subset <- df_social[which(df_social$n_obs == n_obs & df_social$c == c_index), ]
  # Calculate mean number of demographic events (matings) across replicates
  average_matings <- round(mean(df_subset$mean_demos, na.rm = TRUE))
  return(average_matings)
}

# Function to calculate 95th percentile of tradition durations from null model
quantile95 <- function(nobs, mutation, survival, df_null) {
  # Filter null model data for matching parameter combination
  df_subset <- df_null[which(df_null$n_obs == nobs & df_null$mutation_rate == mutation & df_null$survival == survival), ]
  # Extract tradition duration values from all replicates
  durations <- df_subset[, "tradition_duration"]
  # Calculate 95th percentile: only 5% of null runs should exceed this duration
  quantile_95 <- quantile(durations, probs = 0.95)
  return(quantile_95)
}

# Function to calculate proportion of social learning runs that form "true traditions"
prop_traditions <- function(n_obs, c_index, mutation_rate, survival, df_null, df_social) {
  # Filter social learning data for specific parameter combination
  df_subset <- df_social[which(df_social$n_obs == n_obs & df_social$c == c_index & df_social$mutation_rate == mutation_rate), ]
  # Get corresponding null model threshold
  quantile_95 <- quantile95(n_obs, mutation_rate, survival, df_null)
  # Calculate proportion of social runs exceeding null threshold
  p_traditions <- nrow(subset(df_subset, tradition_duration > quantile_95)) / nrow(df_subset)
  return(p_traditions)
}

# ------------------------------------------------------------------------------
# TRADITION PROPORTION CALCULATIONS FOR BOTH MODELS
# ------------------------------------------------------------------------------
#### computing the proportions of traditions ####

# Initialize result dataframes to store tradition analysis results
# Each row will represent one parameter combination with calculated tradition proportion
col_names <- c("simu", "n_obs", "c", "mutation_rate", "survival", "n_demos", 
               "quantile_95", "prop_traditions")
prop_traditions_cultural_df <- as.data.frame(matrix(nrow = nrow(param_variation), ncol = length(col_names)))
colnames(prop_traditions_cultural_df) <- col_names
prop_traditions_geneculture_df <- as.data.frame(matrix(nrow = nrow(param_variation), ncol = length(col_names)))
colnames(prop_traditions_geneculture_df) <- col_names

# Calculate tradition proportions for Cultural Model
# Loop through each parameter combination in the experimental grid
for (i in 1 : nrow(param_variation)) {
  prop_traditions_cultural_df[i, "simu"] <- i                           # Simulation index
  prop_traditions_cultural_df[i, colnames(param_variation)] <- param_variation[i,]  # Copy parameter values
  # Calculate average number of matings for demographic context
  prop_traditions_cultural_df[i, "n_demos"] <- average_n_matings(param_variation[i, "n_obs"], 
                                                                 param_variation[i, "c"], 
                                                                 trad_social_cultural)
  # Get null model threshold for tradition definition
  prop_traditions_cultural_df[i, "quantile_95"] <- quantile95(param_variation[i, "n_obs"], 
                                                              param_variation[i, "mutation_rate"], 
                                                              param_variation[i, "survival"], 
                                                              trad_null_cultural)
  # Calculate main result: proportion of runs forming persistent traditions
  prop_traditions_cultural_df[i, "prop_traditions"] <- prop_traditions(param_variation[i, "n_obs"], 
                                                                       param_variation[i, "c"], 
                                                                       param_variation[i, "mutation_rate"], 
                                                                       param_variation[i, "survival"], 
                                                                       trad_null_cultural, 
                                                                       trad_social_cultural)
}

# Calculate tradition proportions for Gene-Culture Model
# Identical analysis pipeline but using Gene-Culture model data
for (i in 1 : nrow(param_variation)) {
  prop_traditions_geneculture_df[i, "simu"] <- i
  prop_traditions_geneculture_df[i, colnames(param_variation)] <- param_variation[i,]
  prop_traditions_geneculture_df[i, "n_demos"] <- average_n_matings(param_variation[i, "n_obs"], 
                                                                    param_variation[i, "c"], 
                                                                    trad_social_geneculture)
  prop_traditions_geneculture_df[i, "quantile_95"] <- quantile95(param_variation[i, "n_obs"], 
                                                                 param_variation[i, "mutation_rate"], 
                                                                 param_variation[i, "survival"], 
                                                                 trad_null_geneculture)
  prop_traditions_geneculture_df[i, "prop_traditions"] <- prop_traditions(param_variation[i, "n_obs"], 
                                                                          param_variation[i, "c"], 
                                                                          param_variation[i, "mutation_rate"], 
                                                                          param_variation[i, "survival"], 
                                                                          trad_null_geneculture, 
                                                                          trad_social_geneculture)
}

# ------------------------------------------------------------------------------
# HEATMAP VISUALIZATION SETUP
# ------------------------------------------------------------------------------
#### plots heatmaps traditions ####

# Load required visualization libraries
library(ggplot2)    # Core plotting functionality
library(ggpattern)  # Extended pattern capabilities for geom_tile (currently unused)
library(viridis)    # Perceptually uniform, colorblind-friendly color palettes  
library(patchwork)  # Grammar for combining multiple ggplot objects

# Create descriptive labels for faceting by exclusion period
# These labels will appear as panel headers in the final figure
survival_labs <- c()
for (i in 1 : length(unique(param_variation$survival))){
  survival_labs[i] <- paste("s =", unique(param_variation$survival)[i])
}
names(survival_labs) <- unique(param_variation$survival)

# Convert numerical variables to factors for proper discrete visualization
# This ensures correct color mapping and axis labeling in heatmaps

# Cultural Model data conversions
prop_traditions_cultural_df$n_obs <- as.factor(prop_traditions_cultural_df$n_obs)
prop_traditions_cultural_df$c <- as.factor(prop_traditions_cultural_df$c)
prop_traditions_cultural_df$mutation_rate <- as.factor(prop_traditions_cultural_df$mutation_rate)
prop_traditions_cultural_df$survival <- as.factor(prop_traditions_cultural_df$survival)

# Gene-Culture Model data conversions
prop_traditions_geneculture_df$n_obs <- as.factor(prop_traditions_geneculture_df$n_obs)
prop_traditions_geneculture_df$c <- as.factor(prop_traditions_geneculture_df$c)
prop_traditions_geneculture_df$mutation_rate <- as.factor(prop_traditions_geneculture_df$mutation_rate)
prop_traditions_geneculture_df$survival <- as.factor(prop_traditions_geneculture_df$survival)

# ------------------------------------------------------------------------------
# CULTURAL MODEL HEATMAP CONSTRUCTION
# ------------------------------------------------------------------------------
## CULTURAL MODEL HEATMAP
# Create heatmap showing tradition formation probability across parameter space
# X-axis: sample size (n_obs), Y-axis: conformity strength (c), Color: tradition proportion
heatmap_traditions_cultural <- ggplot(prop_traditions_cultural_df, aes(x = n_obs, y = c, fill = prop_traditions)) +
  geom_tile() +
  scale_fill_viridis(option = "F", direction = -1, begin = 0.2, end = 1) +
  scale_x_discrete(expression(n[obs])) +
  labs(fill = "Proportion of traditions    ", title = "Cultural Model") +
  facet_wrap(vars(survival), labeller = labeller(survival = survival_labs), ncol = length(survival_labs)) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 12,  margin = margin(t = 15)),  
        axis.title.y = element_text(size = 12, margin = margin(r = 15)),    
        legend.title = element_text(size = 12),                           
        strip.background = element_blank(),                               
        legend.position = "none",                                        
        plot.title = element_text(size = 17))                           

# ------------------------------------------------------------------------------
# GENE-CULTURE MODEL HEATMAP CONSTRUCTION  
# ------------------------------------------------------------------------------
## GENE-CULTURE MODEL HEATMAP
# Identical structure to Cultural Model heatmap for direct comparison
# Shows how genetic preferences interact with social learning to affect tradition formation
heatmap_traditions_geneculture <- ggplot(prop_traditions_geneculture_df, aes(x = n_obs, y = c, fill = prop_traditions)) +
  geom_tile() +
  scale_fill_viridis(option = "F", direction = -1, begin = 0.2, end = 1) +
  scale_x_discrete(expression(n[obs])) +
  labs(fill = "Proportion of traditions    ", title = "Gene-Culture Model") +
  facet_wrap(vars(survival), labeller = labeller(survival = survival_labs), ncol = length(survival_labs)) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 12,  margin = margin(t = 15)), 
        axis.title.y = element_text(size = 12, margin = margin(r = 15)),
        legend.title = element_text(size = 12),  
        strip.background = element_blank(),
        legend.position = "bottom",                                         
        plot.title = element_text(size = 17))

# ------------------------------------------------------------------------------
# FINAL FIGURE ASSEMBLY AND EXPORT
# ------------------------------------------------------------------------------

# Combine both heatmaps vertically using patchwork grammar
# "/" operator creates vertical stacking layout
heatmaps_cultural_and_geneculture <- heatmap_traditions_cultural / heatmap_traditions_geneculture

# Export final combined figure as PDF for supplementary materials
pdf(file = "figure_s6.pdf",
    width = 13, height = 8.5)
heatmaps_cultural_and_geneculture
dev.off()

# Alternative SVG export option (commented out)
# svg(file = "heatmaps_cultural_and_geneculture.svg",
#     width = 13, height = 8.5)
# heatmaps_cultural_and_geneculture
# dev.off()