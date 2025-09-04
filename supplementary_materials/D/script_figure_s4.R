# ==============================================================================
# SCRIPT FOR GENERATING MATING SYSTEM VARIATION ANALYSIS (FIGURE 8)
# ==============================================================================
# This script generates Figure 8 (supplementary material section D) examining
# the effect of restricting the number of matings per male on tradition emergence.
# It explores how different mating systems - from strict monogamy (m=1) to 
# polygyny (m>1) - influence the ability of social learning to 
# generate cultural traditions. This addresses the robustness of 
# results across diverse biological mating systems.
# ==============================================================================

#### open dataframes ####

# ------------------------------------------------------------------------------
# DATA LOADING FOR MATING SYSTEM ANALYSIS
# ------------------------------------------------------------------------------
# Load simulation results focused on Gene-Culture Model only
# This analysis specifically examines how genetic-cultural interactions
# are affected by constraints on male reproductive success

# Gene-Culture Model datasets
trad_null_geneculture <- read.csv("traditions_geneculture_null.csv", sep = ";")
trad_social_geneculture <- read.csv("traditions_geneculture_social.csv", sep = ";")

# Convert experimental parameters to factors for proper categorical analysis

# Social learning scenario factor conversions
trad_social_geneculture$c <- as.factor(trad_social_geneculture$c)                      
trad_social_geneculture$mutation_rate <- as.factor(trad_social_geneculture$mutation_rate) 

# Null model factor conversions
trad_null_geneculture$c <- as.factor(trad_null_geneculture$c)                       
trad_null_geneculture$mutation_rate <- as.factor(trad_null_geneculture$mutation_rate)

# ------------------------------------------------------------------------------
# SIMULATION PARAMETERS DOCUMENTATION
# ------------------------------------------------------------------------------
#### the parameters used in those simulations ####

## Global simulation parameters for mating system analysis
# These parameters are held constant while systematically varying the maximum
# number of matings per male, and the mutation rate, and c 
p <- list(
  nb_class = c(1,1),              # Age structure: 1 juvenile + 1 adult class per sex
  K = 1000,                       # Population carrying capacity
  female_strategy = "conformity",  # Female mate choice strategy (fixed)
  survival = 0.9,                 # Per-timestep survival probability
  ageing = 0.1,                   # Age class transition probability
  initial_trait_frequency = 0.5,  # Balanced initial male trait distribution
  # initial_pref_frequency = 0.5, # Initial female preference (commented - not used in Gene-Culture Model)
  n_matings = 100,                # Baseline maximum matings (varied in param_variation)
  mutation_rate = 0.4,            # Baseline genetic mutation rate (varied in analysis)
  n_obs = 7,                      # Number of couples observed per female (fixed)
  c = 0.7,                        # Baseline conformity strength (varied in analysis)
  n_steps_removed = 100)          # Time steps excluded from tradition analysis

## Parameter grid for mating system analysis
# This design systematically varies mating constraints alongside key transmission parameters
# to understand their interactive effects on tradition formation
param_variation <- expand.grid(n_obs = 7,                             
                               mutation_rate = c(1, 0.2, 0.05),       
                               c = c(0.6, 0.7, 0.8),                 
                               n_matings = seq(1, 20, 1))              # Mating system variation: monogamy to polygyny

## Parameter grid for null model (random choice baseline)
param_variation_null <- expand.grid(n_obs = 7,
                                    mutation_rate = c(1, 0.2, 0.05),
                                    c = 0.5,                            # Fixed at random choice
                                    n_matings = seq(1, 20, 1))

## Other parameters
n_replicates <- 100                                        # Replications per parameter combination
n_steps <- 1000                                           # Full simulation length for tradition analysis

# ------------------------------------------------------------------------------
# ANALYTICAL FUNCTIONS FOR MATING-CONSTRAINED TRADITION ANALYSIS
# ------------------------------------------------------------------------------
#### functions to compute the proportion of traditions ####

# Function to calculate average number of observable matings
average_n_matings <- function(n_obs, c_index, nb_matings_males, df_social) {
  # Filter data for exact parameter match including mating constraints
  df_subset <- df_social[which(df_social$n_obs == n_obs & 
                                 df_social$c == c_index & 
                                 df_social$n_matings == nb_matings_males), ]
  # Calculate mean demographic events (matings) across replicates
  average_matings <- round(mean(df_subset$mean_demos, na.rm = TRUE))
  return(average_matings)
}

# Function to calculate 95th percentile threshold from null model
quantile95 <- function(nobs, mutation, nb_matings_males, df_null) {
  # Filter null data for exact parameter match including mating system
  df_subset <- df_null[which(df_null$n_obs == nobs & 
                               df_null$mutation_rate == mutation & 
                               df_null$n_matings == nb_matings_males), ]
  durations <- df_subset[, "tradition_duration"]
  quantile_95 <- quantile(durations, probs = 0.95)
  return(quantile_95)
}

# Function to calculate tradition proportion exceeding null expectations
prop_traditions <- function(n_obs, c_index, mutation_rate, nb_matings_males, df_null, df_social) {
  # Filter social learning data for exact parameter match including mating constraints
  df_subset <- df_social[which(df_social$n_obs == n_obs & 
                                 df_social$c == c_index & 
                                 df_social$mutation_rate == mutation_rate & 
                                 df_social$n_matings == nb_matings_males), ]
  # Get corresponding null threshold for same mating system
  quantile_95 <- quantile95(n_obs, mutation_rate, nb_matings_males, df_null)
  # Calculate proportion of social runs exceeding mating-system-specific null threshold
  p_traditions <- nrow(subset(df_subset, tradition_duration > quantile_95)) / nrow(df_subset)
  return(p_traditions)
}

# ------------------------------------------------------------------------------
# TRADITION PROPORTION CALCULATIONS ACROSS MATING SYSTEMS
# ------------------------------------------------------------------------------
#### computing the proportions of traditions ####

# Initialize result dataframe for mating system analysis
# Each row represents one combination of transmission and mating system parameters
col_names <- c("simu", "n_obs", "c", "mutation_rate", "n_matings", "n_demos", 
               "quantile_95", "prop_traditions")
prop_traditions_geneculture_df <- as.data.frame(matrix(nrow = nrow(param_variation), ncol = length(col_names)))
colnames(prop_traditions_geneculture_df) <- col_names

# Calculate tradition proportions across all parameter combinations
for (i in 1 : nrow(param_variation)) {
  prop_traditions_geneculture_df[i, "simu"] <- i                              # Simulation index
  prop_traditions_geneculture_df[i, colnames(param_variation)] <- param_variation[i,]  # Copy all parameter values
  
  # Calculate average number of matings for demographic context
  # This shows how mating constraints affect actual reproductive events
  prop_traditions_geneculture_df[i, "n_demos"] <- average_n_matings(param_variation[i, "n_obs"], 
                                                                    param_variation[i, "c"], 
                                                                    param_variation[i, "n_matings"], 
                                                                    trad_social_geneculture)
  
  # Get null model threshold specific to this mating system
  # Different mating systems may have different baseline tradition durations
  prop_traditions_geneculture_df[i, "quantile_95"] <- quantile95(param_variation[i, "n_obs"], 
                                                                 param_variation[i, "mutation_rate"], 
                                                                 param_variation[i, "n_matings"], 
                                                                 trad_null_geneculture)
  
  # Calculate main result: tradition proportion for this specific mating system
  prop_traditions_geneculture_df[i, "prop_traditions"] <- prop_traditions(param_variation[i, "n_obs"], 
                                                                          param_variation[i, "c"], 
                                                                          param_variation[i, "mutation_rate"], 
                                                                          param_variation[i, "n_matings"], 
                                                                          trad_null_geneculture, 
                                                                          trad_social_geneculture)
}

# ------------------------------------------------------------------------------
# VISUALIZATION SETUP AND PLOT CONSTRUCTION
# ------------------------------------------------------------------------------
#### plot ####

# Load required visualization libraries
library(ggplot2)  # Core plotting functionality  
library(viridis)  # Perceptually uniform, colorblind-friendly color palettes

# Create descriptive labels for mutation rate facets
# These will appear as panel headers showing different genetic transmission fidelities
mutation_rate_labs <- c()
for (i in 1 : length(unique(param_variation$mutation_rate))){
  mutation_rate_labs[i] <- paste("µ =", unique(param_variation$mutation_rate)[i])
}
names(mutation_rate_labs) <- unique(param_variation$mutation_rate)

# Convert variables to factors for proper discrete color mapping
prop_traditions_geneculture_df$c <- as.factor(prop_traditions_geneculture_df$c)                      # Conformity strength (discrete colors)
prop_traditions_geneculture_df$mutation_rate <- as.factor(prop_traditions_geneculture_df$mutation_rate)  # Mutation rate (faceting)

# ------------------------------------------------------------------------------
# MAIN PLOT: MATING SYSTEM EFFECTS ON TRADITION FORMATION
# ------------------------------------------------------------------------------
# Create line plot showing how tradition probability changes with mating system constraints
# X-axis: maximum matings per male (mating system), Y-axis: tradition proportion
# Color: conformity strength, Facets: mutation rate (genetic transmission fidelity)
plot_nb_matings <- ggplot(data = prop_traditions_geneculture_df, aes(x = n_matings, y = prop_traditions, color = c)) +
  geom_point() +                                                 
  geom_line(lwd = 0.5) +                                        
  facet_grid(mutation_rate ~., labeller = labeller(mutation_rate = mutation_rate_labs)) +
  scale_y_continuous("Proportion of traditions in simulations", 
                     breaks = seq(0, 1, 0.2), limits = c(0,1)) +
  scale_x_continuous("Maximal number of matings per male", 
                     breaks = seq(1, 20, 2), limits = c(1,20)) +
  scale_color_viridis(option = "D", discrete = TRUE, direction = -1, begin = 0.4, end = 0.8) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 15)),  
        axis.title.y = element_text(size = 12, margin = margin(r = 15)),  
        axis.text.x = element_text(size = 12),                           
        axis.text.y = element_text(size = 12),                           
        legend.title = element_text(size = 12),                          
        legend.text = element_text(size = 12),                         
        strip.text = element_text(size = 12))                           


# ------------------------------------------------------------------------------
# FIGURE EXPORT FOR SUPPLEMENTARY MATERIALS
# ------------------------------------------------------------------------------
# Export final plot as PDF for supplementary material publication
pdf(file = "figure_s4.pdf",
    width = 7, height = 7)
plot_nb_matings
dev.off()

# Alternative SVG export option (commented out)
# svg(file = "nb_matings_males.svg",
#     width = 7, height = 7)
# plot_nb_matings
# dev.off()