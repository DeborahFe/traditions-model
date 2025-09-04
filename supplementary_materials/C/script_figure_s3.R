# ==============================================================================
# SCRIPT TO STUDY INITIAL CONDITIONS EFFECTS AND GENE-CULTURE INTERACTIONS
# ==============================================================================
# This script generates Figure 7 (supplementary material section C) exploring 
# the effects of varying initial conditions on tradition formation. It examines
# how the initial frequencies of both male traits and female preferences interact
# to determine which trait becomes dominant in emergent traditions.
# ==============================================================================

#### opening dataframes from simulations ####

# ------------------------------------------------------------------------------
# DATA LOADING AND PREPROCESSING FOR INITIAL CONDITIONS ANALYSIS
# ------------------------------------------------------------------------------
# Load simulation results focused on the value of the first traditions for different
# initial conditions (different initial frequency for the male trait and for the
# female preference).

# Cultural Model datasets
trad_social_cultural <- read.csv("first_tradition_cultural_social.csv", sep = ";")
trad_social_cultural <- subset(trad_social_cultural, mutation_rate == 0.05) 
trad_null_cultural <- read.csv("first_tradition_cultural_null.csv", sep = ";")
trad_null_cultural <- subset(trad_null_cultural, mutation_rate == 0.05)

# Gene-Culture Model datasets
trad_social_geneculture <- read.csv("first_tradition_geneculture_social.csv", sep = ";")
trad_social_geneculture <- subset(trad_social_geneculture, mutation_rate == 0.05)
trad_null_geneculture <- read.csv("first_tradition_geneculture_null.csv", sep = ";")
trad_null_geneculture <- subset(trad_null_geneculture, mutation_rate == 0.05)


# Convert parameters to factors for proper categorical analysis

# Cultural Model factor conversions
trad_social_cultural$c <- as.factor(trad_social_cultural$c)                      
trad_social_cultural$mutation_rate <- as.factor(trad_social_cultural$mutation_rate)  

trad_null_cultural$c <- as.factor(trad_null_cultural$c)
trad_null_cultural$mutation_rate <- as.factor(trad_null_cultural$mutation_rate)

# Gene-Culture Model factor conversions
trad_social_geneculture$c <- as.factor(trad_social_geneculture$c)
trad_social_geneculture$mutation_rate <- as.factor(trad_social_geneculture$mutation_rate)

trad_null_geneculture$c <- as.factor(trad_null_geneculture$c)
trad_null_geneculture$mutation_rate <- as.factor(trad_null_geneculture$mutation_rate)

# ------------------------------------------------------------------------------
# SIMULATION PARAMETERS DOCUMENTATION
# ------------------------------------------------------------------------------
#### the parameters used in those simulations ####

## Global simulation parameters for initial conditions experiments
# These parameters differ from main analysis by including systematic variation
# of both initial trait and preference frequencies
p <- list(
  nb_class = c(1,1),              # Age structure: 1 juvenile + 1 adult class for each sex
  K = 1000,                       # Population carrying capacity
  female_strategy = "conformity",  # Female mate choice strategy (vs "koinophilia")
  survival = 0.9,                 # Per-timestep survival probability
  ageing = 0.1,                   # Probability to advance to next age class
  initial_trait_frequency = 0.5,  # Baseline trait frequency (varied in param_variation)
  initial_pref_frequency = 0.5,   # Baseline preference frequency (varied in param_variation)
  n_matings = 100,                # Maximum matings per male per timestep
  mutation_rate = 0.4,            # Baseline mutation rate (filtered to 0.05 in analysis)
  n_obs = 7,                      # Number of couples observed per female
  c = 0.7,                        # Baseline conformity strength (varied in param_variation)
  n_steps_removed = 1)            # Minimal exclusion period for tradition detection

## Parameter grid for social learning scenarios
# This design systematically explores the interaction between initial conditions
# and conformity strength to understand tradition emergence patterns
param_variation <- expand.grid(n_obs = 7,                                       
                               mutation_rate = c(1, 0.2, 0.05),                 
                               c = c(0.6, 0.7, 0.8),                            
                               initial_trait_frequency = c(0.1, 0.25, 0.5, 0.75, 0.9),      # Male trait frequencies
                               initial_pref_frequency = c(0.1, 0.25, 0.5, 0.75, 0.9))       # Female preference frequencies

## Parameter grid for null model (random choice baseline)
param_variation_null <- expand.grid(n_obs = 7,
                                    mutation_rate = c(1, 0.2, 0.05),
                                    c = 0.5,                                      # Fixed at random choice
                                    initial_trait_frequency = c(0.1, 0.25, 0.5, 0.75, 0.9),
                                    initial_pref_frequency = c(0.1, 0.25, 0.5, 0.75, 0.9))

## Other parameters
n_replicates <- 100                                        # Replications per parameter combination
n_steps <- 500                                            # Shorter simulation length for initial dynamics

# ------------------------------------------------------------------------------
# ANALYTICAL FUNCTIONS FOR TRADITION ANALYSIS
# ------------------------------------------------------------------------------
#### functions to compute the proportion of traditions ####

# Function to calculate average number of observable matings
average_n_matings <- function(n_obs, c_index, df_social) {
  df_subset <- df_social[which(df_social$n_obs == n_obs & df_social$c == c_index), ]
  average_matings <- round(mean(df_subset$mean_demos, na.rm = TRUE))
  return(average_matings)
}

# Function to calculate 95th percentile threshold from null model
quantile95 <- function(nobs, mutation, initial_trait_frequency, initial_pref_frequency, df_null) {
  # Filter null data for exact parameter match including initial conditions
  df_subset <- df_null[which(df_null$n_obs == nobs & 
                               df_null$mutation_rate == mutation & 
                               df_null$initial_trait_frequency == initial_trait_frequency & 
                               df_null$initial_pref_frequency == initial_pref_frequency), ]
  durations <- df_subset[, "tradition_duration"]
  quantile_95 <- quantile(durations, probs = 0.95)
  return(quantile_95)
}

# Function to calculate tradition proportion exceeding null expectations
prop_traditions <- function(n_obs, c_index, mutation_rate, initial_trait_frequency, initial_pref_frequency, df_null, df_social) {
  # Filter social learning data for exact parameter match
  df_subset <- df_social[which(df_social$n_obs == n_obs & 
                                 df_social$c == c_index & 
                                 df_social$mutation_rate == mutation_rate & 
                                 df_social$initial_trait_frequency == initial_trait_frequency & 
                                 df_social$initial_pref_frequency == initial_pref_frequency), ]
  # Get corresponding null threshold for these specific initial conditions
  quantile_95 <- quantile95(n_obs, mutation_rate, initial_trait_frequency, initial_pref_frequency, df_null)
  # Calculate proportion exceeding null expectation
  p_traditions <- nrow(subset(df_subset, tradition_duration > quantile_95)) / nrow(df_subset)
  return(p_traditions)
}

# ------------------------------------------------------------------------------
# TRADITION PROPORTION CALCULATIONS
# ------------------------------------------------------------------------------
#### computing the proportions of traditions ####

# Initialize result dataframes for both models
# Note: initial condition parameters are included in param_variation columns
col_names <- c("simu", "n_obs", "c", "mutation_rate", "n_demos", 
               "quantile_95", "prop_traditions")
prop_traditions_cultural_df <- as.data.frame(matrix(nrow = nrow(param_variation), ncol = length(col_names)))
colnames(prop_traditions_cultural_df) <- col_names
prop_traditions_geneculture_df <- as.data.frame(matrix(nrow = nrow(param_variation), ncol = length(col_names)))
colnames(prop_traditions_geneculture_df) <- col_names

# Calculate tradition proportions for Cultural Model
for (i in 1 : nrow(param_variation)) {
  prop_traditions_cultural_df[i, "simu"] <- i
  prop_traditions_cultural_df[i, colnames(param_variation)] <- param_variation[i,]  # Includes initial conditions
  prop_traditions_cultural_df[i, "n_demos"] <- average_n_matings(param_variation[i, "n_obs"], 
                                                                 param_variation[i, "c"], 
                                                                 trad_social_cultural)
  prop_traditions_cultural_df[i, "quantile_95"] <- quantile95(param_variation[i, "n_obs"], 
                                                              param_variation[i, "mutation_rate"], 
                                                              param_variation[i, "initial_trait_frequency"], 
                                                              param_variation[i, "initial_pref_frequency"], 
                                                              trad_null_cultural)
  prop_traditions_cultural_df[i, "prop_traditions"] <- prop_traditions(param_variation[i, "n_obs"], 
                                                                       param_variation[i, "c"], 
                                                                       param_variation[i, "mutation_rate"], 
                                                                       param_variation[i, "initial_trait_frequency"], 
                                                                       param_variation[i, "initial_pref_frequency"], 
                                                                       trad_null_cultural, 
                                                                       trad_social_cultural)
}

# Calculate tradition proportions for Gene-Culture Model
for (i in 1 : nrow(param_variation)) {
  prop_traditions_geneculture_df[i, "simu"] <- i
  prop_traditions_geneculture_df[i, colnames(param_variation)] <- param_variation[i,]
  prop_traditions_geneculture_df[i, "n_demos"] <- average_n_matings(param_variation[i, "n_obs"], 
                                                                    param_variation[i, "c"], 
                                                                    trad_social_geneculture)
  prop_traditions_geneculture_df[i, "quantile_95"] <- quantile95(param_variation[i, "n_obs"], 
                                                                 param_variation[i, "mutation_rate"], 
                                                                 param_variation[i, "initial_trait_frequency"], 
                                                                 param_variation[i, "initial_pref_frequency"], 
                                                                 trad_null_geneculture)
  prop_traditions_geneculture_df[i, "prop_traditions"] <- prop_traditions(param_variation[i, "n_obs"], 
                                                                          param_variation[i, "c"], 
                                                                          param_variation[i, "mutation_rate"], 
                                                                          param_variation[i, "initial_trait_frequency"], 
                                                                          param_variation[i, "initial_pref_frequency"], 
                                                                          trad_null_geneculture, 
                                                                          trad_social_geneculture)
}

# ------------------------------------------------------------------------------
# TRAIT-SPECIFIC TRADITION ANALYSIS
# ------------------------------------------------------------------------------
#### computing the proportion of first traditions being for trait 1 ####

# This section addresses the key question: when the first tradition emerges, which trait dominates?

# Initialize dataframes for trait-specific analysis
# These will contain the proportion of traditions that favor trait 1 vs trait 0
prop_trad_trait1_cultural <- param_variation
prop_trad_trait1_cultural$quantile_95 <- prop_traditions_cultural_df$quantile_95  # Use same thresholds
prop_trad_trait1_cultural$prop_trad1 <- NA  # Proportion of traditions favoring trait 1

# Calculate trait 1 tradition proportion for Cultural Model
for (i in 1 : nrow(prop_trad_trait1_cultural)) {
  # Filter data for exact parameter match including initial conditions
  sub_df <- trad_social_cultural[which(trad_social_cultural$c == prop_trad_trait1_cultural[i, "c"] & 
                                         trad_social_cultural$mutation_rate == prop_trad_trait1_cultural[i, "mutation_rate"] & 
                                         trad_social_cultural$initial_trait_frequency == prop_trad_trait1_cultural[i, "initial_trait_frequency"] & 
                                         trad_social_cultural$initial_pref_frequency == prop_trad_trait1_cultural[i, "initial_pref_frequency"]), ]
  
  # Use corresponding threshold from tradition analysis
  quantile_95 <- prop_traditions_cultural_df[i, "quantile_95"]
  
  # Calculate conditional proportion: P(tradition for trait 1 | tradition occurred)
  # Numerator: simulations with trait 1 traditions exceeding threshold
  # Denominator: all simulations with traditions (any trait) exceeding threshold
  prop_trad_trait1_cultural[i, "prop_trad1"] <- nrow(sub_df[which(sub_df$tradition_for_trait == 1 & 
                                                                    sub_df$tradition_duration > quantile_95), ]) / 
    nrow(sub_df[which(sub_df$tradition_duration > quantile_95), ])
}

# Calculate trait 1 tradition proportion for Gene-Culture Model
prop_trad_trait1_geneculture <- param_variation
prop_trad_trait1_geneculture$quantile_95 <- prop_traditions_geneculture_df$quantile_95
prop_trad_trait1_geneculture$prop_trad1 <- NA

for (i in 1 : nrow(prop_trad_trait1_geneculture)) {
  sub_df <- trad_social_geneculture[which(trad_social_geneculture$c == prop_trad_trait1_geneculture[i, "c"] & 
                                            trad_social_geneculture$mutation_rate == prop_trad_trait1_geneculture[i, "mutation_rate"] & 
                                            trad_social_geneculture$initial_trait_frequency == prop_trad_trait1_geneculture[i, "initial_trait_frequency"] & 
                                            trad_social_geneculture$initial_pref_frequency == prop_trad_trait1_geneculture[i, "initial_pref_frequency"]), ]
  
  quantile_95 <- prop_traditions_geneculture_df[i, "quantile_95"]
  
  prop_trad_trait1_geneculture[i, "prop_trad1"] <- nrow(sub_df[which(sub_df$tradition_for_trait == 1 & 
                                                                       sub_df$tradition_duration > quantile_95), ]) / 
    nrow(sub_df[which(sub_df$tradition_duration > quantile_95), ])
}

# ------------------------------------------------------------------------------
# HEATMAP VISUALIZATION SETUP
# ------------------------------------------------------------------------------
#### plot ####

# Load required visualization libraries
library(ggplot2)    # Core plotting functionality
library(patchwork)  # For combining multiple plots

# Create descriptive labels for conformity strength facets
# These will appear as panel headers showing different conformity levels
c_labs <- c()
for (i in 1 : length(unique(param_variation$c))){
  c_labs[i] <- paste("c =", unique(param_variation$c)[i])
}
names(c_labs) <- unique(param_variation$c)

# Convert variables to factors for proper discrete heatmap visualization

# Cultural Model data conversions
prop_trad_trait1_cultural$c <- as.factor(prop_trad_trait1_cultural$c)
prop_trad_trait1_cultural$mutation_rate <- as.factor(prop_trad_trait1_cultural$mutation_rate)
prop_trad_trait1_cultural$initial_trait_frequency <- as.factor(prop_trad_trait1_cultural$initial_trait_frequency)
prop_trad_trait1_cultural$initial_pref_frequency <- as.factor(prop_trad_trait1_cultural$initial_pref_frequency)

# Gene-Culture Model data conversions
prop_trad_trait1_geneculture$c <- as.factor(prop_trad_trait1_geneculture$c)
prop_trad_trait1_geneculture$mutation_rate <- as.factor(prop_trad_trait1_geneculture$mutation_rate)
prop_trad_trait1_geneculture$initial_trait_frequency <- as.factor(prop_trad_trait1_geneculture$initial_trait_frequency)
prop_trad_trait1_geneculture$initial_pref_frequency <- as.factor(prop_trad_trait1_geneculture$initial_pref_frequency)

# ------------------------------------------------------------------------------
# CULTURAL MODEL HEATMAP CONSTRUCTION
# ------------------------------------------------------------------------------
# Create heatmap showing interaction between initial trait and preference frequencies
# X-axis: initial trait frequency in males, Y-axis: initial preference frequency in females
# Color: proportion of emergent traditions that favor trait 1 (vs trait 0)
heatmap_traditions_cultural <- ggplot(prop_trad_trait1_cultural, aes(x = initial_trait_frequency, y = initial_pref_frequency, fill = prop_trad1)) +
  geom_tile() +
  facet_wrap(~ c, labeller = labeller(c = c_labs)) +
  labs(x = "Initial frequency of type 1 males",                           
       y = "Initial frequency of females preferring type 1 males",       
       fill = "Proportion of first traditions being for type 1 males",  
       title = "Cultural Model") +                                        
  theme_minimal() +
  theme(legend.position = "none",             
        strip.background = element_blank(),     
        title = element_text(size = 16),       
        axis.title.x = element_text(size = 14, margin = margin(t = 15)),  
        axis.title.y = element_text(size = 14, margin = margin(r = 15)), 
        axis.text.x = element_text(size = 12),  
        axis.text.y = element_text(size = 12),  
        legend.title = element_text(size = 12), 
        strip.text = element_text(size = 12)) 
heatmap_traditions_cultural

# ------------------------------------------------------------------------------
# GENE-CULTURE MODEL HEATMAP CONSTRUCTION
# ------------------------------------------------------------------------------
# Identical structure to Cultural Model heatmap for direct comparison
heatmap_traditions_geneculture <- ggplot(prop_trad_trait1_geneculture, aes(x = initial_trait_frequency, y = initial_pref_frequency, fill = prop_trad1)) +
  geom_tile() +
  facet_wrap(~ c, labeller = labeller(c = c_labs)) +
  labs(x = "Initial frequency of type 1 males", 
       y = "Initial frequency of females preferring type 1 males", 
       fill = "Proportion of first traditions being for type 1 males", 
       title = "Gene-Culture Model") +
  theme_minimal() +
  theme(legend.position = "bottom",     
        strip.background = element_blank(),
        title = element_text(size = 16),    
        axis.title.x = element_text(size = 14, margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, margin = margin(r = 15)), 
        axis.text.x = element_text(size = 12),  
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),  
        strip.text = element_text(size = 12)) 
heatmap_traditions_geneculture

# ------------------------------------------------------------------------------
# FINAL FIGURE ASSEMBLY AND EXPORT
# ------------------------------------------------------------------------------
# Combine both heatmaps vertically using patchwork grammar
# This allows direct comparison between Cultural and Gene-Culture model outcomes
heatmaps_cultural_culturalnd_geneculture <- heatmap_traditions_cultural / heatmap_traditions_geneculture + 
  plot_layout(axes = "collect") 
heatmaps_cultural_culturalnd_geneculture

# Export final combined figure as PDF for supplementary materials
pdf(file = "figure_s3.pdf",
    width = 9, height = 8.5)
heatmaps_cultural_culturalnd_geneculture
dev.off()

# Alternative SVG export option (commented out)
# svg(file = "gc_interaction.svg",
#     width = 9, height = 8.5)
# heatmaps_cultural_culturalnd_geneculture
# dev.off()