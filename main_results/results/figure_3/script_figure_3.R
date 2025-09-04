# ================================================================================
# FIGURE 3: TEMPORAL DYNAMICS AND TRADITION DURATION ANALYSIS
# ================================================================================
#
# This script generates Figure 3 from the publication, which shows:
# 1. Examples of temporal dynamics of female mate choice and male traits over time
# 2. Duration of majority behaviors (traditions) under different social learning parameters
# 3. Comparison between Cultural Model and Gene-Culture Model dynamics
# 4. Effect of majority exaggeration on tradition persistence
#
# The figure consists of four panels arranged in a 2x2 layout:
# - Left panels: Temporal dynamics (null models vs social models)  
# - Right panels: Boxplots of tradition durations across parameter values
# - Top row: Cultural Model
# - Bottom row: Gene-Culture Model
#
# ================================================================================

#### SECTION 1: EXAMPLES OF TEMPORAL DYNAMICS ####

#### Define Model Parameters for Temporal Examples ####

# Base parameter set for Cultural Model temporal dynamics
# These parameters are chosen to demonstrate clear tradition emergence
p <- list(
  nb_class = c(1,1),                    # Age structure: 1 juvenile + 1 adult class  
  K = 1000,                             # Population carrying capacity (individuals)
  female_strategy = "conformity",       # Social learning strategy
  survival = 0.9,                       # Baseline survival probability per time step
  ageing = 0.1,                         # Probability of staying in same age class (low = rapid aging)
  initial_trait_frequency = 0.5,        # Initial frequency of trait 1 (neutral start)
  n_matings = 100,                      # Maximum matings per male per generation
  mutation_rate = 0.4,                  # Moderate genetic mutation rate
  n_obs = 12,                           # Number of matings observed by learning females
  c = 0.7,                              # Copying fidelity (moderate social learning strength)
  n_steps_removed = 50                  # Burn-in period before tradition detection
)

# Create parameter sets for different model scenarios
# Each scenario demonstrates specific aspects of the theoretical framework

# SOCIAL MODEL (Cultural Model)
p_C <- p

# NULL MODEL (Cultural Model)
p_C_null <- p
p_C_null$c = 0.5  # c = 0.5 means random choice (no social learning effect)

# SOCIAL MODEL (Gene-Culture):
p_GC <- p
p_GC$n_obs = 5  # Note: Different n_obs chosen to show interesting dynamics specific to each model 
                # (n_obs = 12 for the Cultural Model)

# NULL MODEL (Gene-Culture):
p_GC_null <- p_GC
p_GC_null$c <- 0.5 # Random choice

# Simulation duration for temporal analysis
n_time_steps <- 800                    # Long enough to observe tradition dynamics

#### Import Pre-computed Temporal Dynamics Data ####

# Import temporal evolution data from previous simulation runs
# These datasets contain time series of choice frequencies and trait frequencies
# Each file represents one simulation run with the parameters defined above

# Cultural Model temporal data
dynamics_cultural <- read.csv("dynamics_cultural.csv")           # Social learning active (c = 0.7)
dynamics_cultural_null <- read.csv("dynamics_cultural_null.csv") # Null model (c = 0.5)

# Gene-Culture Model temporal data  
dynamics_geneculture <- read.csv("dynamics_geneculture.csv")           # Social learning active (c = 0.7)
dynamics_geneculture_null <- read.csv("dynamics_geneculture_null.csv") # Null model (c = 0.5)

#### Create Temporal Dynamics Plots ####

# Load required libraries for advanced plotting
library(ggplot2)
library(patchwork) 

# Define consistent color scheme for biological interpretation
# Red: Female choice (culturally transmitted behavior)
# Green: Male type frequency (genetically inherited trait)
colors_legend <- c("Female choice" = "#ff2e58", "Male type" = "#058d74")

# CULTURAL MODEL - SOCIAL LEARNING ACTIVE
# This plot shows how traditions emerge when social learning operates
plot_dynamics_cultural <- ggplot(data = dynamics_cultural, aes(x = iteration)) + 
  
  # Add horizontal reference line at 50% frequency
  # This line helps identify when one trait/choice becomes the majority
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  
  # Plot female choice frequency over time (red line)
  # This represents the culturally transmitted behavior (mate preferences)
  geom_line(aes(y = choice, color = "Female choice"), linewidth = 0.4) +
  
  # Plot male trait frequency over time (green line)  
  # This represents the genetically inherited trait under selection
  geom_line(aes(y = trait, color = "Male type"), linewidth = 0.4) +
  
  # Configure plot aesthetics
  labs(color = "") +                    # Remove color legend title
  scale_color_manual(values = colors_legend) +
  scale_y_continuous("Frequency", breaks = seq(0, 1, 0.25), limits = c(0,1)) +
  scale_x_continuous("Time steps") +
  
  # Add parameter annotation for reproducibility
  annotate("text", x = 140, y = 0.9, 
           label = expression(paste(c, " = 0.7, ", n[obs], " = 12, ", μ, " = 0.4")), 
           colour = "grey40") +
  
  # Apply clean theme suitable for publication
  theme_minimal() +
  theme(legend.position = "none",      # Hide legend (will be shown in bottom plot)
        axis.title.x = element_text(size = 15, margin = margin(t = 15)),
        axis.title.y = element_text(size = 15, margin = margin(r = 15)), 
        legend.title = element_text(size = 13),  
        legend.text = element_text(size = 13),
        strip.text = element_text(size = 12))

# CULTURAL MODEL - NULL MODEL (NO SOCIAL LEARNING)
# This plot shows random drift when social learning is absent
plot_dynamics_cultural_null <- ggplot(data = dynamics_cultural_null, aes(x = iteration)) + 
  
  # Identical structure to social model for direct comparison
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_line(aes(y = choice, color = "Female choice"), linewidth = 0.4) +
  geom_line(aes(y = trait, color = "Male type"), linewidth = 0.4) +
  
  # Add panel title to distinguish models
  labs(color = "", title = "Cultural Model") +
  scale_color_manual(values = colors_legend) +
  scale_y_continuous("Frequency", breaks = seq(0, 1, 0.25), limits = c(0,1)) +
  scale_x_continuous("Time steps") +
  
  # Parameter annotation shows c = 0.5 (random choice)
  annotate("text", x = 140, y = 0.9, 
           label = expression(paste(c, " = 0.5, ", n[obs], " = 12, ", μ, " = 0.4")), 
           colour = "grey40") +
  
  theme_minimal() +
  theme(legend.position = "none",
        title = element_text(size = 15, margin = margin(t = 20)),    
        axis.title.x = element_text(size = 15, margin = margin(t = 15)),
        axis.title.y = element_text(size = 15, margin = margin(r = 15)), 
        legend.title = element_text(size = 13),  
        legend.text = element_text(size = 13),
        strip.text = element_text(size = 12))

# GENE-CULTURE MODEL - SOCIAL LEARNING ACTIVE  
# This plot demonstrates gene-culture coevolution dynamics
plot_dynamics_geneculture <- ggplot(data = dynamics_geneculture, aes(x = iteration)) + 
  
  # Identical visual structure for direct comparison with Cultural Model
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_line(aes(y = choice, color = "Female choice"), linewidth = 0.4) +
  geom_line(aes(y = trait, color = "Male type"), linewidth = 0.4) +
  
  labs(color = "") +
  scale_color_manual(values = colors_legend) +
  scale_y_continuous("Frequency", breaks = seq(0, 1, 0.25), limits = c(0,1)) +
  scale_x_continuous("Time steps") +
  
  # Note different n_obs parameter (5 vs 12) chosen to demonstrate model-specific dynamics
  annotate("text", x = 140, y = 0.9, 
           label = expression(paste(c, " = 0.7, ", n[obs], " = 5, ", μ, " = 0.4")), 
           colour = "grey40") +
  
  theme_minimal() +
  theme(legend.position = "bottom",    # Show legend in final plot
        axis.title.x = element_text(size = 15, margin = margin(t = 15)),
        axis.title.y = element_text(size = 15, margin = margin(r = 15)), 
        legend.title = element_text(size = 13),  
        legend.text = element_text(size = 13),
        strip.text = element_text(size = 12))

# GENE-CULTURE MODEL - NULL MODEL
# Shows random dynamics under encounter constraints but no social learning bias  
plot_dynamics_geneculture_null <- ggplot(data = dynamics_geneculture_null, aes(x = iteration)) + 
  
  # Identical structure to other plots for comparison
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_line(aes(y = choice, color = "Female choice"), linewidth = 0.4) +
  geom_line(aes(y = trait, color = "Male type"), linewidth = 0.4) +
  
  labs(color = "", title = "Gene-Culture Model") +
  scale_color_manual(values = colors_legend) +
  scale_y_continuous("Frequency", breaks = seq(0, 1, 0.25), limits = c(0,1)) +
  scale_x_continuous("Time steps") +
  
  annotate("text", x = 140, y = 0.9, 
           label = expression(paste(c, " = 0.5, ", n[obs], " = 5, ", μ, " = 0.4")), 
           colour = "grey40") +
  
  theme_minimal() +
  theme(legend.position = "none",
        title = element_text(size = 15, margin = margin(t = 20)),    
        axis.title.x = element_text(size = 15, margin = margin(t = 15)),
        axis.title.y = element_text(size = 15, margin = margin(r = 15)), 
        legend.title = element_text(size = 13),  
        legend.text = element_text(size = 13),
        strip.text = element_text(size = 12))

# Combine Cultural Model plots (null above, social below)
plots_dynamics_cultural <- plot_dynamics_cultural_null / plot_dynamics_cultural + 
  plot_layout(axes = "collect")

# Combine Gene-Culture Model plots (null above, social below)  
plots_dynamics_geneculture <- plot_dynamics_geneculture_null / plot_dynamics_geneculture + 
  plot_layout(axes = "collect")

#### SECTION 2: TRADITION DURATION ANALYSIS (BOXPLOTS) ####

#### Import Tradition Duration Data ####

# These datasets contain tradition duration measurements from systematic parameter explorations
# Each row represents one simulation replicate with specific parameter combinations
# tradition_duration = number of time steps a majority choice persisted

# Gene-Culture Model tradition data (constrained mate choice)
traditions_geneculture <- read.csv("traditions_geneculture.csv", sep = ";")
# Subset to specific parameter combination for focused analysis
traditions_geneculture_sub <- subset(traditions_geneculture, mutation_rate == 0.4 & n_obs == 5)

# Cultural Model tradition data (unconstrained mate choice)
traditions_cultural <- read.csv("traditions_cultural.csv", sep = ";")
# Use identical parameter subset for direct comparison
traditions_cultural_sub <- subset(traditions_cultural, mutation_rate == 0.4 & n_obs == 5)

#### Functions to Detect Majority Exaggeration ####

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

#' Detect Presence of Majority Exaggeration in Conformity Curve
#'
#' Determines whether a given parameter combination produces majority exaggeration. 
#'
#' @param n_obs Number of matings observed by each female  
#' @param c Copying fidelity parameter
#' @param n_matings Total matings available for observation
#' @return Character: "1" if majority exaggeration exists, "0" if not
exaggeration_existence <- function(n_obs, c, n_matings) {
  
  majority_exaggeration <- c()
  
  # Calculate choice probabilities across frequency spectrum
  observable_ratios <- seq(0, 1, (1/n_matings))
  p_choice_type1 <- conformity_function(n_obs, c, n_matings)
  
  # Focus analysis on majority frequencies (≥ 0.5) where exaggeration matters most
  ratios <- which(observable_ratios >= 0.5)
  observable_ratios <- round(observable_ratios[ratios], 3)  # Round to avoid floating-point errors
  p_choice_type1 <- round(p_choice_type1[ratios], 3)
  
  # Test for exaggeration: does choice probability exceed population frequency?
  count <- 1
  while (count <= length(observable_ratios)) {
    if (p_choice_type1[count] <= observable_ratios[count]) {
      # No exaggeration detected at this frequency
      majority_exaggeration <- "0" 
      count <- count + 1
    } else {
      # Exaggeration detected - choice probability exceeds frequency
      majority_exaggeration <- "1"
      count <- length(observable_ratios) + 1  # Exit loop early
    }
  }
  
  return(majority_exaggeration)
}

#### Add Majority Exaggeration Classification to Datasets ####

# CULTURAL MODEL: Classify each parameter combination for majority exaggeration
for (i in 1:nrow(traditions_cultural_sub)) {
  traditions_cultural_sub[i, "exaggeration"] <- exaggeration_existence(
    traditions_cultural_sub[i, "n_obs"], 
    traditions_cultural_sub[i, "c"], 
    round(traditions_cultural_sub[i, "mean_demos"])  # Use observed number of demonstrations
  )
}

# Convert numeric codes to descriptive labels
traditions_cultural_sub$exaggeration <- replace(traditions_cultural_sub$exaggeration, 
                                         traditions_cultural_sub$exaggeration == "0", "no")
traditions_cultural_sub$exaggeration <- replace(traditions_cultural_sub$exaggeration, 
                                         traditions_cultural_sub$exaggeration == "1", "yes")

# GENE-CULTURE MODEL: Apply identical classification
for (i in 1:nrow(traditions_geneculture_sub)) {
  traditions_geneculture_sub[i, "exaggeration"] <- exaggeration_existence(
    traditions_geneculture_sub[i, "n_obs"], 
    traditions_geneculture_sub[i, "c"], 
    round(traditions_geneculture_sub[i, "mean_demos"])
  )
}

# Convert to descriptive labels
traditions_geneculture_sub$exaggeration <- replace(traditions_geneculture_sub$exaggeration, 
                                         traditions_geneculture_sub$exaggeration == "0", "no")
traditions_geneculture_sub$exaggeration <- replace(traditions_geneculture_sub$exaggeration, 
                                         traditions_geneculture_sub$exaggeration == "1", "yes")

#### Create Tradition Duration Boxplots ####

# Load additional libraries for boxplot creation
library(ggplot2)
library(cowplot)      # Additional plot arrangement functions
library(RColorBrewer) # Color palettes for scientific visualization
library(patchwork)    # Advanced plot composition

# Define color scheme for majority exaggeration classification
# Blue: Majority exaggeration present (theoretically favorable for traditions)
# Gray: No majority exaggeration (traditions less likely in Cultural Model)
colors_legend <- c("yes" = "#48bbe8", "no" = "grey60")

# Convert copying fidelity to factor for proper ggplot handling
traditions_cultural_sub$c <- as.factor(traditions_cultural_sub$c)
traditions_geneculture_sub$c <- as.factor(traditions_geneculture_sub$c)

# CULTURAL MODEL BOXPLOT
# Shows tradition duration distribution across copying fidelity values
boxplot_cultural <- ggplot(traditions_cultural_sub, aes(x = c, y = tradition_duration, color = exaggeration)) + 
  
  # Create boxplots with gray fill for readability
  geom_boxplot(fill = "grey90", lwd = 0.65) + 
  
  # Apply color scheme based on majority exaggeration classification
  scale_color_manual(values = colors_legend, labels = c("No", "Yes")) +
  
  # Configure axes and legend
  labs(y = "Duration of majority choices",
       color = "Majority exaggeration") +
  
  # Add parameter annotation for context
  annotate("text", x = 9, y = 50, 
           label = expression(paste(n[obs], " = 5, ", μ, " = 0.4")), 
           colour = "grey40") +
  
  theme_minimal() +
  theme(legend.position = "none",      # Hide legend (shown in Gene-Culture plot)
        axis.title.x = element_text(size = 15, margin = margin(t = 15)),
        axis.title.y = element_text(size = 15, margin = margin(r = 15)), 
        legend.title = element_text(size = 13),  
        legend.text = element_text(size = 13),
        strip.text = element_text(size = 12))

# GENE-CULTURE MODEL BOXPLOT  
# Identical structure for direct comparison with Cultural Model
boxplot_geneculture <- ggplot(traditions_geneculture_sub, aes(x = c, y = tradition_duration, color = exaggeration)) + 
  
  geom_boxplot(fill = "grey90", lwd = 0.65) + 
  scale_color_manual(values = colors_legend, labels = c("No", "Yes")) +
  
  labs(y = "Duration of majority choices",
       color = "Majority exaggeration") +
  
  annotate("text", x = 9, y = 50, 
           label = expression(paste(n[obs], " = 5, ", μ, " = 0.4")), 
           colour = "grey40") +
  
  theme_minimal() +
  theme(legend.position = "bottom",    # Show legend in final plot
        axis.title.x = element_text(size = 15, margin = margin(t = 15)),
        axis.title.y = element_text(size = 15, margin = margin(r = 15)), 
        legend.title = element_text(size = 13),  
        legend.text = element_text(size = 13),
        strip.text = element_text(size = 12))

#### SECTION 3: COMBINE ALL PLOTS INTO PUBLICATION FIGURE ####

# Create final 2x2 layout combining temporal dynamics and tradition analysis
# Left column: Temporal dynamics (Cultural above, Gene-Culture below)
# Right column: Boxplots (Cultural above, Gene-Culture below)
temporal_dynamics_and_boxplots <- (plots_dynamics_cultural | boxplot_cultural) / 
  (plots_dynamics_geneculture | boxplot_geneculture)

# Display combined figure
temporal_dynamics_and_boxplots

# Export to PDF for publication
# Dimensions optimized for two-column academic journal layout  
pdf(file="figure3.pdf", width = 11, height = 9)
temporal_dynamics_and_boxplots
dev.off()

# Alternative SVG export (commented out)
# Provides scalable vector graphics for online publication
# svg(file="figure3.svg", width = 11, height = 9)
# temporal_dynamics_and_boxplots  
# dev.off()