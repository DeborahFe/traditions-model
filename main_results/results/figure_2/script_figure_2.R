# ================================================================================
# FIGURE 2: CONFORMITY AND CHOICE CURVES ANALYSIS
# ================================================================================
#
# This script generates Figure 2 from the publication, which shows the
# relationship between population-level choice frequencies at t and choice
# probabilities at t+1 under different conformist social learning parameters.
# 
# The figure contains two panels:
# a) Conformity curves (identical for both Cultural and Gene-Culture Models)
# b) Choice curves (specific to Gene-Culture Model with encounter constraints)
#
# ================================================================================

#### Load Required Libraries ####

library(ggplot2)
library(viridis)

#### Data Import and Model Parameters ####

# Import pre-computed simulation results for Gene-Culture Model
# These data represent means from 1000 simulation replicates, each lasting one time step
# Each row represents one unique combination of tested parameters
# Don't forget to set the data folder as working directory
df_geneculture <- read.csv("GC_model_means.csv")

# Baseline parameters used across all simulations
parameters <- list(
  nb_class = c(1,1),                    # Age structure: 1 juvenile + 1 adult class
  K = 1000,                             # Population carrying capacity (individuals)
  female_strategy = "conformity",       # Social learning strategy
  survival = 0.9,                       # Baseline survival probability per time step
  ageing = 0.1,                         # Probability of remaining in same age class
  initial_trait_frequency = 0.5,        # Initial frequency of trait 1 in population
  n_matings = 100,                      # Maximum matings per male per generation
  mutation_rate = 0.1,                  # Genetic mutation rate per generation
  n_obs = 10,                           # Number of matings observed by learning females
  c = 0.5                               # Copying fidelity
)

# Theoretical maximum number of observable matings
# This represents the total number of mating demonstrations available when population
# is at carrying capacity with optimal demographic structure
n_observable_matings <- 224 # Pre-computed from demographic simulations with K=1000

# Define parameter space for systematic exploration
# This grid covers the full range of social learning intensities and initial conditions
param_variation <- expand.grid(
  n_obs = c(0, 1, seq(5, 30, 5), n_observable_matings), # Observation capacity range
  mutation_rate = 0,                                    # No mutation (pure cultural dynamics)
  c = seq(0.5, 1, 0.1),                                # Copying fidelity range
  initial_trait_frequency = seq(0, 1, 1/n_observable_matings) # Full frequency spectrum
)

n_rep <- 1000  # Number of simulation replicates per parameter combination

# Extract unique parameter values for systematic analysis
c_vector <- unique(param_variation$c)      # All copying fidelity values tested
n_obs_vector <- unique(param_variation$n_obs)  # All observation capacities tested
nb_of_plots <- length(c_vector)           # Number of panels needed in figure

# Create frequency ratio vector for theoretical calculations
# This represents all possible frequencies of trait 1 among observable matings
ratios <- seq(0, 1, 1/n_observable_matings)
possible_n_matings_1 <- seq(0, n_observable_matings, 1)  # Possible counts of trait 1 matings

#### Analytical Function for Conformity Curves ####

#' Calculate Theoretical Conformity Response
#' 
#' Computes the probability that a female will choose a type 1 male given her
#' learned preference, based on hypergeometric sampling of social demonstrations.
#' This function implements the mathematical model described in the Methods section.
#'
#' @param n_matings Total number of matings available for observation
#' @param n_matings_1 Number of matings involving type 1 males
#' @param n_obs Number of matings each female can observe (observation capacity)
#' @param c Copying fidelity (probability of following learned preference)
#' @return Probability of choosing type 1 male
#'
#' @details
#' The function implements the following biological process:
#' 1. Females sample n_obs matings from n_matings available demonstrations
#' 2. They detect majority choice among their sample (hypergeometric distribution)
#' 3. They develop preference for majority type (or random if tied)
#' 4. They apply learned preference with probability c when encountering males
theoretical_conformity_response <- function(n_matings, n_matings_1, n_obs, c) {
  
  # Ensure observation capacity doesn't exceed available demonstrations
  # This prevents mathematical errors when female density is low
  n_obs <- min(n_obs, n_matings)
  
  # Calculate threshold for majority detection
  # floor() ensures we use discrete counts for hypergeometric distribution
  quorum_needed = floor(n_obs/2)
  
  # Probability of observing majority for type 1 males
  # Uses hypergeometric distribution: sampling n_obs from population of n_matings
  # where n_matings_1 are "successes" (type 1 matings)
  p_maj1 = phyper(quorum_needed, n_matings_1, n_matings - n_matings_1, 
                  n_obs, lower.tail = FALSE)
  
  # Probability of observing majority for type 0 males (symmetric calculation)
  p_maj0 = phyper(quorum_needed, n_matings - n_matings_1, n_matings_1, 
                  n_obs, lower.tail = FALSE)
  
  # Probability of observing no clear majority (tie situation)
  # Rounded to avoid negative infinitesimal probabilities from floating-point arithmetic
  p_no_maj = round(1 - p_maj1 - p_maj0, 5)
  
  # Calculate preference probabilities based on majority detection
  p_pref1 <- p_maj1 + 0.5*p_no_maj  # Prefer type 1 if majority or half of ties
  p_pref0 <- p_maj0 + 0.5*p_no_maj  # Prefer type 0 if majority or half of ties
  
  # Calculate final choice probability incorporating copying fidelity
  # This represents the conditional probability P(C₁|E₁) from the paper
  p_choice1 <- c*p_pref1 + (1-c)*p_pref0
  
  return(p_choice1)
}

#### Generate Analytical Conformity Curves ####

# Initialize storage for analytical conformity curves data
# These curves show the relationship between population frequency and individual choice probability
columns_df <- c("p_choice_1","n_matings_1", "ratio_matings_1","n_obs", "c") 
analytical_curves <- data.frame(matrix(nrow = 0, ncol = length(columns_df))) 
colnames(analytical_curves) = columns_df
plots_list <- list()  # Store individual plot data for potential separate analysis

# Generate conformity curves for each copying fidelity level
count <- 1
while (count <= nb_of_plots) {
  tested_c <- c_vector[count]  # Current copying fidelity being analyzed
  
  # Initialize storage for current copying fidelity level
  one_plot_df <- data.frame(matrix(nrow = 0, ncol = length(columns_df))) 
  colnames(one_plot_df) = columns_df
  
  # Generate curves for each observation capacity (n_obs)
  prop <- 1
  for (nobs in n_obs_vector) {
    
    # Initialize storage for current observation capacity
    one_nobs_df <- data.frame(matrix(nrow = length(possible_n_matings_1), ncol = length(columns_df))) 
    colnames(one_nobs_df) = columns_df
    
    # Calculate choice probabilities across all possible trait frequencies
    for (i in 1:length(possible_n_matings_1)) { 
      one_nobs_df[i, "n_matings_1"] <- possible_n_matings_1[i]  # Count of type 1 matings
      one_nobs_df[i, "n_obs"] <- nobs                          # Observation capacity
      one_nobs_df[i, "ratio_matings_1"] <- ratios[i]           # Frequency of type 1 matings
      one_nobs_df[i, "c"] <- tested_c                          # Copying fidelity
      
      # Calculate theoretical choice probability using analytical model
      one_nobs_df[i, "p_choice_1"] <- theoretical_conformity_response(
        n_observable_matings, possible_n_matings_1[i], nobs, tested_c)
    } 
    
    # Compile data for current observation capacity
    one_plot_df <- rbind(one_plot_df, one_nobs_df)
    plots_list[[count]] <- one_plot_df
    prop <- prop + 1
  }
  
  # Compile data for current copying fidelity level
  analytical_curves <- rbind(analytical_curves, one_plot_df)
  count <- count + 1
}

#### Create Analytical Conformity Curves Plot ####

# Convert numeric variables to factors for proper ggplot handling
analytical_curves$n_obs <- as.factor(analytical_curves$n_obs)
analytical_curves$c <- as.factor(analytical_curves$c)

# Create panel labels for copying fidelity levels
# These will appear as facet headers in the final plot
c_labs <- c()
for (i in 1:length(unique(param_variation$c))){
  c_labs[i] <- paste("c =", unique(param_variation$c)[i])
}
names(c_labs) <- unique(param_variation$c)

# Subset data for cleaner visualization
# Remove n_obs = 0 and intermediate c values for clarity
analytical_curves_sub <- subset(analytical_curves, n_obs != 0 & c != 0.7 & c!= 0.9)

# Create conformity curves plot (Panel A)
plot_analytical <- ggplot(analytical_curves_sub, aes(x = ratio_matings_1, y = p_choice_1, color = n_obs)) +
  
  # Add background regions to highlight majority exaggeration
  # Blue regions: where P(C₁|E₁) > F(C₁) (majority exaggeration occurs)
  # Gray regions: where P(C₁|E₁) ≤ F(C₁) (no majority exaggeration)
  annotate(geom = "polygon", x = c(0, 0.5, 0.5), y = c(0, 0, 0.5), 
           fill = "#48bbe8", alpha = 0.3) +
  annotate(geom = "polygon", x = c(0.5, 0.5, 1), y = c(0.5, 1, 1), 
           fill = "#48bbe8", alpha = 0.3) +
  annotate(geom = "polygon", x = c(0.5, 1, 1), y = c(0.5, 0.5, 1), 
           fill = "gray80", alpha = 0.2) +
  annotate(geom = "polygon", x = c(0, 0, 0.5), y = c(0, 0.5, 0.5), 
           fill = "gray80", alpha = 0.2) +
  
  # Plot conformity curves as lines (different colors for different n_obs)
  geom_line(lwd = 0.75) +
  
  # Add reference lines for interpretation
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +  # Random choice baseline
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey") +  # Equal frequency line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +  # No bias line
  
  # Create separate panels for each copying fidelity level
  facet_grid(~ c, labeller = labeller(c = c_labs)) +
  
  # Configure axes with scientific notation
  scale_x_continuous(expression("Frequency in the population of choices for type 1 males" ~ F(C[1]) ~ "at" ~ t-1), 
                     breaks = seq(0, 1, 0.25), limits = c(0,1), labels = seq(0, 1, 0.25)) +
  scale_y_continuous(expression(P(C[1] ~ "|" ~ E[1]) ~ "at" ~ t), 
                     breaks = seq(0, 1, 0.25), limits = c(0,1), labels = seq(0, 1, 0.25)) +
  
  # Apply color scheme
  scale_color_viridis(option = "D", discrete = TRUE, direction = -1, begin = 0.4, end = 0.8, 
                      labels = c(1, seq(5, 30, 5), expression("">= n[mat]))) +
  
  # Configure legend
  guides(color = guide_legend(title = expression(n[obs]))) +
  
  # Add panel title
  ggtitle("a. Conformity curves (identical in the Cultural and Gene–Culture Models)") +
  
  # Apply theme
  theme_minimal() +
  theme(strip.background = element_blank(),
        title = element_text(size = 15, margin = margin(t = 20)),    
        axis.title.x = element_text(size = 15, margin = margin(t = 15)),
        axis.title.y = element_text(size = 15, margin = margin(r = 15)), 
        axis.text.x = element_text(size = 12),  
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),  
        legend.text = element_text(size = 11),
        strip.text = element_text(size = 12)) 

# Display the analytical plot
plot_analytical

#### Create Choice Curves Plot from Gene-Culture Model Simulations ####

# Prepare simulation data for plotting (same factor conversions)
df_geneculture$n_obs <- as.factor(df_geneculture$n_obs)
df_geneculture$c <- as.factor(df_geneculture$c)

# Use identical panel labels for consistency between plots
c_labs <- c()
for (i in 1:length(unique(param_variation$c))){
  c_labs[i] <- paste("c =", unique(param_variation$c)[i])
}
names(c_labs) <- unique(param_variation$c)

# Subset simulation data for visualization
# Remove boundary conditions and intermediate c values for clarity
df_geneculture_sub <- subset(df_geneculture, n_obs != 0 & 
                               initial_trait_frequency != 0 & 
                               initial_trait_frequency != 1 & 
                               c != 0.7 & c!= 0.9)

# Create choice curves plot (Panel B)
# This uses simulation results rather than analytical calculations
# because encounters in Gene-Culture Model depend on complex frequency dynamics
plot_simulations <- ggplot(df_geneculture_sub, aes(x = female_choice_first, y = female_choice_last, color = n_obs)) +
  
  # Add identical background regions for direct comparison with Panel A
  annotate(geom = "polygon", x = c(0, 0.5, 0.5), y = c(0, 0, 0.5), 
           fill = "#48bbe8", alpha = 0.3) +
  annotate(geom = "polygon", x = c(0.5, 0.5, 1), y = c(0.5, 1, 1), 
           fill = "#48bbe8", alpha = 0.3) +
  annotate(geom = "polygon", x = c(0.5, 1, 1), y = c(0.5, 0.5, 1), 
           fill = "gray80", alpha = 0.2) +
  annotate(geom = "polygon", x = c(0, 0, 0.5), y = c(0, 0.5, 0.5), 
           fill = "gray80", alpha = 0.2) +
  
  # Plot choice curves from simulation results
  geom_line(lwd = 0.75) +
  
  # Add identical reference lines for comparison
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  
  # Use identical faceting structure
  facet_grid(~ c, labeller = labeller(c = c_labs)) +
  
  # Configure axes
  scale_x_continuous(expression("Frequency in the population of choices for type 1 males" ~ F(C[1]) ~ "at" ~ t-1), 
                     breaks = seq(0, 1, 0.25), limits = c(0,1), labels = seq(0, 1, 0.25)) +
  scale_y_continuous(expression(P(C[1]) ~ "at" ~ t), 
                     breaks = seq(0, 1, 0.25), limits = c(0,1), labels = seq(0, 1, 0.25)) +
  
  # Apply color scheme
  scale_color_viridis(option = "D", discrete = TRUE, direction = -1, begin = 0.4, end = 0.8, 
                      labels = c(1, seq(5, 30, 5), expression("">= n[mat]))) +
  
  # Configure legend
  guides(color = guide_legend(title = expression(n[obs]))) +
  
  # Add panel title
  ggtitle("b. Choice curves (Gene–Culture Model)") +
  
  # Apply theme
  theme_minimal() +
  theme(strip.background = element_blank(),
        title = element_text(size = 15, margin = margin(t = 20)),    
        axis.title.x = element_text(size = 15, margin = margin(t = 15)),
        axis.title.y = element_text(size = 15, margin = margin(r = 15)), 
        axis.text.x = element_text(size = 12),  
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),  
        legend.text = element_text(size = 11),
        strip.text = element_text(size = 12)) 

# Display the simulation plot
plot_simulations

#### Combine Plots into Publication Figure ####

# Combine both panels into single publication-ready figure
# Uses patchwork syntax for professional multi-panel layout
figure_conf_curves <- plot_analytical / plot_simulations + 
  plot_layout(axes = "collect") +     # Align axes across panels
  plot_layout(guides = 'collect')     # Use shared legend

# Display combined figure
figure_conf_curves

# Export figure to PDF format for publication
# Dimensions optimized for two-column academic journal format
pdf(file="figure2.pdf", width = 13, height = 8.5)
figure_conf_curves
dev.off()

# Alternative SVG export (commented out)
# SVG format provides scalable vector graphics for online publication
# svg(file="figure2.svg", width = 13, height = 8.5)
# figure_conf_curves
# dev.off()