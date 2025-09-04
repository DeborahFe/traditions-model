# ================================================================================
# AUTHORS: [D. FEDERICO, J-B. FERDY, F-X. DECHAUME-MONCHARMONT, A. POCHEVILLE]
# DATE: [2025]
# ASSOCIATED PUBLICATION: [Gene-Culture Coevolution Favours the Emergence of 
# Traditions in Mate Choice through Conformist Social Learning]
# 
# ================================================================================
# STUDY'S GOALS
# ================================================================================
#
# This script implements individual-based models to study how conformist social 
# learning interacts with genetic inheritance to influence the emergence of 
# traditions in mate choice behaviour.
#
# THEORETICAL FRAMEWORK:
# - Cultural Model: Pure social learning without resource constraints
# - Gene-Culture Model: Social learning with encounter constraints
# - Both models explore how copying fidelity (c) and observation capacity (n_obs)
#   affect tradition formation
#
# MAIN RESEARCH QUESTION:
# Under what conditions do traditions emerge in mate choice, and how does
# gene-culture coevolution facilitate or constrain traditions?
#
# ================================================================================
# SCRIPT ORGANIZATION
# ================================================================================
#
# POPULATION MANAGEMENT:
#   - init_pop()                    : Initialize age-structured population
#
# DYNAMICS FOR THE TWO MODELS:
#   - cultural_model()              : Cultural Model (uncontrained mate choice)
#   - geneculture_model()           : Gene-Culture Model (constrained mate choice)
#
# RESULTS OPTIONS (in the two models functions below):
#   - return_tradition_only = FALSE : Return raw temporal results
#   - return_tradition_only = TRUE  : Return the duration of the first tradition
#
# SIMULATION CONTROL:
#   - parameters_exploration()      : Systematic parameter space exploration (only
#                                     for return_tradition_only = TRUE)
#
# ================================================================================
# KEY PARAMETERS
# ================================================================================
#
# SOCIAL LEARNING:
#   - n_obs          : Number of matings observed by learning females
#   - c              : Copying fidelity (probability of following learned preference)
#   - female_strategy: "conformity" or "koinophilia" learning strategy (only the 
#                       conformity strategy is studied in this article)
#
# GENETIC INHERITANCE:
#   - mutation_rate  : Rate of trait mutation per generation
#
# POPULATION DYNAMICS:
#   - K              : Population carrying capacity
#   - survival       : Baseline survival probability
#   - ageing         : Probability of advancing to next age class
#
# GENE-CULTURE MODEL SPECIFIC:
#   - n_matings      : Maximum matings per male per generation
#
# ================================================================================
# USAGE EXAMPLES
# ================================================================================
#
# # Cultural Model simulation:
# result_cultural <- cultural_model(tmax = 1000, parameters = p, 
#                                   return_tradition_only = TRUE)
#
# # Gene-Culture Model simulation:
# result_geneculture <- geneculture_model(tmax = 1000, parameters = p, 
#                                         return_tradition_only = TRUE)
#
# # Parameter space exploration:
# result_cultural <- parameters_exploration(cultural_model, param_grid, p, 
#                                           "result_cultural.csv", n_replicates = 100)
# result_geneculture <- parameters_exploration(geneculture_model, param_grid, p, 
#                                              "results_geneculture.csv", n_replicates = 100)
#
# ================================================================================
# FUNCTIONS
# ================================================================================

#### Population Initialization ####

#' Initialize Age-Structured Population
#' 
#' Creates a population matrix with realistic age structure based on survival
#' probabilities, initializes genetic traits, and establishes initial mating
#' pairs to provide a foundation for social learning in subsequent generations.
#' 
#' @param parameters List containing model parameters
#' @return List containing males matrix, females matrix, and age classes definitions
#' 
#' @details
#' Population Structure:
#' - Males: age, genetic trait, neutral trait, alive status
#' - Females: age, genetic trait, neutral trait, alive status, observed frequency,
#'   preferred trait type, chosen mate ID
#' 
#' Age Classes:
#' - Juvenile classes: non-reproductive, only last class can observe social information
#' - Adult classes: reproductive, can be observed by juveniles
#' 
#' Initial Conditions:
#' - Traits follow initial_trait_frequency parameter
#' - Age distribution follows geometric survival model
#' - Random initial mating establishes baseline for social learning
initialize_population <- function(parameters) {
  
  # ================================================================
  # SECTION 1: Initialize population matrices
  # ================================================================
  # Create separate matrices for males and females with different column requirements
  # Females need additional columns for social learning (observed, preferred, ID_father)
  
  males <- matrix(NA, nrow = parameters$K, ncol = 4)
  females <- matrix(NA, nrow = parameters$K, ncol = 7)
  colnames(males) <- c("age", "trait", "neutral", "alive")
  colnames(females) <- c("age", "trait", "neutral", "alive", "observed", "preferred", "ID_father")
  
  # ================================================================
  # SECTION 2: Define age-class behavioral characteristics
  # ================================================================
  # These vectors determine which individuals can reproduce and observe social information
  # This creates the generational structure necessary for cultural transmission
  
  # Reproductive capacity: juveniles cannot reproduce, adults can mate
  reproductive <- c(rep(FALSE, parameters$nb_class[1]), rep(TRUE, parameters$nb_class[2]))
  
  # Observing capacity: only last juvenile class can learn socially
  # This models cognitive development - very young juveniles lack learning ability
  # Last adult class doesn't observe (focused on reproduction, not learning)
  observing <- c(rep(FALSE, parameters$nb_class[1]-1), TRUE, rep(TRUE, parameters$nb_class[2]-1), FALSE)
  
  # ================================================================
  # SECTION 3: Calculate realistic age structure
  # ================================================================
  # Use geometric survival model to create natural age pyramid
  # Older individuals are progressively rarer due to cumulative mortality
  
  n <- sum(parameters$nb_class)
  pr <- with(parameters, survival^(1:n) / ((survival-survival^(n+1)) / (1-survival))) 
  
  # Random sex assignment with binomial sampling
  n_males <- rbinom(1, size = parameters$K, prob = 1/2)
  n_females <- parameters$K - n_males
  
  # Distribute individuals across age classes using multinomial sampling
  n_males <- as.numeric(rmultinom(1, size = n_males, prob = pr))
  n_females <- as.numeric(rmultinom(1, size = n_females, prob = pr))
  
  n_tot_females <- sum(n_females)
  n_tot_males <- sum(n_males)
  
  # ================================================================
  # SECTION 4: Initialize individual traits and characteristics
  # ================================================================
  
  # Initialize females with genetic traits and demographic characteristics
  if (n_tot_females > 0) {
    pos <- 1:n_tot_females
    females[pos, "age"] <- rep(1:n, times = n_females)
    
    # Genetic trait: subject to sexual selection and cultural transmission
    females[pos, "trait"] <- rbinom(n_tot_females, size = 1, prob = parameters$initial_trait_frequency)
    
    # Neutral trait: control for genetic drift, not subject to selection
    females[pos, "neutral"] <- rbinom(n_tot_females, size = 1, prob = parameters$initial_trait_frequency)
    
    females[pos, "alive"] <- 1
  }
  
  # Initialize males with identical genetic traits and demographic characteristics
  if (n_tot_males > 0) {
    pos <- 1:n_tot_males
    males[pos, "age"] <- rep(1:n, times = n_males)
    males[pos, "trait"] <- rbinom(n_tot_males, size = 1, prob = parameters$initial_trait_frequency)
    males[pos, "neutral"] <- rbinom(n_tot_males, size = 1, prob = parameters$initial_trait_frequency)
    males[pos, "alive"] <- 1
  }
  
  # ================================================================
  # SECTION 5: Establish initial mating pairs
  # ================================================================
  # Create initial mate choice pairs that juveniles can observe
  # This enables social learning from the first generation
  
  pos_f <- which(reproductive[females[, 1]]) # Reproductive females
  pos_m <- which(reproductive[males[, 1]])   # Reproductive males
  
  if (length(pos_f) > 0 & length(pos_m) > 0) {
    if (length(pos_m) == 1) {
      # Single male available: all females mate with him
      females[pos_f, "ID_father"] <- pos_m
    } else {
      # Multiple males: random assignment with replacement (polygyny possible)
      females[pos_f, "ID_father"] <- sample(pos_m, size = length(pos_f), replace = TRUE)
    }
  }
  
  return(list(males = males, females = females, reproductive = reproductive, observing = observing))
}


#### Cultural Model: Unconstrained Mate Choice ####

#' Cultural Model - Tradition Emergence Through Pure Conformist Social Learning
#' 
#' Implements mate choice where females have direct access to their preferred
#' male types, allowing pure cultural dynamics without resource constraints.
#' This serves as a baseline to isolate the effects of social transmission
#' on tradition formation.
#' 
#' @param tmax Maximum number of generations to simulate
#' @param parameters List of model parameters
#' @param return_tradition_only Logical; if TRUE, returns only tradition duration,
#'   if FALSE, returns complete temporal dynamics
#' @return Either tradition statistics or complete simulation results
#' 
#' @details
#' Key Features:
#' - Females encounter males of preferred type with probability c
#' - Mate choice independent of male type frequencies
#' - Pure cultural evolution without gene-culture coevolution
#' - Traditions require majority exaggeration in conformity curve
cultural_model <- function(tmax, parameters, return_tradition_only = TRUE) {
  
  # ================================================================
  # SECTION 1: Initialize simulation
  # ================================================================
  
  pop <- initialize_population(parameters)
  nb <- sum(parameters$nb_class) # Total age classes
  
  # Calculate initial demographic structure
  n_females <- n_females_old <- table(factor(pop$females[, "age"], levels = 1:nb))
  n_males <- n_males_old <- table(factor(pop$males[, "age"], levels = 1:nb))
  
  # ================================================================
  # SECTION 2: Prepare output storage based on analysis mode
  # ================================================================
  
  if (!return_tradition_only) {
    # Full temporal analysis mode: store complete evolutionary trajectories
    
    # Male trait evolution vectors
    sequence_trait_m <- sequence_neutral_m <- rep(NA, tmax + 1)
    sequence_trait_m[1] <- mean(pop$males[, "trait"], na.rm = TRUE)
    sequence_neutral_m[1] <- mean(pop$males[, "neutral"], na.rm = TRUE)
    
    # Female trait and behaviour evolution vectors  
    sequence_trait_f <- sequence_neutral_f <- rep(NA, tmax + 1)
    sequence_observed <- sequence_preferred <- sequence_choice <- rep(NA, tmax + 1)
    sequence_trait_f[1] <- mean(pop$females[, "trait"], na.rm = TRUE)
    sequence_neutral_f[1] <- mean(pop$females[, "neutral"], na.rm = TRUE)
    sequence_choice[1] <- mean(pop$males[pop$females[, "ID_father"], "trait"], na.rm = TRUE)
    sequence_observed[1] <- NA  # No observations at initialization
    sequence_preferred[1] <- NA # No learned preferences at initialization
    
    # Demographic tracking matrix (age class distributions)
    pop_size <- matrix(NA, nrow = tmax + 1, ncol = 2*nb)
    pop_size[1,] <- c(n_females, n_males)
    
  } else {
    # Tradition detection mode: minimal storage for tradition analysis
    mean_choice <- NA
    trad_duration <- NA
  }
  
  # Activity tracking vectors (always needed for final statistics)
  sequence_nb_demos <- sequence_nb_observations <- rep(NA, tmax + 1)
  sequence_nb_demos[1] <- NA
  sequence_nb_observations[1] <- NA
  
  # ================================================================
  # SECTION 3: Main simulation loop - generation by generation
  # ================================================================
  
  for(g in 1:tmax) {
    
    # ----------------------------------------------------------------
    # STEP 1: Identify demonstrator females from previous generation
    # ----------------------------------------------------------------
    # These are reproductive females who successfully mated and can be
    # observed by juvenile females for social learning
    pos_fR <- which(pop$reproductive[pop$females[, "age"]] & !is.na(pop$females[, "ID_father"]))
    
    # ----------------------------------------------------------------
    # STEP 2: Social Learning - Conformist Preference Acquisition
    # ----------------------------------------------------------------
    # Juvenile females observe matings and develop preferences 
    # for the observed majority choice)
    
    pos_fO <- which(pop$observing[pop$females[, "age"]]) # Observing females
    
    if (parameters$n_obs == 0) {
      # No social learning possible
      pop$females[pos_fO, "observed"] <- NA
      pop$females[pos_fO, "preferred"] <- rbinom(length(pos_fO), size = 1, prob = 1/2)
      sequence_nb_observations[g + 1] <- parameters$n_obs
      
    } else {
      # SOCIAL LEARNING MODEL: Conformist observation and preference formation
      
      # Determine what can be observed based on learning strategy
      if (parameters$female_strategy == "conformity") {
        # Observe actual mating outcomes from previous generation
        traits_observed_males <- pop$males[pop$females[pos_fR, "ID_father"], "trait"]
      }
      if (parameters$female_strategy == "koinophilia") {
        # Observe trait frequencies among available males
        traits_observed_males <- pop$males[which(pop$reproductive[pop$males[, "age"]]), "trait"]
      }
      
      # Handle limited observational capacity scenarios
      if (parameters$n_obs >= length(traits_observed_males)) {
        # LIMITED INFORMATION: All females observe same (limited) data
        sequence_nb_observations[g + 1] <- length(traits_observed_males)
        
        if (length(traits_observed_males) == 0) {
          # No demonstrations available: random preference assignment
          pop$females[pos_fO, "observed"] <- NA
          pop$females[pos_fO, "preferred"] <- rbinom(length(pos_fO), size = 1, prob = 1/2)
        } else {
          # Apply majority rule to observed frequency
          fr_obs <- mean(traits_observed_males)
          pop$females[pos_fO, "observed"] <- fr_obs
          
          if (fr_obs > 1/2) {
            pop$females[pos_fO, "preferred"] <- 1
          } else if (fr_obs == 1/2) {
            # Tie-breaking: random assignment when no clear majority
            pop$females[pos_fO,"preferred"] <- rbinom(length(pos_fO), size = 1, prob = 1/2)
          } else {
            pop$females[pos_fO,"preferred"] <- 0
          }
        }
      } else {
        # INDIVIDUAL SAMPLING: Each female gets unique sample via hypergeometric distribution
        # This creates individual variation in observations and preferences
        sequence_nb_observations[g + 1] <- parameters$n_obs
        
        nb_1 <- sum(traits_observed_males) # Type 1 males available to observe
        nb_0 <- length(traits_observed_males) - nb_1 # Type 0 males available
        
        # Hypergeometric sampling: each female samples without replacement
        nb_1_obs <- rhyper(length(pos_fO), nb_1, nb_0, parameters$n_obs)
        fr_obs <- nb_1_obs/parameters$n_obs
        
        pop$females[pos_fO, "observed"] <- fr_obs
        pop$females[pos_fO, "preferred"] <- ifelse(fr_obs > 1/2, 1, 0)
        
        # Handle individual ties in observations
        tied_observers <- fr_obs == 1/2
        if (any(tied_observers)) {
          pop$females[pos_fO[tied_observers], "preferred"] <- rbinom(sum(tied_observers), size = 1, prob = 1/2)
        }
      }
    }
    
    # ----------------------------------------------------------------
    # STEP 3: Demographic processes - Survival, reproduction, mutation
    # ----------------------------------------------------------------
    
    # SURVIVAL PHASE: Age-dependent mortality
    pos <- which(pop$females[, "age"] < nb) # Non-final age classes
    pop$females[pos, "alive"] <- with(parameters, rbinom(length(pos), size = 1, prob = survival))
    pop$females[which(pop$females[, "age"] == nb), "alive"] <- 0 # Final age class dies
    
    pos <- which(pop$males[, "age"] < nb)
    pop$males[pos, "alive"] <- with(parameters, rbinom(length(pos), size = 1, prob = survival))
    pop$males[which(pop$males[, "age"] == nb), "alive"] <- 0
    
    # POPULATION REGULATION: Maintain constant population size
    nb_newborns <- parameters$K - sum(pop$males[, "alive"], na.rm = TRUE) - sum(pop$females[, "alive"], na.rm = TRUE)
    
    # Safety checks for demographic balance
    if (nb_newborns < 0) stop("Number of newborns cannot be negative!")
    if (nb_newborns > parameters$K) stop("Number of newborns cannot exceed max population size!")
    
    # REPRODUCTION AND GENETIC TRANSMISSION
    if (length(pos_fR) > 0) { # If reproductive females exist
      if (nb_newborns > 0) { # And population has space for offspring
        
        # Select mothers: random sampling with replacement allows variance in reproductive success
        if (length(pos_fR) == 1) {
          who_is_the_mother <- rep(pos_fR, nb_newborns)
        } else {
          who_is_the_mother <- sample(pos_fR, size = nb_newborns, replace = TRUE)
        }
        who_is_the_father <- pop$females[who_is_the_mother, "ID_father"]
        
        # GENETIC INHERITANCE: Each trait inherited from random parent
        newborn_trait <- ifelse(runif(nb_newborns) < 0.5, 
                                pop$females[who_is_the_mother, "trait"], 
                                pop$males[who_is_the_father, "trait"])
        newborn_neutral <- ifelse(runif(nb_newborns) < 0.5, 
                                  pop$females[who_is_the_mother, "neutral"], 
                                  pop$males[who_is_the_father, "neutral"])
        
        # MUTATION PROCESS: Introduces genetic variation
        mute <- runif(nb_newborns) < parameters$mutation_rate
        if (any(mute)) newborn_trait[mute] <- rbinom(sum(mute), size = 1, prob = 1/2)
        mute <- runif(nb_newborns) < parameters$mutation_rate
        if (any(mute)) newborn_neutral[mute] <- rbinom(sum(mute), size = 1, prob = 1/2)
        
        # SEX DETERMINATION: Random assignment
        newborn_female <- runif(nb_newborns) < 0.5
      }
    }
    
    # ----------------------------------------------------------------
    # STEP 4: Aging and demographic transitions
    # ----------------------------------------------------------------
    
    # Remove dead individuals (set to NA)
    pop$males[pop$males[, "alive"] == 0, ] <- NA
    pop$females[pop$females[, "alive"] == 0, ] <- NA
    
    # Age advancement: individuals may advance to next age class
    pos_males_age <- which(!is.na(pop$males[, "age"]))
    pop$males[pos_males_age, "age"] <- ifelse(runif(length(pos_males_age)) > parameters$ageing, 
                                              pop$males[pos_males_age, "age"] + 1, 
                                              pop$males[pos_males_age, "age"])
    
    pos_females_age <- which(!is.na(pop$females[, "age"]))
    pop$females[pos_females_age, "age"] <- ifelse(runif(length(pos_females_age)) > parameters$ageing, 
                                                  pop$females[pos_females_age, "age"] + 1, 
                                                  pop$females[pos_females_age, "age"])
    
    # ----------------------------------------------------------------
    # STEP 5: Add newborns to population
    # ----------------------------------------------------------------
    
    if (length(pos_fR) > 0 && exists("nb_newborns") && nb_newborns > 0) {
      
      # Add female newborns
      if (exists("newborn_female") && length(newborn_female) > 0 && any(newborn_female)) {
        pos <- which(is.na(pop$females[, "age"]))
        pos_f <- pos[1:sum(newborn_female)]
        pop$females[pos_f, "trait"] <- newborn_trait[newborn_female]
        pop$females[pos_f, "neutral"] <- newborn_neutral[newborn_female]
        pop$females[pos_f, "age"] <- 1
        pop$females[pos_f, "alive"] <- 1
      }
      
      # Add male newborns
      if (exists("newborn_female") && length(newborn_female) > 0 && any(!newborn_female)) {
        pos <- which(is.na(pop$males[, "age"]))
        pos_m <- pos[1:sum(!newborn_female)]
        pop$males[pos_m, "trait"] <- newborn_trait[!newborn_female]
        pop$males[pos_m, "neutral"] <- newborn_neutral[!newborn_female]
        pop$males[pos_m, "age"] <- 1
        pop$males[pos_m, "alive"] <- 1
      }
    }
    
    # ----------------------------------------------------------------
    # STEP 6: Update reproductive status after demographic changes
    # ----------------------------------------------------------------
    
    pos_fR <- which(pop$reproductive[pop$females[, "age"]])
    pos_mR <- which(pop$reproductive[pop$males[, "age"]])
    
    # Check for population viability
    if (length(pos_fR) == 0) {
      warning(paste("No reproductive females at generation", g))
      nb_repro <- 0
      break
    }
    if (length(pos_mR) == 0) {
      warning(paste("No reproductive males at generation", g))
      nb_repro <- 0
      break
    }
    
    # ----------------------------------------------------------------
    # STEP 7: CULTURAL MODEL MATE CHOICE - Unconstrained Access
    # ----------------------------------------------------------------
    # Key difference from Gene-Culture Model: females can directly access
    # their preferred male type regardless of frequency in population
    
    pop$females[, "ID_father"] <- NA # Reset previous matings
    pos_f_non_mated <- which(pop$reproductive[pop$females[, "age"]] & is.na(pop$females[, "ID_father"]))
    pos_m_non_mated <- pos_mR
    
    # Categorize individuals by preference and trait
    pos_f_0 <- which(pop$females[pos_f_non_mated, "preferred"] == 0) # Females preferring type 0
    pos_f_1 <- which(pop$females[pos_f_non_mated, "preferred"] == 1) # Females preferring type 1
    pos_m_0 <- which(pop$males[pos_m_non_mated, "trait"] == 0)      # Type 0 males
    pos_m_1 <- which(pop$males[pos_m_non_mated, "trait"] == 1)      # Type 1 males
    
    # Initialize encounter tracking
    pos_m_met_0 <- pos_m_met_1 <- pos_m_error_0 <- pos_m_error_1 <- NULL
    
    # SCENARIO 1: Only type 1 males available
    if (length(pos_m_0) == 0 & length(pos_m_1) > 0) {
      if (parameters$c == 1) {
        # Perfect copying: only type 1 preferring females can mate
        if (length(pos_f_1) > 0) {
          if (length(pos_m_1) == 1) {
            pos_m_met_1 <- rep(pos_m_non_mated[pos_m_1], length(pos_f_1))
          } else {
            pos_m_met_1 <- sample(pos_m_non_mated[pos_m_1], length(pos_f_1), replace = TRUE)
          }
          pop$females[pos_f_non_mated[pos_f_1], "ID_father"] <- pos_m_met_1
        }
        pop$females[pos_f_non_mated[pos_f_0], "ID_father"] <- NA # Type 0 preferring females don't mate
      } else {
        # Imperfect copying: all females mate with available type 1 males
        if (length(pos_m_1) == 1) {
          pos_m_met <- rep(pos_m_non_mated[pos_m_1], length(pos_f_non_mated))
        } else {
          pos_m_met <- sample(pos_m_non_mated[pos_m_1], length(pos_f_non_mated), replace = TRUE)
        }
        pop$females[pos_f_non_mated, "ID_father"] <- pos_m_met
      }
    }
    
    # SCENARIO 2: Only type 0 males available (symmetric to scenario 1)
    if (length(pos_m_1) == 0 & length(pos_m_0) > 0) {
      if (parameters$c == 1) {
        if (length(pos_f_0) > 0) {
          if (length(pos_m_0) == 1) {
            pos_m_met_0 <- rep(pos_m_non_mated[pos_m_0], length(pos_f_0))
          } else {
            pos_m_met_0 <- sample(pos_m_non_mated[pos_m_0], length(pos_f_0), replace = TRUE)
          }
          pop$females[pos_f_non_mated[pos_f_0], "ID_father"] <- pos_m_met_0
        }
        pop$females[pos_f_non_mated[pos_f_1], "ID_father"] <- NA
      } else {
        if (length(pos_m_0) == 1) {
          pos_m_met <- rep(pos_m_non_mated[pos_m_0], length(pos_f_non_mated))
        } else {
          pos_m_met <- sample(pos_m_non_mated[pos_m_0], length(pos_f_non_mated), replace = TRUE)
        }
        pop$females[pos_f_non_mated, "ID_father"] <- pos_m_met
      }
    }
    
    # SCENARIO 3: Both male types available - full preference expression
    if (length(pos_m_0) > 0 & length(pos_m_1) > 0) {
      
      # Mate choice for females preferring type 0
      if (length(pos_f_0) > 0) {
        if (length(pos_m_0) == 1) {
          pos_m_met_0 <- rep(pos_m_non_mated[pos_m_0], length(pos_f_0))
        } else {
          pos_m_met_0 <- sample(pos_m_non_mated[pos_m_0], length(pos_f_0), replace = TRUE)
        }
        
        # Apply copying fidelity: probability c of accepting preferred type
        chosen_0 <- rbinom(length(pos_f_0), size = 1, prob = parameters$c)
        occurences_0 <- table(factor(chosen_0, levels = 0:1))
        
        # Females who accept their preferred type
        if (occurences_0["1"] > 0) {
          pop$females[pos_f_non_mated[pos_f_0[chosen_0 == 1]], "ID_father"] <- pos_m_met_0[chosen_0 == 1]
        }
        
        # Females who reject preferred type get alternative type
        if (occurences_0["0"] > 0) {
          if (length(pos_m_1) == 1) {
            pos_m_error_0 <- rep(pos_m_non_mated[pos_m_1], occurences_0["0"])
          } else {
            pos_m_error_0 <- sample(pos_m_non_mated[pos_m_1], occurences_0["0"], replace = TRUE)
          }
          pop$females[pos_f_non_mated[pos_f_0[chosen_0 == 0]], "ID_father"] <- pos_m_error_0
        }
      }
      
      # Mate choice for females preferring type 1 (symmetric logic)
      if (length(pos_f_1) > 0) {
        if (length(pos_m_1) == 1) {
          pos_m_met_1 <- rep(pos_m_non_mated[pos_m_1], length(pos_f_1))
        } else {
          pos_m_met_1 <- sample(pos_m_non_mated[pos_m_1], length(pos_f_1), replace = TRUE)
        }
        
        chosen_1 <- rbinom(length(pos_f_1), size = 1, prob = parameters$c)
        occurences_1 <- table(factor(chosen_1, levels = 0:1))
        
        if (occurences_1["1"] > 0) {
          pop$females[pos_f_non_mated[pos_f_1[chosen_1 == 1]], "ID_father"] <- pos_m_met_1[chosen_1 == 1]
        }
        
        if (occurences_1["0"] > 0) {
          if (length(pos_m_0) == 1) {
            pos_m_error_1 <- rep(pos_m_non_mated[pos_m_0], occurences_1["0"])
          } else {
            pos_m_error_1 <- sample(pos_m_non_mated[pos_m_0], occurences_1["0"], replace = TRUE)
          }
          pop$females[pos_f_non_mated[pos_f_1[chosen_1 == 0]], "ID_father"] <- pos_m_error_1
        }
      }
    }
    
    # ----------------------------------------------------------------
    # STEP 8: Record biological activity and update output
    # ----------------------------------------------------------------
    
    nb_repro <- length(which(pop$reproductive[pop$females[, "age"]] & !is.na(pop$females[, "ID_father"])))
    sequence_nb_demos[g + 1] <- nb_repro
    
    if (!return_tradition_only) {
      # Full temporal analysis: record all evolutionary variables
      sequence_trait_m[g + 1] <- mean(pop$males[, "trait"], na.rm = TRUE)
      sequence_neutral_m[g + 1] <- mean(pop$males[, "neutral"], na.rm = TRUE)
      sequence_trait_f[g + 1] <- mean(pop$females[, "trait"], na.rm = TRUE)
      sequence_neutral_f[g + 1] <- mean(pop$females[, "neutral"], na.rm = TRUE)
      sequence_observed[g + 1] <- mean(pop$females[, "observed"], na.rm = TRUE)
      sequence_preferred[g + 1] <- mean(pop$females[, "preferred"], na.rm = TRUE)
      sequence_choice[g + 1] <- mean(pop$males[pop$females[, "ID_father"], "trait"], na.rm = TRUE)
      
      # Update demographic tracking
      n_females <- table(factor(pop$females[, "age"], levels = 1:nb))
      n_males  <- table(factor(pop$males[, "age"], levels = 1:nb))      
      pop_size[g + 1, ] <- c(n_females, n_males)
      
    } else {
      # Tradition detection mode: monitor for tradition termination
      if (g == parameters$n_steps_removed) {
        # Initialize tradition tracking
        mean_choice <- mean(pop$males[pop$females[, "ID_father"], "trait"], na.rm = TRUE)
        trad_duration <- 0
      }
      
      if (g > parameters$n_steps_removed) {
        # Active tradition monitoring: detect changes in majority choice
        mean_choice_old <- mean_choice # Previous generation's majority
        mean_choice <- mean(pop$males[pop$females[, "ID_father"], "trait"], na.rm = TRUE)
        
        # XOR logic detects tradition termination (majority switch)
        if (xor(mean_choice > 0.5, mean_choice_old > 0.5)) {
          return(c(last_time_step = g, mean_choice_last_step = mean_choice_old, 
                   n_repro_last_step = nb_repro, mean_demos = mean(sequence_nb_demos, na.rm = T),
                   mean_n_obs = mean(sequence_nb_observations, na.rm = T),
                   tradition_duration = trad_duration))
        }
        trad_duration <- trad_duration + 1
      }
    }
  } # End of main simulation loop
  
  # ================================================================
  # SECTION 4: Return results based on analysis mode
  # ================================================================
  
  if (!return_tradition_only) {
    # Complete evolutionary analysis results
    results <- list(population = pop, 
                    male_trait = sequence_trait_m, male_neutral = sequence_neutral_m,
                    female_trait = sequence_trait_f, female_neutral = sequence_neutral_f,
                    observed = sequence_observed, preference = sequence_preferred,
                    choice = sequence_choice, effective_nobs = sequence_nb_observations,
                    n_demos = sequence_nb_demos,
                    size = pop_size)
    return(results)
  } else {
    # Tradition analysis results
    return(c(last_time_step = g, mean_choice_last_step = mean_choice, 
             n_repro_last_step = nb_repro, mean_demos = mean(sequence_nb_demos, na.rm = T),
             mean_n_obs = mean(sequence_nb_observations, na.rm = T),
             tradition_duration = trad_duration))
  }
}


#### Gene-Culture Model: Constrained Mate Choice ####

#' Gene-Culture Coevolution Model - Constrained Encounter-Based Mate Choice
#' 
#' Implements mate choice where female encounters depend on male type frequencies,
#' creating gene-culture coevolution between socially transmitted preferences
#' and genetically inherited traits. This enables traditions to emerge even
#' without majority exaggeration in conformity curves.
#' 
#' @param tmax Maximum number of generations to simulate
#' @param parameters List of model parameters
#' @param return_tradition_only Logical; if TRUE, returns only tradition duration
#' @return Either tradition statistics or complete simulation results
#' 
#' @details
#' Key Differences from Cultural Model:
#' - Sequential encounters based on male type frequencies
#' - Males have limited mating capacity (n_matings parameter)
#' - Gene-culture coevolution creates reinforcing dynamics
#' - Traditions can emerge without majority exaggeration
#' - Implements "sexy son" effect through frequency-dependent mating
geneculture_model <- function(tmax, parameters, return_tradition_only = TRUE) {
  
  # ================================================================
  # SECTION 1: Initialize simulation (identical to Cultural Model)
  # ================================================================
  
  pop <- initialize_population(parameters)
  nb <- sum(parameters$nb_class)
  
  n_females <- n_females_old <- table(factor(pop$females[, "age"], levels = 1:nb))
  n_males <- n_males_old <- table(factor(pop$males[, "age"], levels = 1:nb))
  
  # ================================================================
  # SECTION 2: Prepare output storage (identical structure)
  # ================================================================
  
  if (!return_tradition_only) {
    sequence_trait_m <- sequence_neutral_m <- rep(NA, tmax + 1)
    sequence_trait_m[1] <- mean(pop$males[, "trait"], na.rm = TRUE)
    sequence_neutral_m[1] <- mean(pop$males[, "neutral"], na.rm = TRUE)
    
    sequence_trait_f <- sequence_neutral_f <- rep(NA, tmax + 1)
    sequence_observed <- sequence_preferred <- sequence_choice <- rep(NA, tmax + 1)
    sequence_trait_f[1] <- mean(pop$females[, "trait"], na.rm = TRUE)
    sequence_neutral_f[1] <- mean(pop$females[, "neutral"], na.rm = TRUE)
    sequence_choice[1] <- mean(pop$males[pop$females[, "ID_father"], "trait"], na.rm = TRUE)
    sequence_observed[1] <- NA
    sequence_preferred[1] <- NA
    
    pop_size <- matrix(NA, nrow = tmax + 1, ncol = 2*nb)
    pop_size[1,] <- c(n_females, n_males)
  } else {
    mean_choice <- NA
    trad_duration <- NA
  }
  
  sequence_nb_demos <- sequence_nb_observations <- rep(NA, tmax + 1)
  sequence_nb_demos[1] <- NA
  sequence_nb_observations[1] <- NA
  
  # ================================================================
  # SECTION 3: Main simulation loop
  # ================================================================
  
  for(g in 1:tmax) {
    
    # NB: STEPS 1-6 = Identical to Cultural Model
    # ----------------------------------------------------------------
    # STEP 1: Identify demonstrator females from previous generation
    # ----------------------------------------------------------------
    # These are reproductive females who successfully mated and can be
    # observed by juvenile females for social learning
    pos_fR <- which(pop$reproductive[pop$females[, "age"]] & !is.na(pop$females[, "ID_father"]))
    
    # ----------------------------------------------------------------
    # STEP 2: Social Learning - Conformist Preference Acquisition
    # ----------------------------------------------------------------
    # Juvenile females observe matings and develop preferences 
    # for the observed majority choice)
    
    pos_fO <- which(pop$observing[pop$females[, "age"]]) # Observing females
    
    if (parameters$n_obs == 0) {
      # No social learning possible
      pop$females[pos_fO, "observed"] <- NA
      pop$females[pos_fO, "preferred"] <- rbinom(length(pos_fO), size = 1, prob = 1/2)
      sequence_nb_observations[g + 1] <- parameters$n_obs
      
    } else {
      # SOCIAL LEARNING MODEL: Conformist observation and preference formation
      
      # Determine what can be observed based on learning strategy
      if (parameters$female_strategy == "conformity") {
        # Observe actual mating outcomes from previous generation
        traits_observed_males <- pop$males[pop$females[pos_fR, "ID_father"], "trait"]
      }
      if (parameters$female_strategy == "koinophilia") {
        # Observe trait frequencies among available males
        traits_observed_males <- pop$males[which(pop$reproductive[pop$males[, "age"]]), "trait"]
      }
      
      # Handle limited observational capacity scenarios
      if (parameters$n_obs >= length(traits_observed_males)) {
        # LIMITED INFORMATION: All females observe same (limited) data
        sequence_nb_observations[g + 1] <- length(traits_observed_males)
        
        if (length(traits_observed_males) == 0) {
          # No demonstrations available: random preference assignment
          pop$females[pos_fO, "observed"] <- NA
          pop$females[pos_fO, "preferred"] <- rbinom(length(pos_fO), size = 1, prob = 1/2)
        } else {
          # Apply majority rule to observed frequency
          fr_obs <- mean(traits_observed_males)
          pop$females[pos_fO, "observed"] <- fr_obs
          
          if (fr_obs > 1/2) {
            pop$females[pos_fO, "preferred"] <- 1
          } else if (fr_obs == 1/2) {
            # Tie-breaking: random assignment when no clear majority
            pop$females[pos_fO,"preferred"] <- rbinom(length(pos_fO), size = 1, prob = 1/2)
          } else {
            pop$females[pos_fO,"preferred"] <- 0
          }
        }
      } else {
        # INDIVIDUAL SAMPLING: Each female gets unique sample via hypergeometric distribution
        # This creates individual variation in observations and preferences
        sequence_nb_observations[g + 1] <- parameters$n_obs
        
        nb_1 <- sum(traits_observed_males) # Type 1 males available to observe
        nb_0 <- length(traits_observed_males) - nb_1 # Type 0 males available
        
        # Hypergeometric sampling: each female samples without replacement
        nb_1_obs <- rhyper(length(pos_fO), nb_1, nb_0, parameters$n_obs)
        fr_obs <- nb_1_obs/parameters$n_obs
        
        pop$females[pos_fO, "observed"] <- fr_obs
        pop$females[pos_fO, "preferred"] <- ifelse(fr_obs > 1/2, 1, 0)
        
        # Handle individual ties in observations
        tied_observers <- fr_obs == 1/2
        if (any(tied_observers)) {
          pop$females[pos_fO[tied_observers], "preferred"] <- rbinom(sum(tied_observers), size = 1, prob = 1/2)
        }
      }
    }
    
    # ----------------------------------------------------------------
    # STEP 3: Demographic processes - Survival, reproduction, mutation
    # ----------------------------------------------------------------
    
    # SURVIVAL PHASE: Age-dependent mortality
    pos <- which(pop$females[, "age"] < nb) # Non-final age classes
    pop$females[pos, "alive"] <- with(parameters, rbinom(length(pos), size = 1, prob = survival))
    pop$females[which(pop$females[, "age"] == nb), "alive"] <- 0 # Final age class dies
    
    pos <- which(pop$males[, "age"] < nb)
    pop$males[pos, "alive"] <- with(parameters, rbinom(length(pos), size = 1, prob = survival))
    pop$males[which(pop$males[, "age"] == nb), "alive"] <- 0
    
    # POPULATION REGULATION: Maintain constant population size
    nb_newborns <- parameters$K - sum(pop$males[, "alive"], na.rm = TRUE) - sum(pop$females[, "alive"], na.rm = TRUE)
    
    # Safety checks for demographic balance
    if (nb_newborns < 0) stop("Number of newborns cannot be negative!")
    if (nb_newborns > parameters$K) stop("Number of newborns cannot exceed max population size!")
    
    # REPRODUCTION AND GENETIC TRANSMISSION
    if (length(pos_fR) > 0) { # If reproductive females exist
      if (nb_newborns > 0) { # And population has space for offspring
        
        # Select mothers: random sampling with replacement allows variance in reproductive success
        if (length(pos_fR) == 1) {
          who_is_the_mother <- rep(pos_fR, nb_newborns)
        } else {
          who_is_the_mother <- sample(pos_fR, size = nb_newborns, replace = TRUE)
        }
        who_is_the_father <- pop$females[who_is_the_mother, "ID_father"]
        
        # GENETIC INHERITANCE: Each trait inherited from random parent
        newborn_trait <- ifelse(runif(nb_newborns) < 0.5, 
                                pop$females[who_is_the_mother, "trait"], 
                                pop$males[who_is_the_father, "trait"])
        newborn_neutral <- ifelse(runif(nb_newborns) < 0.5, 
                                  pop$females[who_is_the_mother, "neutral"], 
                                  pop$males[who_is_the_father, "neutral"])
        
        # MUTATION PROCESS: Introduces genetic variation
        mute <- runif(nb_newborns) < parameters$mutation_rate
        if (any(mute)) newborn_trait[mute] <- rbinom(sum(mute), size = 1, prob = 1/2)
        mute <- runif(nb_newborns) < parameters$mutation_rate
        if (any(mute)) newborn_neutral[mute] <- rbinom(sum(mute), size = 1, prob = 1/2)
        
        # SEX DETERMINATION: Random assignment
        newborn_female <- runif(nb_newborns) < 0.5
      }
    }
    
    # ----------------------------------------------------------------
    # STEP 4: Aging and demographic transitions
    # ----------------------------------------------------------------
    
    # Remove dead individuals (set to NA)
    pop$males[pop$males[, "alive"] == 0, ] <- NA
    pop$females[pop$females[, "alive"] == 0, ] <- NA
    
    # Age advancement: individuals may advance to next age class
    pos_males_age <- which(!is.na(pop$males[, "age"]))
    pop$males[pos_males_age, "age"] <- ifelse(runif(length(pos_males_age)) > parameters$ageing, 
                                              pop$males[pos_males_age, "age"] + 1, 
                                              pop$males[pos_males_age, "age"])
    
    pos_females_age <- which(!is.na(pop$females[, "age"]))
    pop$females[pos_females_age, "age"] <- ifelse(runif(length(pos_females_age)) > parameters$ageing, 
                                                  pop$females[pos_females_age, "age"] + 1, 
                                                  pop$females[pos_females_age, "age"])
    
    # ----------------------------------------------------------------
    # STEP 5: Add newborns to population
    # ----------------------------------------------------------------
    
    if (length(pos_fR) > 0 && exists("nb_newborns") && nb_newborns > 0) {
      
      # Add female newborns
      if (exists("newborn_female") && length(newborn_female) > 0 && any(newborn_female)) {
        pos <- which(is.na(pop$females[, "age"]))
        pos_f <- pos[1:sum(newborn_female)]
        pop$females[pos_f, "trait"] <- newborn_trait[newborn_female]
        pop$females[pos_f, "neutral"] <- newborn_neutral[newborn_female]
        pop$females[pos_f, "age"] <- 1
        pop$females[pos_f, "alive"] <- 1
      }
      
      # Add male newborns
      if (exists("newborn_female") && length(newborn_female) > 0 && any(!newborn_female)) {
        pos <- which(is.na(pop$males[, "age"]))
        pos_m <- pos[1:sum(!newborn_female)]
        pop$males[pos_m, "trait"] <- newborn_trait[!newborn_female]
        pop$males[pos_m, "neutral"] <- newborn_neutral[!newborn_female]
        pop$males[pos_m, "age"] <- 1
        pop$males[pos_m, "alive"] <- 1
      }
    }
    
    # ----------------------------------------------------------------
    # STEP 6: Update reproductive status after demographic changes
    # ----------------------------------------------------------------
    
    pos_fR <- which(pop$reproductive[pop$females[, "age"]])
    pos_mR <- which(pop$reproductive[pop$males[, "age"]])
    
    # Check for population viability
    if (length(pos_fR) == 0) {
      warning(paste("No reproductive females at generation", g))
      nb_repro <- 0
      break
    }
    if (length(pos_mR) == 0) {
      warning(paste("No reproductive males at generation", g))
      nb_repro <- 0
      break
    }
    
    # ----------------------------------------------------------------
    # STEP 7: GENE-CULTURE MODEL MATE CHOICE - Constrained Encounters
    # ----------------------------------------------------------------
    # KEY DIFFERENCE with Cultural Model: This is where gene-culture coevolution occurs
    # Females encounter males based on their frequencies in the population
    
    pop$females[, "ID_father"] <- NA # Reset mating assignments
    pos_f_non_mated <- which(pop$reproductive[pop$females[, "age"]] & is.na(pop$females[, "ID_father"]))
    pos_m_non_mated <- pos_mR
    
    # MALE MATING BUDGET: Each male can mate maximum n_matings times
    # This creates male-male competition and realistic reproductive constraints
    nb_matings_males <- rep(parameters$n_matings, length(pos_mR))
    names(nb_matings_males) <- pos_mR
    
    # SEQUENTIAL ENCOUNTER LOOP: Continue until all females mated or no males available
    while (length(pos_f_non_mated) > 0 && length(pos_m_non_mated) > 0) {
      
      # Safety check: prevent infinite loops with perfect copying and impossible preferences
      if (parameters$c == 1 & !any(pop$females[pos_f_non_mated, "preferred"] %in% pop$males[pos_m_non_mated, "trait"])) {  
        break 
      }
      
      # ----------------------------------------------------------------
      # ENCOUNTER ORGANIZATION: Determine who meets whom in this round
      # ----------------------------------------------------------------
      # This creates realistic constraints where not all individuals
      # can interact simultaneously
      
      pos_f_met <- pos_m_met <- NULL
      
      if (length(pos_f_non_mated) > length(pos_m_non_mated)) {
        # More females than available males: random female selection
        pos_f_met <- sample(pos_f_non_mated, length(pos_m_non_mated), replace = FALSE)
        pos_m_met <- pos_m_non_mated
      } else if (length(pos_m_non_mated) > length(pos_f_non_mated)) {
        # More males than unmated females: random male selection
        pos_m_met <- sample(pos_m_non_mated, length(pos_f_non_mated), replace = FALSE)
        pos_f_met <- pos_f_non_mated
      } else if (length(pos_m_non_mated) == length(pos_f_non_mated) && length(pos_m_non_mated) > 1) {
        # Equal numbers: shuffle male order for random pairing
        pos_m_met <- sample(pos_m_non_mated)
        pos_f_met <- pos_f_non_mated
      } else if (length(pos_m_non_mated) == length(pos_f_non_mated) && length(pos_m_non_mated) == 1) {
        # Final pair: last female meets last available male
        pos_m_met <- pos_m_non_mated
        pos_f_met <- pos_f_non_mated
      }
      
      # ----------------------------------------------------------------
      # PREFERENCE-BASED ACCEPTANCE DECISIONS
      # ----------------------------------------------------------------
      # Core gene-culture coevolution mechanism: acceptance probability
      # depends on match between female preference and encountered male trait
      
      pr <- ifelse(pop$males[pos_m_met, "trait"] == pop$females[pos_f_met, "preferred"], 
                   parameters$c,        # High acceptance for preferred type
                   1 - parameters$c)    # Low acceptance for non-preferred type
      
      # Validation checks for data integrity
      if(any(is.na(pop$males[pos_m_met, "trait"]))) stop("NA in encountered male traits")
      if (any(pr < 0) | any(pr > 1)) stop("Invalid acceptance probability!")
      
      # Execute acceptance decisions
      chosen <- rbinom(length(pos_f_met), size = 1, prob = pr)
      
      # ----------------------------------------------------------------
      # UPDATE POPULATION STATE AFTER SUCCESSFUL MATINGS
      # ----------------------------------------------------------------
      
      if (any(chosen == 1)) {
        # Record successful matings
        pop$females[pos_f_met[chosen == 1], "ID_father"] <- pos_m_met[chosen == 1]
        
        # Update list of unmated females
        pos_f_non_mated <- pos_fR[is.na(pop$females[pos_fR, "ID_father"])]
        
        # MALE MATING BUDGET MANAGEMENT: Reduce remaining capacity
        successful_males <- pos_m_met[chosen == 1]
        nb_matings_males[as.character(successful_males)] <- nb_matings_males[as.character(successful_males)] - 1
        
        # Remove males who have exhausted their mating capacity
        pos_m_non_mated <- pos_m_non_mated[nb_matings_males > 0]
        nb_matings_males <- nb_matings_males[nb_matings_males > 0]
      }
    } # End of sequential encounter loop
    
    # ----------------------------------------------------------------
    # STEP 8: Record activity and update output (identical to Cultural Model)
    # ----------------------------------------------------------------
    
    nb_repro <- length(which(pop$reproductive[pop$females[, "age"]] & !is.na(pop$females[, "ID_father"])))
    sequence_nb_demos[g + 1] <- nb_repro
    
    if (!return_tradition_only) {
      # Complete temporal tracking
      sequence_trait_m[g + 1] <- mean(pop$males[, "trait"], na.rm = TRUE)
      sequence_neutral_m[g + 1] <- mean(pop$males[, "neutral"], na.rm = TRUE)
      sequence_trait_f[g + 1] <- mean(pop$females[, "trait"], na.rm = TRUE)
      sequence_neutral_f[g + 1] <- mean(pop$females[, "neutral"], na.rm = TRUE)
      sequence_observed[g + 1] <- mean(pop$females[, "observed"], na.rm = TRUE)
      sequence_preferred[g + 1] <- mean(pop$females[, "preferred"], na.rm = TRUE)
      sequence_choice[g + 1] <- mean(pop$males[pop$females[, "ID_father"], "trait"], na.rm = TRUE)
      
      n_females <- table(factor(pop$females[, "age"], levels = 1:nb))
      n_males  <- table(factor(pop$males[, "age"], levels = 1:nb))      
      pop_size[g + 1, ] <- c(n_females, n_males)
      
    } else {
      # Tradition detection
      if (g == parameters$n_steps_removed) {
        mean_choice <- mean(pop$males[pop$females[, "ID_father"], "trait"], na.rm = TRUE)
        trad_duration <- 0
      }
      if (g > parameters$n_steps_removed) {
        mean_choice_old <- mean_choice
        mean_choice <- mean(pop$males[pop$females[, "ID_father"], "trait"], na.rm = TRUE)
        
        if (xor(mean_choice > 0.5, mean_choice_old > 0.5)) {
          return(c(last_time_step = g, mean_choice_last_step = mean_choice_old, 
                   n_repro_last_step = nb_repro, mean_demos = mean(sequence_nb_demos, na.rm = T),
                   mean_n_obs = mean(sequence_nb_observations, na.rm = T),
                   tradition_duration = trad_duration))
        }
        trad_duration <- trad_duration + 1
      }
    }
  } # End of main simulation loop
  
  # ================================================================
  # SECTION 4: Return results (identical to Cultural Model)
  # ================================================================
  
  if (!return_tradition_only) {
    results <- list(population = pop, 
                    male_trait = sequence_trait_m, male_neutral = sequence_neutral_m,
                    female_trait = sequence_trait_f, female_neutral = sequence_neutral_f,
                    observed = sequence_observed, preference = sequence_preferred,
                    choice = sequence_choice, effective_nobs = sequence_nb_observations,
                    n_demos = sequence_nb_demos,
                    size = pop_size)
    return(results)
  } else {
    return(c(last_time_step = g, mean_choice_last_step = mean_choice, 
             n_repro_last_step = nb_repro, mean_demos = mean(sequence_nb_demos, na.rm = T),
             mean_n_obs = mean(sequence_nb_observations, na.rm = T),
             tradition_duration = trad_duration))
  }
}


#### Parameter Exploration Function ####

#' Systematic Parameter Space Exploration for Tradition Studies
#' 
#' Executes comprehensive parameter sweeps across social learning parameters
#' (n_obs, c) and genetic transmission fidelity (mutation_rate) -- or any other 
#' parameter -- to identify conditions promoting tradition emergence. 
#' Designed specifically for return_tradition_only = TRUE mode.
#' 
#' @param dynamics_type Function to execute (cultural_model or geneculture_model)
#' @param param_grid Data.frame with parameter combinations to explore
#' @param parameters Base parameter list to be modified
#' @param file_name Output file for results storage
#' @param n_replicates Number of replicates per parameter combination
#' @param tmax Maximum generations per simulation
#' 
#' @details
#' Output Structure:
#' - simu: unique simulation identifier
#' - parameter columns: n_obs, mutation_rate, c values
#' - replicate: replicate number within parameter set
#' - tradition_duration: primary outcome variable
#' - Additional statistics: last_time_step, mean_choice_last_step, etc.
#' 
#' Error Handling:
#' - Failed simulations recorded as NA values
#' - Try-catch prevents single failures from stopping entire exploration
#' 
#' Performance Considerations:
#' - Large parameter grids may require substantial computational time
parameters_exploration <- function(dynamics_type, param_grid, parameters, file_name, n_replicates = 10, tmax = 200) {
  
  # ================================================================
  # SECTION 1: Initialize exploration tracking
  # ================================================================
  
  simu_index <- 1 # Global simulation counter for unique identification
  
  # ================================================================
  # SECTION 2: Main exploration loop
  # ================================================================
  
  for (i in 1:nrow(param_grid)) {
    
    # ----------------------------------------------------------------
    # SUBSECTION 2A: Configure parameters for current combination
    # ----------------------------------------------------------------
    # Update base parameter set with current grid values
    # This approach preserves all other parameters while varying only those in the grid
    parameters[colnames(param_grid)] <- param_grid[i,]
    
    # ----------------------------------------------------------------
    # SUBSECTION 2B: Execute replicate simulations
    # ----------------------------------------------------------------
    # Multiple replicates capture stochastic variation in tradition emergence
    
    for (r in 1:n_replicates) {
      
      # Execute simulation with error handling
      simu <- try(dynamics_type(tmax, parameters, return_tradition_only = TRUE))
      
      # ----------------------------------------------------------------
      # SUBSECTION 2C: Initialize output file on first simulation
      # ----------------------------------------------------------------
      
      if (simu_index == 1) {
        # Create header row with all relevant columns
        cat("", "simu", colnames(param_grid), "replicate", names(simu), 
            sep = ";", file = file_name, append = FALSE)
        cat("\n", file = file_name, append = TRUE)
      }
      
      # ----------------------------------------------------------------
      # SUBSECTION 2D: Handle simulation outcomes and record results
      # ----------------------------------------------------------------
      
      if ("try-error" %in% class(simu)) {
        # SIMULATION FAILURE: Record attempt with NA values
        # This maintains complete record of parameter space exploration
        cat(simu_index, i, unlist(param_grid[i,]), r, rep(NA, length(names(simu))), 
            sep = ";", file = file_name, append = TRUE)
        cat("\n", file = file_name, append = TRUE)
      } else {
        # SIMULATION SUCCESS: Record complete results
        cat(simu_index, i, unlist(param_grid[i,]), r, simu, 
            sep = ";", file = file_name, append = TRUE)
        cat("\n", file = file_name, append = TRUE)
      }
      
      simu_index <- simu_index + 1
    }
  }
  
  # Function completes silently, having written all results to file
  # Post-processing and analysis performed by separate scripts
}


# ================================================================================
# THE PARAMETERS AND PARAMETERS VARIATION
# ================================================================================

# Global parameters
# NB: the parameters mutation_rate, n_obs and c are varied in param_variation 
p <- list(
  nb_class = c(1,1),                    # Age structure: 1 juvenile + 1 adult class
  K = 1000,                             # Population carrying capacity
  female_strategy = "conformity",       # Learning strategy: "conformity" or "koinophilia"
  survival = 0.9,                       # Baseline survival probability
  ageing = 0.9,                         # Probability of advancing to next age class
  initial_trait_frequency = 0.5,        # Initial frequency of trait 1 in population
  n_matings = 100,                      # Maximum matings per male (gene-culture model only)
  mutation_rate = 0.1,                  # Genetic mutation rate on male trait per generation
  n_obs = 10,                           # Number of matings observed by learning females
  c = 0.7,                              # Copying fidelity (probability of following preference)
  n_steps_removed = 100                 # Time steps to exclude before tradition detection
)

# Number of generations in simulation
n_steps <- 1000

# Define parameter variations for systematic exploration
param_variation <- expand.grid(
  n_obs = c(seq(1, 30, 2), 224),                       # Different observation capacities
  mutation_rate = c(1, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05), # Different genetic transmission fidelities
  c = seq(0.5, 1, 0.05)                                # Different copying fidelities
)

# Null model parameters (no social learning)
param_variation_null <- expand.grid(
  n_obs = c(seq(1, 30, 2), 224),                       
  mutation_rate = c(1, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05), 
  c = 0.5                                              # Random mate choice
)

# Number of replicates per parameter combination
n_replicates <- 100

