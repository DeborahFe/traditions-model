#### parameters ####

## list of global parameters
p <- list(
  nb_class = c(1,1), ## i classes of juveniles; non reproductive; only the last observes
  ## j classes of adults; reproductive; only the last does not observe
  K = 200, ## pop size
  female_strategy = "conformity", ## alternatives: "koinophilia" or "conformity"
  survival = 0.9, ## must be < 1
  # ageing = 0.1, ## probability to stay in the same age class from one time step to the next
  initial_trait_frequency = 0.5, # initial probability that males have trait 1 (and neutral trait 1)
  # initial_choice_frequency = 0.1, # initial probability that females have chosen trait 1
  n_matings = 10, ## maximum number of matings for males at every time step, must be > 0
  heritability = 1, ## with mutation rate = 1 - heritability
  n_obs = 10, ## maximal number of couples observed per female
  c = 0.6, ## copying rate
  n_steps_removed = 50) ## must be > 0

## dataframe for parameters variation
param_variation <- expand.grid(n_obs = c(4, 6, 8, 10, 15),
                               heritability = 1,
                               c = seq(0.6, 0.8, 0.05))
param_variation_null <- expand.grid(n_obs = 0,
                               heritability = 1,
                               c = seq(0.6, 0.8, 0.05))

## number of replicates per combination of parameters
n_replicates <- 100

## total number of simulations to run
n_simus <- nrow(param_variation) * n_replicates

## number of time steps
n_steps <- 200


#### functions to run simulations ####

# parameters <- p

init_pop <- function(parameters) {
  males <- matrix(NA, nrow = parameters$K, ncol = 4) ## age / trait / neutral trait / alive / nb of matings    
  females <- matrix(NA, nrow = parameters$K, ncol = 7) ## age / trait / neutral trait / alive / observed / preferred / ID male chosen
  colnames(males) <- c("age", "trait", "neutral", "alive")
  colnames(females) <- c("age", "trait", "neutral", "alive", "observed", "preferred", "ID_father")
  ## description of classes
  reproductive = c(rep(FALSE, parameters$nb_class[1]), rep(TRUE, parameters$nb_class[2]))
  observing = c(rep(FALSE, parameters$nb_class[1]-1), TRUE, rep(TRUE, parameters$nb_class[2]-1), FALSE)
  ## initial repartition among age classes
  ## geometric decrease in frequency, based on survival chances
  n <- with(parameters,sum(nb_class))
  pr <- with(parameters, survival^(1:n) / ((survival-survival^(n+1)) / (1-survival))) 
  ## randomly draw number of males
  n_males <- rbinom(1, size = parameters$K, prob = 1/2)
  n_females <- parameters$K - n_males
  ## random repartition of individuals among age classes
  n_males <- as.numeric(rmultinom(1, size = n_males, prob = pr)) # number of males per age class
  n_females <- as.numeric(rmultinom(1, size = n_females, prob = pr)) # number of females per age class
  
  n_tot_females <- sum(n_females) # total number of females
  n_tot_males <- sum(n_males) # total number of males
  
  ## Initialization of age, trait, neutral trait and "who is alive" for females and males
  if (n_tot_females > 0) {
    pos <- 1 : n_tot_females
    females[pos, "age"] <- rep(1:n, times = n_females) # initialization of age in the female matrix
    females[pos, "trait"] <- rbinom(n_tot_females, size = 1, prob = parameters$initial_trait_frequency) # initialization of trait
    females[pos, "neutral"] <- rbinom(n_tot_females, size = 1, prob = parameters$initial_trait_frequency) # initialization of neutral trait
    females[pos, "alive"] <- 1 # initialization of "who is alive"
  }
  if (n_tot_males > 0) {
    pos <- 1 : n_tot_males 
    males[pos, "age"] <- rep(1:n, times = n_males) # initialization of age in the male matrix
    males[pos, "trait"] <- rbinom(n_tot_males, size = 1, prob = parameters$initial_trait_frequency) # initialization of trait
    males[pos, "neutral"] <- rbinom(n_tot_males, size = 1, prob = parameters$initial_trait_frequency) # initialization of neutral trait
    males[pos, "alive"] <- 1 # initialization of "who is alive"
  }
  
  ## Initialization of male chosen for reproductive females
  pos_f <- which(reproductive[females[, 1]]) # positions of reproductive females
  pos_m <- which(reproductive[males[, 1]]) # positions of reproductive males
  if (length(pos_f) > 0 & length(pos_m) > 0) {
    if (length(pos_m) == 1) {
      females[pos_f, "ID_father"] <- pos_m
    } else {
      females[pos_f, "ID_father"] <- sample(pos_m, size = length(pos_f), replace = TRUE)
    }
  }
  
  return(list(males = males, females = females, reproductive = reproductive, observing = observing))
}
# 
# parameters <- p
# tmax <- 10
# return_tradition_only = FALSE

dynamics <- function(tmax, parameters, return_tradition_only = FALSE) { # if return_tradition_only = FALSE, store raw temporal results, if not, store only one tradition (duration) per simulation
  ## initialization
  pop <- init_pop(parameters)
  
  nb <- sum(parameters$nb_class) # number of age classes
  n_females <- n_females_old <- table(factor(pop$females[, "age"], levels = 1 : nb)) # number of females per age class
  n_males <- n_males_old <- table(factor(pop$males[, "age"], levels = 1 : nb)) # number of males per age class
  
  ## prepare vectors and matrices to store output (only if we store raw results)
  if (!return_tradition_only) {
    # for males
    sequence_trait_m <- sequence_neutral_m <- rep(NA, tmax + 1)
    sequence_trait_m[1] <- mean(pop$males[, "trait"], na.rm = TRUE)
    sequence_neutral_m[1] <- mean(pop$males[, "neutral"], na.rm = TRUE)
    # for females
    sequence_trait_f <- sequence_neutral_f <- rep(NA, tmax + 1)
    sequence_observed <- sequence_preferred <- sequence_choice <- rep(NA, tmax + 1)
    sequence_trait_f[1] <- mean(pop$females[, "trait"], na.rm = TRUE)
    sequence_neutral_f[1] <- mean(pop$females[, "neutral"], na.rm = TRUE)
    sequence_choice[1] <- mean(pop$males[pop$females[, "ID_father"], "trait"], na.rm = TRUE)
    sequence_observed[1] <- NA
    sequence_preferred[1] <- NA
    # number of individuals in each age class
    pop_size <- matrix(NA, nrow = tmax + 1, ncol = 2*nb)
    pop_size[1,] <- c(n_females, n_males)
  } else {
    mean_choice <- NA
    trad_duration <- NA
  }
  sequence_nb_demos <- sequence_nb_observations <- rep(NA, tmax + 1)
  sequence_nb_demos[1] <- NA
  sequence_nb_observations[1] <- NA
  
  ## start loop on time steps
  ## as females are "initialized" with a mate, reproduction must come first
  for(g in 1 : tmax) {
    
    ## position of females who reproduced in the previous step
    pos_fR <- which(pop$reproductive[pop$females[, "age"]] & !is.na(pop$females[, "ID_father"])) # position of the females who reproduced (demonstrator females)
    # if (length(pos_fR) == 0) { # if no mate choice at the previous step, we want to end the simulation
    #   print(g)
    #   print("no mate choice in the previous step, stop time, return results")
    #   break # go out of the loop on time steps, to store results
    # }
    # print("the next 2 lines are the position of females who mated and of males that have been chosen in the previous step")
    # print(pos_fR)
    # print(pop$females[pos_fR, "ID_father"])
    
    ## observations / preference
    pos_fO <- which(pop$observing[pop$females[, "age"]]) # position of observing females
    if (length(pos_fO) == 0) {
      # print(g)
      # print("no observing females")
    }
    if (parameters$n_obs == 0) { # preference in the null model, for n_obs = 0
      pop$females[pos_fO, "observed"] <- NA # nothing is observed
      pop$females[pos_fO, "preferred"] <- rbinom(length(pos_fO), size = 1, prob = 1/2) # the preference is random
      sequence_nb_observations[g + 1] <- parameters$n_obs
    } else { # preference in the social model, for n_obs > 0
      # we need the trait of the observed males
      if (parameters$female_strategy == "conformity") { # if female strategy is conformity
        traits_observed_males <- pop$males[pop$females[pos_fR, "ID_father"], "trait"] # traits of males who reproduced with the demonstrator females at g-1
      }
      if (parameters$female_strategy == "koinophilia") { # if female strategy is koinophilia
        traits_observed_males <- pop$males[which(pop$reproductive[pop$males[, "age"]]), "trait"] # traits of reproductive males who are present at g
      }
      # if (!return_tradition_only) sequence_nb_demos[g + 1] <- length(traits_observed_males) # store the number of social info available
      
      if (parameters$n_obs >= length(traits_observed_males)) { # case where there is less to observe than what females could: all females observe the same matings
        sequence_nb_observations[g + 1] <- length(traits_observed_males) # the effective number of observations is the number of available matings to observe
        if (length(traits_observed_males) == 0) { # if there is nothing to observe at this time step, the preference is random
          # print(g)
          # print("nothing to observe, because no repro females or repro males at the previous step")
          pop$females[pos_fO, "observed"] <- NA
          pop$females[pos_fO, "preferred"] <- rbinom(length(pos_fO), size = 1, prob = 1/2)
        } else {
          fr_obs <- mean(traits_observed_males) # frequency for trait 1 observed
          pop$females[pos_fO, "observed"] <- fr_obs
          if (fr_obs > 1/2) {
            pop$females[pos_fO, "preferred"] <- 1
          } 
          else {
            if (fr_obs == 1/2) {
              pop$females[pos_fO,"preferred"] <- rbinom(length(pos_fO), size = 1, prob = 1/2)
            } 
            else {
              pop$females[pos_fO,"preferred"] <- 0
            }
          }
        }
      }
      if (parameters$n_obs < length(traits_observed_males)) { # could be replaced by a else?
        # case where there are more matings to observe than what females can observe
        # use a hypergeometric sampling for each female
        # and determine which male trait value is most frequent among the observations
        sequence_nb_observations[g + 1] <- parameters$n_obs
        nb_1 <- sum(traits_observed_males) # number of 1 in all possible observations
        nb_0 <- length(traits_observed_males) - nb_1 # number of 0 in all possible observations
        nb_1_obs <- rhyper(length(pos_fO), nb_1, nb_0, parameters$n_obs) # number of 1 observed by each observer female
        fr_obs <- nb_1_obs/parameters$n_obs # frequency of 1 observed for each observer female
        pop$females[pos_fO, "observed"] <- fr_obs
        pop$females[pos_fO, "preferred"] <- ifelse(fr_obs > 1/2, 1, 0)
        l <- fr_obs == 1/2
        if (any(l)) {
          pop$females[pos_fO[l], "preferred"] <- rbinom(sum(l), size = 1, prob = 1/2)
        } 
      }
    } ## end of acquisition of the preference in the social model
    
    
    ## who will survive?
    pos <- which(pop$females[, "age"] < nb) # position of females of age class in which survival occurs (not the last age class)
    pop$females[pos, "alive"] <- with(parameters, rbinom(length(pos), size = 1, prob = survival)) # survival (1) or death (0) for every female
    pop$females[which(pop$females[, "age"] == nb), "alive"] <- 0 # females of the last age class will all die
    pos <- which(pop$males[, "age"] < nb) # position of males of age class in which survival occurs (not the last age class)
    pop$males[pos, "alive"] <- with(parameters, rbinom(length(pos), size = 1, prob = survival)) # survival (1) or death (0) for every male
    pop$males[which(pop$males[, "age"] == nb), "alive"] <- 0 # males of the last age class will all die
    
    ## how many newborns to keep pop size constant?
    nb_newborns <- parameters$K - sum(pop$males[, "alive"], na.rm = TRUE) - sum(pop$females[, "alive"], na.rm = TRUE)
    if (nb_newborns < 0) stop("Number of newborns cannot be negative!") # check that everything is OK, in principle this could be removed...
    if (nb_newborns > parameters$K) stop("Number of newborns cannot excess max population size!")
    
    ## randomly draw the identity of the mother / transmission of trait
    if (length(pos_fR) > 0) { # if there are mothers...
      if (nb_newborns > 0) { # ... and space for newborns
        if (length(pos_fR) == 1) {
          who_is_the_mother <- pos_fR
        } else {
          who_is_the_mother <- sample(pos_fR, size = nb_newborns, replace = TRUE)
        }
        who_is_the_father <- pop$females[who_is_the_mother, "ID_father"] # position of the father for each newborn
        # compute trait values for newborns, as transmitted by the two parents
        # print("the next 4 lines are who_is_the_mother, who_is_the_father and theirs relative traits")
        # print(who_is_the_mother)
        # print(who_is_the_father)
        # print(pop$females[who_is_the_mother, "trait"])
        # print(pop$males[who_is_the_father, "trait"])
        newborn_trait <- ifelse(runif(nb_newborns) < 0.5, pop$females[who_is_the_mother, "trait"], pop$males[who_is_the_father, "trait"])
        newborn_neutral <- ifelse(runif(nb_newborns) < 0.5, pop$females[who_is_the_mother, "neutral"], pop$males[who_is_the_father, "neutral"])
        # print("next line = generation")
        # print(g)
        # print("next line = newborn_trait")
        # print(newborn_trait)
        # print("next line = newborn_neutral")
        # print(newborn_neutral)
        # print("next 2 lines = mother and father traits")
        # print(pop$females[who_is_the_mother, "trait"])
        # print(pop$males[who_is_the_father, "trait"])
        # mutation
        mutation_rate <- 1 - parameters$heritability
        mute <- runif(nb_newborns) < mutation_rate
        if (any(mute)) newborn_trait[mute] <- rbinom(sum(mute), size = 1, prob = 1/2)
        mute <- runif(nb_newborns) < mutation_rate
        if (any(mute)) newborn_neutral[mute] <- rbinom(sum(mute), size = 1, prob = 1/2)
        # print("the next line is newborns traits, and male newborn traits")
        # print(newborn_trait)
        # are newborns males or females ?
        newborn_female <- runif(nb_newborns) < 0.5 # TRUE if the newborn is a female
      } else {print("no space for newborns at this step")}
    } else {print("no newborns because no repo females at this step")}

    ## ageing / death
    pop$males[pop$males[, "alive"] == 0, ] <- NA # rows of NA for dead males
    # pos_males_age <- which(!is.na(pop$males[, "age"]))
    # pop$males[pos_males_age, "age"] <- ifelse(runif(length(pos_males_age)) > parameters$ageing, pop$males[pos_males_age, "age"] + 1, pop$males[pos_males_age, "age"])
    pop$males[, "age"] <- pop$males[, "age"] + 1
    pop$females[pop$females[, "alive"] == 0, ] <- NA # rows of NA for dead females
    # pos_females_age <- which(!is.na(pop$females[, "age"]))
    # pop$females[pos_females_age, "age"] <- ifelse(runif(length(pos_females_age)) > parameters$ageing, pop$females[pos_females_age, "age"] + 1, pop$females[pos_females_age, "age"])
    pop$females[, "age"] <- pop$females[, "age"] + 1

    
    # add newborns to tables
    if (length(pos_fR) > 0) { # if there are mothers...
      if (nb_newborns > 0) { # ... and space for newborns: there are newborns to add to tables
        if (length(newborn_female) > 0 & any(newborn_female)) { # female newborns
          pos <- which(is.na(pop$females[, "age"])) # all empty spaces in the matrix
          pos_f <- pos[1 : sum(newborn_female)] # the empty spaces that will be used to store newborns
          pop$females[pos_f, "trait"] <- newborn_trait[newborn_female]
          pop$females[pos_f, "neutral"] <- newborn_neutral[newborn_female]
          pop$females[pos_f, "age"] <- 1
          pop$females[pos_f, "alive"] <- 1
        }
        if (length(newborn_female) > 0 & any(!newborn_female)) { # male newborns
          pos <- which(is.na(pop$males[, "age"])) # all empty spaces in the matrix
          pos_m <- pos[1 : sum(!newborn_female)] # the empty spaces that will be used to store newborns
          # print(newborn_trait[!newborn_female])
          pop$males[pos_m, "trait"] <- newborn_trait[!newborn_female]
          pop$males[pos_m, "neutral"] <- newborn_neutral[!newborn_female]
          pop$males[pos_m, "age"] <- 1
          pop$males[pos_m, "alive"] <- 1
        }
      }
    }

    ## update who is reproductive (who will have to choose a mate)
    pos_fR <- which(pop$reproductive[pop$females[, "age"]])
    pos_mR <- which(pop$reproductive[pop$males[, "age"]])
    
    if (length(pos_fR) == 0) {
      print(g)
      print("no more reproductive females in the pop")
      nb_repro <- 0
      break
    }
    if (length(pos_mR) == 0) {
      print(g)
      print("no more reproductive males in the pop, stop time, return results")
      nb_repro <- 0
      break
    }
    
    ## mating based on preferences with competition among males
    pop$females[, "ID_father"] <- NA # reset the column of chosen mates from the previous time step
    pos_f_non_mated <- which(pop$reproductive[pop$females[, "age"]] & is.na(pop$females[, "ID_father"])) # females who have to choose a mate
    pos_m_non_mated <- pos_mR # males who still have matings to make
    nb_matings_males <- rep(parameters$n_matings, length(pos_mR)) # number of matings to make for each male
    names(nb_matings_males) <- pos_mR
    
    while (length(pos_f_non_mated) > 0 && length(pos_m_non_mated) > 0) { # as long as there are males and females available to mate
      # stop the while loop if c = 1 and all females preferences are different from males traits (no choice possible, the while loop would never end)
      if (parameters$c == 1 & !any(pop$females[pos_f_non_mated, "preferred"] %in% pop$males[pos_m_non_mated, "trait"])) {  
        break }
      
      pos_f_met <- pos_m_met <- NULL # vectors of females and males who will meet each other
      if (length(pos_f_non_mated) > length(pos_m_non_mated)) { # if there are more females to mate than males available
        pos_f_met <- sample(pos_f_non_mated, length(pos_m_non_mated), replace = FALSE) # we randomly select the females who will have a mate
        pos_m_met <- pos_m_non_mated
      }
      if (length(pos_m_non_mated) > length(pos_f_non_mated)) { # if there are more males available than females
        pos_m_met <- sample(pos_m_non_mated, length(pos_f_non_mated), replace = FALSE) # we randomly select one male for each female to meet with
        pos_f_met <- pos_f_non_mated
      }
      if (length(pos_m_non_mated) == length(pos_f_non_mated) && length(pos_m_non_mated) > 1) { # if there are as many males as females available (and > 1)
        pos_m_met <- sample(pos_m_non_mated) # all males available are met, we just reorder the vector of males
        pos_f_met <- pos_f_non_mated
      }
      if (length(pos_m_non_mated) == length(pos_f_non_mated) && length(pos_m_non_mated) == 1) { # if there are as many males as females available (1 of each)
        pos_m_met <- pos_m_non_mated # the last male available is met by the last female 
        pos_f_met <- pos_f_non_mated
      }
      
      pr <- ifelse(pop$males[pos_m_met, "trait"] == pop$females[pos_f_met, "preferred"], # probability to accept the sampled (met) male, for each female
                   parameters$c,
                   1 - parameters$c)
      
      if(any(is.na(pop$males[pos_m_met, "trait"]))) stop("NA in met males traits")
      if (any(pr < 0) | any(pr > 1)) stop("Invalid probability value!")
      
      chosen <- rbinom(length(pos_f_met), size = 1, prob = pr) # is the male chosen or not, for every female
      
      if (any(chosen == 1)) { # if matings occur, store the choice of females and update the vectors of males and females still available for mating
        pop$females[pos_f_met[chosen == 1], "ID_father"] <- pos_m_met[chosen == 1]
        pos_f_non_mated <- pos_fR[is.na(pop$females[pos_fR, "ID_father"])] # update positions of females without a mate
        
        nb_matings_males[as.character(pos_m_met[chosen == 1])] <- nb_matings_males[as.character(pos_m_met[chosen == 1])] - 1 # remove a mating for males who mated
        pos_m_non_mated <- pos_m_non_mated[nb_matings_males > 0] # keep only males who still have matings to make
        nb_matings_males <- nb_matings_males[nb_matings_males > 0] # keep only males who still have matings to make
      }
    }
    
    ## number of demonstrations during that time step
    nb_repro <- length(which(pop$reproductive[pop$females[, "age"]] & !is.na(pop$females[, "ID_father"])))
    sequence_nb_demos[g + 1] <- nb_repro
    
    ## done with that time step, store output
    
    # if (!any(!is.na(pop$females[, "age"])) || !any(!is.na(pop$males[, "age"]))) {
    #   g <- g - 1
    #   print("male or female extinction")
    #   break # go out of the for (generations) loop
    # }
    
    if (!return_tradition_only) { # outputs if we store raw temporal results
      # for males
      sequence_trait_m[g + 1] <- mean(pop$males[, "trait"], na.rm = TRUE)
      sequence_neutral_m[g + 1] <- mean(pop$males[, "neutral"], na.rm = TRUE)
      # for females
      sequence_trait_f[g + 1] <- mean(pop$females[, "trait"], na.rm = TRUE)
      sequence_neutral_f[g + 1] <- mean(pop$females[, "neutral"], na.rm = TRUE)
      sequence_observed[g + 1] <- mean(pop$females[, "observed"], na.rm = TRUE)
      sequence_preferred[g + 1] <- mean(pop$females[, "preferred"], na.rm = TRUE)
      sequence_choice[g + 1] <- mean(pop$males[pop$females[, "ID_father"], "trait"], na.rm = TRUE)
      # number of individuals in each age class
      n_females <- table(factor(pop$females[, "age"], levels = 1 : nb))
      n_males  <- table(factor(pop$males[, "age"], levels = 1 : nb))      
      pop_size[g + 1, ] <- c(n_females, n_males)
    }
    else { # if we store results for one tradition only, compute the duration of the tradition and stops the simulation if the tradition ends
      if (g == parameters$n_steps_removed) {
        mean_choice <- mean(pop$males[pop$females[, "ID_father"], "trait"], na.rm = TRUE)
        trad_duration <- 0
      }
      if (g > parameters$n_steps_removed) { # start computing the duration of the tradition and return results if the tradition ends
        mean_choice_old <- mean_choice # mean choice at g-1
        mean_choice <- mean(pop$males[pop$females[, "ID_father"], "trait"], na.rm = TRUE) # mean choice at g
        # print("the next 2 lines are mean_choice_old and mean_choice")
        # print(mean_choice_old)
        # print(mean_choice)
        # XOR operator (exclusive or): TRUE if one input is FALSE and the other is TRUE; FALSE if the two are TRUE, or if the two are FALSE
        if (xor(mean_choice > 0.5, mean_choice_old > 0.5)) return(c(last_time_step = g, mean_choice_last_step = mean_choice_old, 
                                                                    n_repro_last_step = nb_repro, mean_demos = mean(sequence_nb_demos, na.rm = T),
                                                                    mean_n_obs = mean(sequence_nb_observations, na.rm = T),
                                                                    tradition_duration = trad_duration))
        trad_duration <- trad_duration + 1
      }
    }
  } ## end of for loop (generations)
  
  ## done with all time steps, return results
  if (!return_tradition_only) {
    results <- list(population = pop, 
                    male_trait = sequence_trait_m, male_neutral = sequence_neutral_m,
                    female_trait = sequence_trait_f, female_neutral = sequence_neutral_f,
                    observed = sequence_observed, preference = sequence_preferred,
                    choice = sequence_choice, effective_nobs = sequence_nb_observations,
                    n_demos = sequence_nb_demos,
                    size = pop_size)
    return(results)
  } 
  else {
    return(c(last_time_step = g, mean_choice_last_step = mean_choice, 
             n_repro_last_step = nb_repro, mean_demos = mean(sequence_nb_demos, na.rm = T),
             mean_n_obs = mean(sequence_nb_observations, na.rm = T),
             tradition_duration = trad_duration))
  }
}

# p$c <- 0.7
# p$n_matings <- 2
# (test <- dynamics(150, p, return_tradition_only = T))

# test code
p$K <- 200
p$nb_class <- c(1,1)
p$ageing <- 0.1
p$female_strategy <- "conformity"
p$n_obs <- 0
p$c <- 0.5
p$heritability <- 1
p$n_matings <- 10
p$n_steps_removed <- 100
sim <- dynamics(500, p, return_tradition_only = F)

# test plot
plot(sim$male_trait, type = "l", col = "blue", ylim = c(0,1))
lines(sim$male_neutral, col = "gray")
lines(sim$choice, col = "red", lty = 3)
abline(h = 0.5, lty = 2)


# nb <- sum(p$nb_class)
# coul <- NULL
# for(i in 1 : nb) coul <- c(coul, adjustcolor("indianred", alpha.f = i/nb))
# for(i in 1 : nb) coul <- c(coul, adjustcolor("steelblue", alpha.f = i/nb))
# barplot(t(sim$size[1:100, ]), col = coul)

## function inputs:
# param_grid = dataframe for parameters variation
# parameters = list of global parameters
# file_name = name of the file in which we store results
# n_replicates = number of replicates for each parameters combination
# tmax = number of time steps in each simulation

param_grid <- param_variation
parameters <- p
file_name <- "test.csv"
n_replicates = 10
tmax = 100

parameters_exploration <- function(param_grid, parameters, file_name, n_replicates = 10, tmax = 200) {
  n_simus <- nrow(param_grid) * n_replicates
  simu_index <- 1
  for (i in 1 : nrow(param_grid)) {
    parameters[colnames(param_grid)] <- param_grid[i,] # set parameters value for the simulation to run
    # print(param_grid[i,])
    for (r in 1 : n_replicates) {
      cat(simu_index, "/", n_simus, "\n")
      simu <- try(dynamics(tmax, parameters, return_tradition_only = TRUE)) # run the simulation
      if (simu_index == 1) { # set up the file to store results
        cat("", "simu", colnames(param_grid), "replicate", names(simu), sep = ";", file = file_name, append = FALSE) # columns names in the file
        cat("\n", file = file_name, append = TRUE) # skip a line
      }
      if ("try-error" %in% class(simu)) { # try ro run the simulation, return NA in the results if there is an error
        # c(generation = NA, mean_choice = NA, tradition_duration = NA) # to be changed in something more generic
        cat(simu_index, i, unlist(param_grid[i,]), r, rep(NA, length(simu)), sep = ";", file = file_name, append = TRUE)
        cat("\n", file = file_name, append = TRUE) # skip a line
      }
      else {
        cat(simu_index, i, unlist(param_grid[i,]), r, simu, sep = ";", file = file_name, append = TRUE) # store results of the simulation in one row of the file
        cat("\n", file = file_name, append = TRUE) # skip a line
      }
      simu_index <- simu_index + 1
    }
  }
}

# test code
# gr <- expand.grid(heritability = seq(0.1, 0.9, by = 0.1), c = seq(0.5, 1, length = 10))
# parameters_exploration(gr, p, "toto.csv", tmax = 100)
parameters_exploration(param_grid, parameters, file_name, n_replicates = 10, tmax = 200)


