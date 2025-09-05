# Gene-Culture Coevolution and Traditions in Mate Choice

## Overview

This repository contains the simulation models, analysis scripts, and figure generation code for the research article:

**"Gene-Culture Coevolution Favours the Emergence of Traditions in Mate Choice through Conformist Social Learning"**

*Authors: D. Federico, F.X. Dechaume-Moncharmont, J.B. Ferdy, A. Pocheville*

## Research Question

This theoretical study investigates the conditions under which traditions emerge in mate choice when female preferences are socially transmitted through conformist learning while male traits are genetically inherited. We challenge the assumption that sophisticated social learning mechanisms are required for tradition formation by demonstrating how gene-culture coevolutionary dynamics can fundamentally alter these conditions.

## Models

We compare two individual-based models to understand how gene-culture interactions facilitate tradition formation:

### Cultural Model
- **Pure social learning** without resource constraints
- Females observe mating demonstrations and can always find their preferred male type
- Represents classical cultural evolution scenarios where social transmission is the primary evolutionary force
- Traditions emerge only when conformity exhibits majority exaggeration

### Gene-Culture Model  
- **Social learning with encounter constraints** based on male type availability
- Female choice depends on both preference and demographic composition
- Creates gene-culture coevolutionary dynamics through positive frequency-dependent interactions
- Traditions can emerge even without majority exaggeration in conformist learning

## Repository Structure

```
traditions-model/
├── model/
│   └── core_model.R             # Main simulation models (Cultural & Gene-Culture)
├── main_results/
│   ├── figure_2/                 
│   ├── figure_3/                 
│   └── figure_4/                 
├── supplementary_materials/
│   ├── A/                    
│   ├── B/                 
│   ├── C/                  
│   ├── D/                  
│   └── E/                  
└── README.md
```

## Key Parameters

The models are controlled by three key transmission fidelity parameters:

- **n_obs**: *n_obs$ Number of mating observations per learning female (1 to 30, plus unlimited)
- **c**: Copying fidelity - probability of following learned preference (0.5 to 1.0)
- **mutation_rate**: Genetic mutation rate affecting male trait heritability (0.05 to 1.0)

Additional demographic parameters:
- **K**: Population size (default: 1000)
- **survival**: Survival probability (default: 0.9)
- **ageing**: Aging rate (default: 0.9)
- **n_matings**: Maximum matings per male (default: 100)

## Installation

### Requirements

- **R version**: 4.4.0 or higher
- **Required packages**:
  ```r
  install.packages(c(
    "ggplot2",
    "viridis", 
    "patchwork",
    "dplyr",
    "tidyr"
  ))
  ```

### System Requirements

- **Memory**: Minimum 8 GB RAM for full parameter exploration
- **Storage**: ~2 GB for complete simulation datasets
- **Compute time**: 24-48 hours for full parameter exploration (Figure 4 data generation)

## Usage

### Quick Start

```r
# Load the main model
source("model/core_model.R")

# Set basic parameters
parameters <- list(
  nb_class = c(1,1),
  K = 1000,
  female_strategy = "conformity",
  survival = 0.9,
  ageing = 0.9,
  initial_trait_frequency = 0.5,
  n_matings = 100,
  mutation_rate = 0.1,
  n_obs = 10,
  c = 0.7,
  n_steps_removed = 100
)

# Run single simulation (Cultural Model)
result_cultural <- dynamics_exploration(
  parameters, 
  tmax = 1000, 
  constrained_mating = FALSE
)

# Run single simulation (Gene-Culture Model)
result_gene_culture <- dynamics_exploration(
  parameters, 
  tmax = 1000, 
  constrained_mating = TRUE
)
```

### Reproducing Figures

Each figure can be reproduced independently using the provided data and scripts:

```r
# Example: Reproduce Figure 2 (Conformity Curves)
source("main_results/scripts/script_figure2.R")

# Example: Reproduce Supplementary Figure A
source("supplementary_materials/scripts/script_figure5.R")
```

### Full Parameter Exploration

⚠️ **Warning**: Complete parameter exploration requires significant computational resources.

```r
# Load simulation framework
source("model/core_model.R")

# Define parameter grid (as used in the study)
param_variation <- expand.grid(
  n_obs = c(seq(1, 30, 2), 224),
  mutation_rate = c(1, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05),
  c = seq(0.5, 1, 0.05)
)

# Run exploration (this will take 24-48 hours)
parameters_exploration(
  param_variation, 
  parameters, 
  "results.csv", 
  n_replicates = 100, 
  tmax = 1000
)
```

## Main Results

Our key findings demonstrate:

1. **Cultural Model**: Traditions require majority exaggeration in conformist learning
2. **Gene-Culture Model**: Traditions can emerge without majority exaggeration due to positive frequency-dependent interactions
3. **Parameter Space**: Even minimal conformist capacities (n_obs=3, c=0.6) can generate persistent traditions when traits and choices coevolve
4. **Single Copying**: Can produce long-lasting traditions under gene-culture coevolution, contradicting previous theoretical predictions

## Data Availability

All simulation data used to generate figures are provided in the repository:
- Raw simulation outputs in `/data/` folders
- Pre-processed analysis datasets
- Parameter specifications for each figure

## Computational Notes

**Original Simulations**: Performed on the Migale computing cluster
**Reproduction**: Can be run on standard desktop computers for individual figures
**Performance**: Use optimized functions for large-scale parameter explorations

## Citation

If you use this code or data, please cite:

```
Federico, D., Dechaume-Moncharmont, F.X., Ferdy, J.B., & Pocheville, A. (2025).
Gene-Culture Coevolution Favours the Emergence of Traditions in Mate Choice 
through Conformist Social Learning. [Journal], [Volume], [Pages].
```

## License

This work is licensed under [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

## Contact

For questions about the models or data:
- **Déborah Federico**: [email]
- **Repository issues**: [GitHub issues link]

## Acknowledgments

Simulations were performed using the Migale bioinformatics facility.
