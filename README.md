# Gene-Culture Coevolution and Traditions in Mate Choice

## Overview

This repository contains the simulation models, simulations data, and scripts for data analysis and figure generation for the research article:

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
  
## Key Parameters

The models are controlled by three key transmission fidelity parameters:
- **n_obs**: Number of mating observations per learning female
- **c**: Copying fidelity, i.e. probability of following learned preference at each male encounter 
- **mutation_rate**: Genetic mutation rate affecting male trait heritability 

Additional demographic parameters:
- **K**: Population size (default: 1000)
- **survival**: Survival probability (default: 0.9)
- **ageing**: Aging rate (default: 0.9)
- **n_matings**: Maximum matings per male (default: 100)

## Content and replicability

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

- model/ contains the script for our two Individual-Based Models. 
  
- main_results/ contains, for each figure of the main text:
    - simulation data used to generate the figure
    - the script analyzing data and generating the figure
    - the corresponding figure in .pdf format
      
- supplementary_materials/ contains, for each figure of the appendix:
    - simulation data used to generate the figure
    - the script analyzing data and generating the figure
    - the corresponding figure in .pdf format

Simulations data in main_results/ and supplementary_materials/ can be replicated using `core_model.R` with appropriate parameter values. These parameter values can be found in the script generating each figure.

## Installation

### Requirements
- **R version**: 4.4.0 or higher

#### Warning
⚠️ **Compute time**: 24-48 hours for full parameter exploration (Figure 4 data generation)
Computations were performed on the [Migale cluster](https://migale.inrae.fr/cluster). 

## Contact
For questions about the models or data, contact the corresponding author in the paper.
