# SAlient Multivariate BetA Mixture Model (SAMBA-Mix)

Implementation for the paper:

**Salient Bayesian Clustering for Proportional Data via a Multivariate Beta Mixture Model**  
Rodríguez & Stephenson (2026)

---

## Overview

This repository implements the **SAMBA-Mix** a Bayesian framework for clustering multivariate proportional data.

Proportional data (values in (0,1)) are common in applications such as social determinants of health (SDoH), genomics, and ecology, but present challenges due to:
- bounded support,
- skewness,
- and heterogeneous feature relevance.

SAMBA-Mix addresses these by combining:
- a **multivariate beta mixture model** (no transformations required),
- an **overfitted Bayesian mixture** to infer the number of clusters,
- and **feature saliency indicators** to identify variables that drive clustering.

---

## Key Features

- Bayesian clustering for **bounded continuous data**
- Learns the **number of clusters automatically**
- Performs **joint feature selection + clustering**
- Robust to **irrelevant/noisy variables**
- Fully implemented in **R using NIMBLE (MCMC)**

---


## Repository Structure
SAMBA-Mix/

- Application/ # Real data analyses (AHRQ SDoH)
  - figures_tabs/ # Figures and tables from paper
  
  
-  Simulations/ # Simulation studies from paper

- R/ # Core model implementation
  - bayesmbmm_OFMM_SF_v2.R # Model definition
  - sampler_Utils.R # MCMC sampling utilities
  - helpers.R # Supporting functions

- FASRC/ # Scripts for cluster computing (HPC)

---

## Installation

### Requirements

- R (≥ 4.3)
- Recommended packages:
```r
install.packages(c(
  "nimble",
  "label.switching",
  "coda",
  "ggplot2",
  "dplyr"
))
```

## Quick Start Example

Below is a minimal example showing how to load the model and fit SAMBA-Mix.

### Load dependencies and source code

```r
library(nimble)

# Update path as needed
source("~/bayesmbmm/R/helpers.R")
source("~/bayesmbmm/R/bayesmbmm_OFMM_SF_v2.R")

fit <- Ombmm_FeatSaliency(
  x = Xdat,        # n x p matrix with values in (0,1)
  K = 30,           # Upper bound on number of clusters
  niter = 20000,    # Total MCMC iterations
  nburnin = 10000,  # Burn-in
  nchain = 1,

  # Cluster concentration (controls separation)
  g = 15, h = 1,    # Higher values → tighter clusters

  # Prior for non-informative features
  a_mu0 = 1, b_mu0 = 1, 
  # Flat prior prevents bias toward assigning features as irrelevant

  # Feature saliency prior
  a_rho = 2, b_rho = 5, 
  # Encourages moderate feature inclusion (~28% prior expectation)

  # Mixture sparsity (controls number of active clusters)
  e0 = 0.01,        
  # Smaller values → fewer active components

  priorOnAlpha = "no",
  plot_trace = FALSE
)
```
---

## Notes

**Input data (`Xmat`) must:**
- Be numeric
- Have all values strictly in **(0,1)** (no 0 or 1)

**Tuning tips:**
- Increase `g` → tighter clusters (less overlap)
- Decrease `e0` → fewer clusters retained
- Adjust (`a_rho`, `b_rho`) → controls number of selected features

**Additional resources:**
- `Application/` — real data examples  
- `Simulations/` — reproducible experiments

## Citation (to add)

Contact information:  **crodriguezcabrera@g.harvard.edu**
