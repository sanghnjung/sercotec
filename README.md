# Replication Package: Tailoring Government Support for SME's Investment and Financial Health: Balancing Advisory Services and Training Programs.

## Project Overview
This repository contains the data and code required to replicate the analysis of the causal impact of **Advisory vs. Training services** on SME financial outcomes. The study addresses potential selection bias by utilizing a high-dimensional covariate set (100+ variables) and implementing two primary econometric strategies:

1.  **Ridge-Regularized Propensity Score Weighting**: Including Average Treatment Effects (ATE), Treatment on the Treated (ATT), Treatment on the Control (ATC), and Overlap Weighting (ATO).
2.  **Post-Double Selection (PDS) Lasso**: Utilizing weighted selection for binary outcomes to satisfy orthogonality conditions and ensure valid post-selection inference.

## Data Availability
The dataset provided (`Synthetic_Replication_Data.csv`) is a **synthetic version** of the original administrative records.

* **Methodology**: Generated via a **Gaussian Copula Synthesizer** to maintain the joint distribution, correlation structure, and marginal properties of the original Chilean SME dataset.

### R Environment
Required packages can be installed via:
```r
install.packages(c("dplyr", "glmnet", "PSweight", "doParallel", "fastDummies", "sandwich", "lmtest"))
