# Bayesian Inference via Gibbs Sampling

University coursework (Bayesian Statistics) implementing a Gibbs sampler from scratch in R to estimate the posterior distribution of the mean and variance of a Normal model, then comparing it against a Metropolis-Hastings (accept/reject) sampler on the same data.

## What it does

- **Simulated data**: draws 50 samples from a known N(5, 2) distribution as the working dataset.
- **Gibbs sampler**: implements the full conditional updates for θ (mean) and φ (precision) by hand — 1000 iterations from a Normal-Gamma prior — and plots the trace and posterior histograms for both parameters (discarding the first half as burn-in).
- **Joint posterior**: visualizes the joint distribution of θ and φ as an interactive 3D surface (Plotly).
- **Metropolis-Hastings comparison**: re-implements the same estimation with an accept/reject step on each parameter update, to compare convergence and posterior shape against the Gibbs sampler.
- **Report** (`PDF_BayesianStatistics.pdf`): write-up with the full derivation and interpretation of results.

## Tech stack

R, ggplot2, plotly, dplyr, pastecs, plot3D.

## Project structure

```
Jupyter_BayesianStatistics.r   # Gibbs sampler + Metropolis-Hastings, from scratch
PDF_BayesianStatistics.pdf     # Report with derivations and interpretation
```

## How to run

```r
install.packages(c("ggplot2", "dplyr", "plotly", "gridExtra", "pastecs", "plot3D"))
```

Update the hardcoded `setwd()` path at the top of `Jupyter_BayesianStatistics.r`, then run it top to bottom in R/RStudio.
