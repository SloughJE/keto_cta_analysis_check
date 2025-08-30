# KETO-CTA Study – Audit & Reproduction

This repository contains an audit of the published KETO-CTA paper [Longitudinal Data From the KETO-CTA Study: Plaque Predicts Plaque, ApoB Does Not; **Soto-Mota et al. (2025)**](https://www.jacc.org/doi/10.1016/j.jacadv.2025.101686) using the released dataset and reproducible R/Quarto code. 

## Data source and variable scope

CT angiography–derived plaque-burden metrics (PAV, TPS, CAC, NCPV) at baseline and follow-up were obtained from the 
Citizen Science Foundation [keto-CTA repository](https://citizensciencefoundation.org/keto-cta/). These are the same per-participant plaque metrics used in the publication; 
basic levels and change checks match the reported values. No other variables (e.g., ApoB, LDL or other lipids, demographics) are 
available. Analyses are therefore restricted to specifications that rely solely on plaque-burden metrics; 
models requiring non-plaque variables could not be examined.

## Very brief summary of the analysis
- Recreates key figures from the public dataset and compares them to the paper; several panels/captions do not match the data.
- Objective diagnostics (RESET, Breusch–Pagan, Shapiro–Wilk) show violations of linear-model assumptions used in the paper; robust regression does not fix nonlinearity/heteroskedasticity.
- Some stated “no association” conclusions (e.g., TPS vs CAC) are not reproduced; reliance on univariable change-score models is fragile.
- Baseline-adjusted models (ANCOVA) are shown as the appropriate framework for plaque progression.

## Directory structure
```text
repo-root/
├── index.qmd                 # Slide deck (Quarto, Reveal.js)
├── index.html                # Rendered slides (self-contained)
├── styles.css                # Slide styles
├── code/                     # Reusable R scripts
│   ├── charts.R
│   ├── charts_extended.R
│   ├── descriptive_statistics.R
│   ├── linear_models.R
│   ├── modeling_work.R
│   └── tables.R
├── data/                     # Study dataset(s)
│   └── ...
├── figures/                  # Derived figures and diagnostics
│   ├── diagnostics/
│   ├── plots_plaque_metrics/
│   └── core panels (e.g., Figure1A/B, Figure2F)
├── assets/                   # Reference images extracted from paper
│   └── ...
├── papers/                   # Preprint and related materials
│   └── ...
└── keto_cta_data_analysis.Rproj
```

