# Air Pollution and Tuberculosis: A Systematic Review, Meta-Analysis, and Burden Estimation Study

This repository contains the code used in the manuscript **"Air Pollution and Tuberculosis: A Systematic Review, Meta-Analysis, and Burden Estimation Study"**.

## Contents

- **meta_analysis.R** — Performs all random-effects meta-analyses and generates corresponding forest plots.
- **burden_calculations_hap.R**, **burden_calculations_oap.R**, **burden_calculations_pm2.5.R** — Estimate the burden of TB deaths attributable to PM2.5 exposure from:
  - Household air pollution (HAP)
  - Outdoor air pollution (OAP)
  - Overall PM2.5 exposure  
  These scripts also generate uncertainty intervals and tables.
- **descriptive_tables.R** — Generates Tables E1 and E2.
- **other_analyses.R** — Produces the maps presented in the manuscript.

## Related Resources

- **GBD modeling code** (MR-BRT, ST-GPR, DisMod-MR) is available via the [Institute for Health Metrics and Evaluation (IHME) GitHub page](https://github.com/ihmeuw).
- **Data and supplementary materials** are available on the [Open Science Framework (OSF) project page](https://osf.io/s3tyb/overview?view_only=aaf4f00051ad435783486ab7c59ee51e).

## Corresponding Author Information

Siobhan Carroll, PhD Candidate in Epidemiology at McGill University

siobhan.carroll@mail.mcgill.ca

## Requirements

These analyses were developed using **R version 4.5.1**.

**Required packages**
```r
# Install required packages
install.packages(c("tidyverse", "here", "meta", "readODS", "readxl", "kableExtra", "webshot2", "chromote",
"rnaturalearth", "tmap", "xtable", "PrettyCols", "data.table", "officer", "flextable"))
```

## Usage
**Reproducibility note:**
All file paths in the code are relative to the project root directory.
To reproduce the analyses, clone this repository and open the .Rproj file in RStudio.
Required data files (available via the OSF link above) should be saved in the same directory as the .Rproj file.
Scripts can be run in any order.
