# Interplay_Top-Down_Bottom-up_Biodiversity

Analysis code for paper: "Complex Interplay Between Top-down and Bottom-up Drivers of Biodiversity Maintenance"

## Overview

This repository contains R scripts for analyzing the effects of apex-predator presence/absence and nutrient addition on freshwater eukaryotic biodiversity in experimental lake ecosystems.

**Experimental Design:**
- 12 artificial lakes (CNRS-ENS platform, Nemours, France)
- 2×2 factorial design: Predator (present/absent) × Nutrients (enriched/control)
- Monthly sampling from April 2022 to April 2023
- Environmental DNA (eDNA) metabarcoding approach

## Repository Structure

```
.
├── R_codes/
│   ├── 12_24_ChapI_Treatment_ALPHAdiversity_DNAcDNA.R
│   └── 0125-beta_div_FAMILY_Clean.R
└── input_files/
    ├── 1_Euka_DNA_Taxo_verified_pres_abs.csv
    ├── sample_mapping_names.csv
    ├── sample_names.csv
    └── Samples_info_Planaqua.csv
```

## Analyses Performed

### 1. Alpha Diversity Analysis
**Script:** `12_24_ChapI_Treatment_ALPHAdiversity_DNAcDNA.R`

Calculates taxonomic richness at multiple hierarchical levels:
- MOTU (Molecular Operational Taxonomic Unit) richness
- Phylum richness
- Class richness
- Order richness
- Family richness

**Statistical approach:**
- Generalized Linear Mixed Models (GLMMs) using `glmmTMB`
- Fixed effects: Predator presence × Nutrient addition
- Random effect: Month of collection (to account for temporal autocorrelation)
- Model validation with `DHARMa` package

**Output:**
- Treatment effect estimates
- Interaction effects
- Model diagnostics
- Formatted tables and figures

### 2. Beta Diversity Analysis
**Script:** `0125-beta_div_FAMILY_Clean.R`

Examines community composition changes in response to experimental treatments.

**Analyses include:**
- Canonical Correspondence Analysis (CCA) using `vegan`
- Family-level community composition
- Treatment effects on community structure
- Removal of sporadic families (presence threshold analysis)
- Biplot visualization of CCA results

**Key features:**
- Focuses on aquatic families only
- Monthly data (excludes daily monitoring samples)
- Condition on collection month to account for temporal effects
- Permutation tests (999 permutations) for statistical significance

**Visualizations:**
- CCA biplots showing species-environment relationships
- Family presence distribution plots
- Treatment effect visualizations

## Input Files

### Required Data Files
All input files should be placed in the `input_files/` directory:

1. **1_Euka_DNA_Taxo_verified_pres_abs.csv**
   - Taxonomically verified presence/absence matrix
   - MOTUs × samples
   - Includes taxonomic hierarchy (phylum, class, order, family)
   - Aquatic/terrestrial classification

2. **sample_mapping_names.csv**
   - Maps sample IDs to simplified names
   - Links sequencing IDs to biological samples

3. **sample_names.csv**
   - Sample metadata from sequencing
   - PCR plate information
   - Sampling dates

4. **Samples_info_Planaqua.csv**
   - Detailed sample metadata
   - Treatment information
   - Collection dates and times
   - Lake IDs
   - Sampling type (Monthly/Daily)

## Requirements

### R Version
R ≥ 4.3.3

### Required Packages
```r
# Data manipulation
library(tidyverse)

# Statistical modeling
library(glmmTMB)
library(lme4)
library(DHARMa)
library(car)
library(effectsize)
library(emmeans)

# Community ecology
library(vegan)

# Visualization
library(cowplot)
library(patchwork)
library(ggrepel)

# Tables and reporting
library(broom.mixed)
library(flextable)
library(officer)
library(kableExtra)
library(knitr)
```

## Usage

### 1. Set Working Directory
Update the working directory path at the beginning of each script:
```r
setwd("path/to/your/Chapter_1_directory")
```

### 2. Run Alpha Diversity Analysis
```r
source("R_codes/12_24_ChapI_Treatment_ALPHAdiversity_DNAcDNA.R")
```

### 3. Run Beta Diversity Analysis
```r
source("R_codes/0125-beta_div_FAMILY_Clean.R")
```

## Data Processing Notes

- **Aquatic taxa only**: Analysis focuses on aquatic families to match the study environment
- **Monthly samples**: Daily monitoring samples excluded as they belong to a separate experimental study
- **Family threshold**: Rare families present in <5% of samples removed from beta diversity analysis
- **Percidae excluded**: Apex predator family (European perch) removed from community analyses

## Related Repositories

- **Data Processing Pipeline:** [Lake_biodiversity_sequencingData_processing](https://github.com/JulieGuenat/Lake_biodiversity_sequencingData_processing/)
  - OBITools 4 processing scripts
  - metabaR cleaning workflows
  - Taxonomic verification protocols

## Citation

If you use this code or data, please cite:

[paper citation - to be added upon publication]

**Data Availability:**
- Raw sequencing data: NCBI BioProject [PRJNA######]
- Processed data: Available in this repository

## Acknowledgments

This research was conducted at the CNRS-ENS artificial lake platform (CEREEP-Ecotron) in Nemours, France, in collaboration with the Laboratory of Conservation Biology (LBC), University of Lausanne, Switzerland.

**Code Development:**
Part of this analysis code was developed with assistance from Claude (Anthropic AI), particularly for data formatting, visualization optimization, and statistical modeling workflows.


## Contact

Julie Guenat   
University of Lausanne, Switzerland  
julie.guenat@ik.me

## Version History

- **v1.0** (January 2025) - Initial release for manuscript submission

---

**Last updated:** January 2025
