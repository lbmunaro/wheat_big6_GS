# Multi-trait GBLUP ----

# Objective ----
# - Predict GEBVs

rm(list=ls()) # clean workspace

# Packages ----

library(tidyverse) # R packages for data science
library(asreml) # ASReml-R package

# Load data ----
## BLUES & Gmat ----
load('Data/Step5_MT-GBLUP_fa3.RData')

# Predict GEBVs ----
GEBVs_MT <- predict(mod0_MT.GBLUP, classify='germplasm:group', pworkspace='16gb')$pvals

# Save ----
save.image('Data/Step6_GEBVs_MT-GBLUP_fa3.RData')
