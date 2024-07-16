# Multi-trait GBLUP ----

# Objective ----
# - Run Multi-trait GBLUP using the combination of trait and trial as groups
# - Get GEBVs

rm(list=ls()) # clean workspace

# Packages ----

library(tidyverse) # R packages for data science
library(asreml) # ASReml-R package

# Load data ----
## BLUES & Gmat ----
load('Data/Pheno_and_Gmat.RData')

# Check data structure and variable names
blues_group <- blues_filter |>
  mutate(group = as.factor(paste(trait,study,sep = '_'))) |>
  glimpse()

# MT-GBLUP ----

mod0_MT.GBLUP <- asreml(fixed = predicted.value ~ group,
                        random = ~fa(group,3):vm(germplasm, Ginv.sparse),
                        weights = weight,
                        data = blues_group,
                        family = asr_gaussian(dispersion = 1),
                        na.action = na.method(y='omit',x='omit'),
                        workspace = '6gb')

save.image('Data/Step5_MT-GBLUP.RData')

GEBVs_MT <- predict(mod0_MT.GBLUP, classify='germplasm:group', pworkspace='8gb')$pvals

save.image('Data/Step5_MT-GBLUP.RData')

load('Data/Step5_MT-GBLUP.RData')

summary(mod0_MT.GBLUP, coef = T)$coef.random


