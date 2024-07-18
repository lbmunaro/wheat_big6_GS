# Multi-trait GBLUP ----

# Objective ----
# - Run Multi-trait GBLUP using the combination of trait and trial as groups

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
                        na.action = na.method(y='include',x='include'),
                        workspace = '64gb', maxit= 20)

mod0_MT.GBLUP <- update.asreml(mod0_MT.GBLUP)
mod0_MT.GBLUP <- update.asreml(mod0_MT.GBLUP)
mod0_MT.GBLUP <- update.asreml(mod0_MT.GBLUP)

# Save ----
save.image('Data/Step5_MT-GBLUP_fa3.RData')


