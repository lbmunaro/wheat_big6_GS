# Single-trait GBLUP ----

# Objective ----
# - Run Single-trait GBLUP
# - Get GEBVs

rm(list=objects()) # clean workspace

# Packages ----
library(tidyverse) # R packages for data science
library(asreml) # ASReml-R package

# Load data ----
## BLUES & Gmat ----
load('Data/Pheno_and_Gmat.RData')

# Check data structure and variable names
blues_filter |>
  glimpse()


# ST-GBLUP ----

### Grain yield ----
# model
mod2_YLD <- asreml(fixed=predicted.value~1,
                  random=~study+vm(germplasm, Ginv.sparse):fa(study,2),
                  weights=weight,
                  family=asr_gaussian(dispersion = 1),
                  data=blues_filter |> filter(trait == 'grain_yield_bu_ac'),
                  na.action = na.method(y='include', x='include'),
                  workspace='4gb', maxit= 20)
mod2_YLD <- update.asreml(mod2_YLD)
mod2_YLD <- update.asreml(mod2_YLD)
mod2_YLD <- update.asreml(mod2_YLD)

# GEBVs
GEBVs_mod2_YLD <- summary(mod2_YLD, coef = T)$coef.random |>
  as.data.frame() |>
  rownames_to_column() |>
  glimpse()

hist(GEBVs_mod2_YLD$solution)

### Test weight----
# model
mod2_TW <- asreml(fixed=predicted.value~1,
                  random=~study+vm(germplasm, Ginv.sparse):fa(study,2),
                  weights=weight,
                  family=asr_gaussian(dispersion = 1),
                  data=blues_filter |> filter(trait == 'test_weight_lb_bu'),
                  na.action = na.method(y='include', x='include'),
                  workspace='4gb', maxit= 20)
mod2_TW <- update.asreml(mod2_TW)
mod2_TW <- update.asreml(mod2_TW)
mod2_TW <- update.asreml(mod2_TW)

# GEBVs
GEBVs_mod2_TW <- summary(mod2_TW, coef = T)$coef.random |>
  as.data.frame() |>
  rownames_to_column() |>
  glimpse()

hist(GEBVs_mod2_TW$solution)

### Heading time----
# model
mod2_HD <- asreml(fixed=predicted.value~1,
                  random=~study+vm(germplasm, Ginv.sparse):fa(study,2),
                  weights=weight,
                  family=asr_gaussian(dispersion = 1),
                  data=blues_filter |> filter(trait == 'heading_time_jd'),
                  na.action = na.method(y='include', x='include'),
                  workspace='4gb', maxit= 20)
mod2_HD <- update.asreml(mod2_HD)
mod2_HD <- update.asreml(mod2_HD)
mod2_HD <- update.asreml(mod2_HD)

# GEBVs
GEBVs_mod2_HD <- summary(mod2_HD, coef = T)$coef.random |>
  as.data.frame() |>
  rownames_to_column() |>
  glimpse()

hist(GEBVs_mod2_HD$solution)

### Plant height----
# model
mod2_HT <- asreml(fixed=predicted.value~1,
                  random=~study+vm(germplasm, Ginv.sparse):fa(study,2),
                  weights=weight,
                  family=asr_gaussian(dispersion = 1),
                  data=blues_filter |> filter(trait == 'plant_height_in'),
                  na.action = na.method(y='include', x='include'),
                  workspace='4gb', maxit= 20)
mod2_HT <- update.asreml(mod2_HT)
mod2_HT <- update.asreml(mod2_HT)
mod2_HT <- update.asreml(mod2_HT)

# GEBVs
GEBVs_mod2_HT <- summary(mod2_HT, coef = T)$coef.random |>
  as.data.frame() |>
  rownames_to_column() |>
  glimpse()

hist(GEBVs_mod2_HT$solution)

# ### Maturity----
# # model
# mod2_MAT <- asreml(fixed=predicted.value~1,
#                   random=~study+vm(germplasm, Ginv.sparse):fa(study,2),
#                   weights=weight,
#                   family=asr_gaussian(dispersion = 1),
#                   data=,
#                   na.action = na.method(y='omit', x='omit'),
#                   workspace='4gb', maxit= 20)
# mod2_MAT <- update.asreml(mod2_MAT)
# mod2_MAT <- update.asreml(mod2_MAT)
# mod2_MAT <- update.asreml(mod2_MAT)
# 
# # GEBVs
# GEBVs_mod2_MAT <- summary(mod2_MAT, coef = T)$coef.random |>
#   as.data.frame() |>
#   rownames_to_column() |>
#   glimpse()
# 
# hist(GEBVs_mod2_MAT$solution)

save.image('Data/Step4_ST-GBLUP_fa2.RData')

# Combine GEBVs ----
GEBVs_ST_GBLUP_fa2_l <- data.frame() |>
  bind_rows(GEBVs_mod2_YLD |> mutate(trait='grain_yield_bu_ac'),
            GEBVs_mod2_TW |> mutate(trait='test_weight_lb_bu'),
            GEBVs_mod2_HD |> mutate(trait='heading_time_jd'),
            GEBVs_mod2_HT |> mutate(trait='plant_height_in'),
#            GEBVs_mod2_MAT |> mutate(trait='maturity'),
            ) |>
  mutate(rowname=str_remove(rowname,'vm\\(germplasm, Ginv\\.sparse\\)_')) |>
  mutate(rowname=str_remove(rowname,'fa\\(study, 2\\)_')) |>
  filter(!str_detect(rowname, c('Comp1|Comp2'))) |>
  filter(!str_detect(rowname, 'study_')) |>
  separate(rowname, into = c('germplasm', 'study'), sep=':') |>
  glimpse()

# Save ----
save.image('Data/Step4_ST-GBLUP_fa2.RData')