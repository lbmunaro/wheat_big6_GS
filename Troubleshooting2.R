# troubleshooting2

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
blues_filter <- blues_filter |>
  filter(!study%in%c('Big6_Vin_23','Big6_Ith_23','Big6_Neo_24')) |>
  arrange(study, germplasm) |>
  droplevels() |>
  filter(!is.na(predicted.value)) |>
  glimpse()

blues_filter |>
  filter(is.na(predicted.value))

blues_filter |>
  ggplot(aes(x=germplasm,y=study)) +
  geom_tile()


# ST-GBLUP ----
blues_filter |> filter(trait == 'grain_yield_bu_ac') |> glimpse()

### Grain yield ----
# pheno_model
pheno_mod0_gy <- asreml(fixed=predicted.value~1,
                  random=~germplasm + study,
                  weights=weight,
                  family=asr_gaussian(dispersion = 1),
                  data=blues_filter |> filter(trait == 'grain_yield_bu_ac'),
                  na.action = na.method(y='include', x='include'),
                  workspace='4gb')
pheno_mod0_gy <- update.asreml(pheno_mod0_gy)
pheno_mod0_gy <- update.asreml(pheno_mod0_gy)

# GEBVs
GEBVs_pheno_mod0_gy <- predict(pheno_mod0_gy, classify='germplasm', pworkspace='2gb')$pval

hist_pheno_mod0_gy <- hist(GEBVs_pheno_mod0_gy$predicted.value)

save.image('Data/troubleshooting2.RData')

### Grain yield fa2----
blues_filter |> glimpse()
# pheno_model
pheno_mod1_gy <- asreml(fixed=predicted.value~1,
                  random=~study+fa(study,2):germplasm,
                  weights=weight,
                  family=asr_gaussian(dispersion = 1),
                  data=blues_filter |> filter(trait == 'grain_yield_bu_ac'),
                  na.action = na.method(y='include', x='include'),
                  workspace='8gb')
pheno_mod1_gy <- update.asreml(pheno_mod1_gy)
pheno_mod1_gy <- update.asreml(pheno_mod1_gy)

# # GEBVs
 GEBVs_pheno_mod1_gy <- predict(pheno_mod1_gy, classify='germplasm:study', pworkspace='2gb')$pval
 
 hist_pheno_mod1_gy <- hist(GEBVs_pheno_mod1_gy$predicted.value)

save.image('Data/troubleshooting2.RData')

### Grain yield corgh2----
# pheno_model
pheno_mod2_gy <- asreml(fixed=predicted.value~1,
                  random=~study+germplasm:corgh(study),
                  weights=weight,
                  family=asr_gaussian(dispersion = 1),
                  data=blues_filter |> filter(trait == 'grain_yield_bu_ac'),
                  na.action = na.method(y='include', x='include'),
                  workspace='4gb', maxit= 20)
pheno_mod2_gy <- update.asreml(pheno_mod2_gy)
pheno_mod2_gy <- update.asreml(pheno_mod2_gy)

# # GEBVs
# GEBVs_pheno_mod2_gy <- predict(pheno_mod2_gy, classify='germplasm', pworkspace='2gb')$pval
# 
# hist_pheno_mod2_gy <- hist(GEBVs_pheno_mod2_gy$predicted.value)

save.image('Data/troubleshooting2.RData')

### Grain yield diag()----
# pheno_model
pheno_mod3_gy <- asreml(fixed=predicted.value~1,
                  random=~study+germplasm:diag(study),
                  weights=weight,
                  family=asr_gaussian(dispersion = 1),
                  data=blues_filter |> filter(trait == 'grain_yield_bu_ac'),
                  na.action = na.method(y='include', x='include'),
                  workspace='4gb', maxit= 20)
pheno_mod3_gy <- update.asreml(pheno_mod3_gy)
pheno_mod3_gy <- update.asreml(pheno_mod3_gy)

save.image('Data/troubleshooting2.RData')


