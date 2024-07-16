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
blues_filter |> filter(trait == 'grain_yield_bu_ac') |> glimpse()

### Grain yield ----
# model
mod0_gy <- asreml(fixed=predicted.value~1,
                  random=~vm(germplasm, Ginv.sparse)+study,
                  weights=weight,
                  family=asr_gaussian(dispersion = 1),
                  data=blues_filter |> filter(trait == 'grain_yield_bu_ac'),
                  na.action = na.method(y='include', x='include'),
                  workspace='4gb', maxit= 20)
mod0_gy <- update.asreml(mod0_gy)
mod0_gy <- update.asreml(mod0_gy)
mod0_gy <- update.asreml(mod0_gy)

# GEBVs
GEBVs_mod0_gy <- predict(mod0_gy, classify='germplasm', pworkspace='2gb')$pval

hist(GEBVs_mod0_gy$predicted.value)

### Test weight----
# model
mod0_tw <- asreml(fixed=predicted.value~1,
                  random=~vm(germplasm, Ginv.sparse)+study,
                  weights=weight,
                  family=asr_gaussian(dispersion = 1),
                  data=blues_filter |> filter(trait == 'test_weight_lb_bu'),
                  na.action = na.method(y='include', x='include'),
                  workspace='4gb', maxit= 20)
mod0_tw <- update.asreml(mod0_tw)
mod0_tw <- update.asreml(mod0_tw)
mod0_tw <- update.asreml(mod0_tw)

# GEBVs
GEBVs_mod0_tw <- predict(mod0_tw, classify='germplasm', pworkspace='2gb')$pval

hist(GEBVs_mod0_tw$predicted.value)

### Heading time----
# model
mod0_ht <- asreml(fixed=predicted.value~1,
                  random=~vm(germplasm, Ginv.sparse)+study,
                  weights=weight,
                  family=asr_gaussian(dispersion = 1),
                  data=blues_filter |> filter(trait == 'heading_time_jd'),
                  na.action = na.method(y='include', x='include'),
                  workspace='4gb', maxit= 20)
mod0_ht <- update.asreml(mod0_ht)
mod0_ht <- update.asreml(mod0_ht)
mod0_ht <- update.asreml(mod0_ht)

# GEBVs
GEBVs_mod0_ht <- predict(mod0_ht, classify='germplasm', pworkspace='2gb')$pval

hist(GEBVs_mod0_ht$predicted.value)

### Plant height----
# model
mod0_ph <- asreml(fixed=predicted.value~1,
                  random=~vm(germplasm, Ginv.sparse)+study,
                  weights=weight,
                  family=asr_gaussian(dispersion = 1),
                  data=blues_filter |> filter(trait == 'plant_height_in'),
                  na.action = na.method(y='include', x='include'),
                  workspace='4gb', maxit= 20)
mod0_ph <- update.asreml(mod0_ph)
mod0_ph <- update.asreml(mod0_ph)
mod0_ph <- update.asreml(mod0_ph)

# GEBVs
GEBVs_mod0_ph <- predict(mod0_ph, classify='germplasm', pworkspace='2gb')$pval

hist(GEBVs_mod0_ph$predicted.value)

### Maturity----
# model
mod0_mt <- asreml(fixed=predicted.value~1,
                  random=~vm(germplasm, Ginv.sparse)+study,
                  weights=weight,
                  family=asr_gaussian(dispersion = 1),
                  data=blues_filter |> filter(trait == 'maturity_jd'),
                  na.action = na.method(y='include', x='include'),
                  workspace='4gb', maxit= 20)
mod0_mt <- update.asreml(mod0_mt)
mod0_mt <- update.asreml(mod0_mt)
mod0_mt <- update.asreml(mod0_mt)

# GEBVs
GEBVs_mod0_mt <- predict(mod0_mt, classify='germplasm', pworkspace='2gb')$pval


hist(GEBVs_mod0_mt$predicted.value)

# Combine GEBVs in a single dataset
GEBVs_ST_GBLUP <- data.frame() |>
  bind_rows(GEBVs_mod0_gy |> mutate(trait='grain_yield'),
            GEBVs_mod0_tw |> mutate(trait='test_weight'),
            GEBVs_mod0_ht |> mutate(trait='heading_time'),
            GEBVs_mod0_ph |> mutate(trait='plant_height'),
            GEBVs_mod0_mt |> mutate(trait='maturity')) |>
  glimpse()

write.csv(GEBVs_ST_GBLUP,'Data/ST-GBLUP_gebvs_long.csv')

write.csv(
  GEBVs_ST_GBLUP |>
    select(-status) |>
    pivot_wider(names_from = trait, values_from = c(predicted.value, std.error)) |>
    arrange(desc(predicted.value_grain_yield )) |>
    glimpse(),
  'Data/ST-GBLUP_gebvs_wide.csv')

save.image('Data/ST-GBLUP.RData')