# Single trial analysis ----

# Objective ----
# - Run single trial BLUP models to check reliability
# - Run single trial BLUE models to estimate genotype's blues

rm(list=objects()) # clean workspace

# Packages ----
library(tidyverse) # R packages for data science
library(asreml) # ASReml-R package

# Load data ----
load('Data/YT_pheno_data.RData')

## Check dataset ----
YTpheno_l |>
  glimpse()


# Toy dataset

dat <- YTpheno_l |>
  filter(study=='Big6_Urb_23', trait=='grain_yield_bu_ac') |>
  glimpse()

# BLUP - single trial ----

# Single trial BLUP to calculate reliability

## ASReml-R options for spatial analysis

asreml.options(gammaPar= T) # gamma parameterization

## BLUP function ----
mod.blup_single_trial <- function(dat) {
  mod <- asreml(value~1+block,
                random = ~germplasm,
                residual = ~ar1v(row):ar1(col),
                na.action = na.method(y='include',x='include'),
                data=dat, maxit = 20)
  # Update model
  mod <- update.asreml(mod)
  mod <- update.asreml(mod)
  mod <- update.asreml(mod)
  mod <- update.asreml(mod)
  # Mean reliabilty
  m_rel<- predict(mod, classify='germplasm', ignore=c('(Intercept)', 'block'))$pvals |>
    as.data.frame() |>
    mutate(pev = std.error^2,
           Vg = summary(mod)$varcomp['germplasm','component'],
           rel = 1-(pev/Vg),
           m_rel=mean(rel)) |>
    summarise(m_rel=mean(rel))
  return(m_rel)
}

## Run BLUP ----

# Use map() to run model for each individual trial and return the trial mean reliability
reliability <- YTpheno_l |>
  group_by(study, trait) |>
  nest() |>
  mutate(r2=map(data,~mod.blup_single_trial(.x))) |>
  dplyr::select(-c(data)) |>
  unnest(r2) |>
  glimpse()

## Plot reliability ----
reliability |>
  ggplot(aes(x=study, y=m_rel)) +
  geom_point() +
  facet_wrap(~trait, scales = 'free_x') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# BLUE - single trial ----
# Individual genotype blues, SE, and weights to be used in other models

## BLUE function ----
mod.blue_single_trial <- function(dat) {
  mod <- asreml(value~1+block + germplasm,
                residual = ~ar1v(row):ar1(col),
                na.action = na.method(y='include',x='include'),
                data=dat, maxit = 20)
  # Update model
  mod <- update.asreml(mod)
  mod <- update.asreml(mod)
  mod <- update.asreml(mod)
  mod <- update.asreml(mod)
  # Predict blues
  blues <- predict(mod, classify='germplasm', pworkspace='1gb')$pvals
  # Result that will be printed
  return(blues)
}

## Run BLUE ----

# Use map() to run model for each individual trial and return individual genotype's blues, SE, and weight
blues <- YTpheno_l |>
  group_by(study,trait) |>
  nest() |>
  mutate(blues=map(data,~mod.blue_single_trial(.x))) |>
  dplyr::select(-data) |>
  unnest(blues) |>
  dplyr::select(-status) |>
  mutate(weight = 1/(std.error^2)) |>
  glimpse()


blues <- YTpheno_l |>
  group_by(study) |>
  summarise_at(vars(year, location, trial), ~unique(.x)) |>
  left_join(blues) |>
  filter(germplasm!='FILL') |>
  glimpse()


# Save blues ----
saveRDS(blues, file = 'Data/Single_Trial_Blues.RDS')