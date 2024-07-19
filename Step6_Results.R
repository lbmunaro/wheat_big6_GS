# Results from ST-GBLUP & MT-GBLUP ----

# Objective ----
# - Get GEBVs
# - Estimate genetic correlations

rm(list=ls()) # clean workspace

# Packages ----

library(tidyverse) # R packages for data science
library(corrplot) # A visual exploratory tool on correlation matrix
library(gplots) # Various R Programming Tools for Plotting Data

# Load data ----
load('Data/Step4_ST-GBLUP_fa2.RData') # ST-GBLUP
load('Data/Step5_MT-GBLUP_fa3.RData') # MT-GBLUP

# Metadata ----
# Trial information
metadata <- blues_group |>
  select(group,study,year,location,trial,trait) |>
  group_by(group) |>
  summarise_at(vars(study:trait),unique) |>
  glimpse()

# Germplasm evaluation years
eval_years <- blues_group |>
  group_by(germplasm) |>
  mutate(eval_years = paste(sort(unique(year)), collapse = "&")) |>
  summarise(eval_years = unique(eval_years)) |>
  glimpse()

# Number of study each genotype was evaluated within year
load('Data/Big6_pheno_data.RData')
core <- YTpheno_l |>
  filter(trait=='grain_yield_bu_ac') |>
  group_by(year, germplasm) |>
  summarise(n=n_distinct(study)) |>
  pivot_wider(names_from = year, values_from = n, names_prefix = 'n.study_') |>
  mutate_at(vars(2,3),~as.factor(.x)) |>
  glimpse()

# Names replacement
replacements <- c("grain_yield_bu_ac" = "YLD",
                  "heading_time_jd" = "HD",
                  'maturity_jd' = 'MAT',
                  "plant_height_in" = "HT",
                  "test_weight_lb_bu" = "TW")

# ST-GBLUP ----
results_ST_mod_fa2_l <- GEBVs_ST_GBLUP_fa2_l |>
  select(study, trait, germplasm, solution) |>
  rename(modfa2_ST_gebvs=solution) |>
  mutate_at(vars(study:germplasm),~as.factor(.x)) |>
  left_join(metadata) |>
  left_join(eval_years) |>
  select(group, trait, study, trial, year, location, germplasm, eval_years, modfa2_ST_gebvs) |>
  glimpse()

results_ST_mod_fa2_l |>
  ggplot(aes(x=modfa2_ST_gebvs)) +
  geom_histogram() +
  facet_wrap(~trait, scales = 'free')

GEBVs_ST_GBLUP_fa2_w <- results_ST_mod_fa2_l |>
  select(group, eval_years, germplasm, modfa2_ST_gebvs) |>
  mutate(group=str_replace_all(group, replacements)) |>
  pivot_wider(names_from = group, values_from = modfa2_ST_gebvs) |>
  rowwise() |>
  mutate(YLD_mean = mean(c_across(starts_with('YLD_'))),
         TW_mean = mean(c_across(starts_with('TW_'))),
         HD_mean = mean(c_across(starts_with('HD_'))),
         HT_mean = mean(c_across(starts_with('HT_')))) |>
  arrange(desc(YLD_mean)) |>
  relocate(YLD_mean,.before = 'YLD_Big6_Fra_23') |>
  relocate(TW_mean, .before = 'TW_Big6_Fra_23') |>
  relocate(HD_mean, .before = 'HD_Big6_Laf_23') |>
  relocate(HT_mean, .before = 'HT_Big6_Laf_24') |>
  left_join(core) |>
  glimpse()
write.csv(GEBVs_ST_GBLUP_fa2_w, 'Results/gebvs_ST-GBLUP_fa2.csv')

# MT-GBLUP ----
## fa(,3) ----
modfa3_MT_gebvs <- summary(modfa3_MT.GBLUP, coef=T)$coef.random |>
  as.data.frame() |>
  rownames_to_column() |>
  mutate(rowname=str_remove(rowname,'fa\\(group, 3\\)_')) |>
  mutate(rowname=str_remove(rowname,'vm\\(germplasm, Ginv\\.sparse\\)_')) |>
  mutate(rowname=str_remove(rowname,'group_')) |>
  filter(!str_detect(rowname, c('Comp1|Comp2|Comp3'))) |>
  separate(rowname, into = c('group', 'germplasm'), sep=':') |>
  select(group,germplasm,solution) |>
  rename(modfa3_MT_gebvs=solution) |>
  mutate_at(vars(group,germplasm),~as.factor(.x)) |>
  left_join(metadata) |>
  left_join(eval_years) |>
  select(group, trait, study, trial, year, location, germplasm, eval_years, modfa3_MT_gebvs) |>
  glimpse()

GEBVs_MT_GBLUP_fa3_w <- modfa3_MT_gebvs |>
  select(group, eval_years, germplasm, modfa3_MT_gebvs) |>
  mutate(group=str_replace_all(group, replacements)) |>
  pivot_wider(names_from = group, values_from = modfa3_MT_gebvs) |>
  rowwise() |>
  mutate(YLD_mean = mean(c_across(starts_with('YLD_'))),
         TW_mean = mean(c_across(starts_with('TW_'))),
         HD_mean = mean(c_across(starts_with('HD_'))),
         HT_mean = mean(c_across(starts_with('HT_'))),
         MAT_mean = mean(c_across(starts_with('MAT_')))) |>
  arrange(desc(YLD_mean)) |>
  relocate(YLD_mean,.before = 'YLD_Big6_Fra_23') |>
  relocate(TW_mean, .before = 'TW_Big6_Fra_23') |>
  relocate(HD_mean, .before = 'HD_Big6_Laf_23') |>
  relocate(HT_mean, .before = 'HT_Big6_Laf_24') |>
  relocate(MAT_mean, .before = 'MAT_Big6_Neo_24') |>
  left_join(core) |>
  glimpse()
write.csv(GEBVs_MT_GBLUP_fa3_w, 'Results/gebvs_MT-GBLUP_fa3.csv')

# Compare models ----

## diag() vs. fa(,3) ----

# diag()
mod0_MT_gebvs <- summary(mod0_MT.GBLUP, coef=T)$coef.random |>
  as.data.frame() |>
  rownames_to_column() |>
  mutate(rowname=str_remove(rowname,'vm\\(germplasm, Ginv\\.sparse\\)_')) |>
  mutate(rowname=str_remove(rowname,'group_')) |>
  separate(rowname, into = c('group', 'germplasm'), sep=':') |>
  select(group,germplasm,solution) |>
  rename(mod0_MT_gebvs=solution) |>
  mutate_at(vars(group,germplasm),~as.factor(.x)) |>
  left_join(metadata) |>
  left_join(eval_years) |>
  select(group, trait, study, trial, year, location, germplasm, eval_years, mod0_MT_gebvs) |>
  glimpse()

# Plot
mod0_MT_gebvs |>
  left_join(modfa3_MT_gebvs) |>
  ggplot(aes(x=mod0_MT_gebvs,y=modfa3_MT_gebvs)) +
  geom_point()
# Correlation
cor(mod0_MT_gebvs$mod0_MT_gebvs,modfa3_MT_gebvs$modfa3_MT_gebvs)

## ST fa(2) & MT fa(3) ----
STvsMT <- results_ST_mod_fa2_l |>
  select(group, germplasm, modfa2_ST_gebvs) |>
  left_join(modfa3_MT_gebvs) |>
  glimpse()
# Plot
ggplot(STvsMT,aes(x=modfa2_ST_gebvs,y=modfa3_MT_gebvs)) +
  geom_point()
# Correlation
cor(STvsMT$modfa2_ST_gebvs,STvsMT$modfa3_MT_gebvs)

# Genetic correlations ----

## YLD among trials ----

# ST
summary(mod2_YLD)$varcomp |>
  as.data.frame() |>
  rownames_to_column()

# Calculate genetic correlation among trials
vc <- summary(mod2_YLD)$varcomp # Variance components
R <- vc$component[2:18] # Variance
L1 <- vc$component[19:35] # Loading 1st order
L2 <- vc$component[36:52] # Loading 2nd order
L <- cbind(L1, L2) # Combine loadings from all three orders

VG <- L%*%t(L) + diag(R) # Variance covariance matrix
rownames(VG) <- levels(blues_filter$study)
colnames(VG) <- levels(blues_filter$study)
CG_ST <- cov2cor(VG) # Correlation matrix

CG_ST |> glimpse()

### Plot ----
png('Results/gen_corr_YLD_ST-GBLUP_fa2.png', width = 16, height = 16, units = 'in', res = 320)
par(mar = c(10, 4, 4, 2) + 0.1)  # Increase the bottom margin to avoid cutting off labels
heatmap.2(CG_ST, trace = "none", margins = c(10, 10))
dev.off()

## Groups as combination of trait  & trials ----

# Check variance components
summary(modfa3_MT.GBLUP)$varcomp |>
  as.data.frame() |>
  rownames_to_column()

# Calculate genetic correlation among groups ----
vc <- summary(modfa3_MT.GBLUP)$varcomp # Variance components
R <- vc$component[1:53] # Variance
L1 <- vc$component[54:106] # Loadings 1st order
L2 <- vc$component[107:159] # Loading 2nd order
L3 <- vc$component[160:212] # Loading 3rd order
L <- cbind(L1, L2, L3) # Combine loadings from all three orders

VG <- L%*%t(L) + diag(R) # Variance covariance matrix
rownames(VG) <- levels(blues_group$group) 
colnames(VG) <- levels(blues_group$group)
CG <- cov2cor(VG) # Correlation matrix

CG |> glimpse()

### Plot ----
png('Results/gen_corr_groups_MT-GBLUP_fa3.png', width = 16, height = 16, units = 'in', res = 320)
par(mar = c(10, 4, 4, 2) + 0.1)  # Increase the bottom margin to avoid cutting off labels
heatmap.2(CG, trace = "none", margins = c(10, 10))
dev.off()