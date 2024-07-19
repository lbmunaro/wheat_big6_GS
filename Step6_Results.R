# Results from ST-GBLUP & MT-GBLUP ----

# Objective ----
# - Get GEBVs
# - Estimate genetic correlations

rm(list=ls()) # clean workspace

# Packages ----

library(tidyverse) # R packages for data science
library(corrplot) # A visual exploratory tool on correlation matrix

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

# ST-GBLUP ----
GEBVs_ST_GBLUP_fa2_l |>
  select(study, trait, germplasm, solution) |>
  rename(modfa2_ST_gebvs=solution) |>
  mutate_at(vars(study:germplasm),~as.factor(.x)) |>
  left_join(metadata) |>
  left_join(eval_years) |>
  select(group, trait, study, trial, year, location, germplasm, eval_years, modfa2_ST_gebvs) |>
  glimpse()


# MT-GBLUP ----
## diag() ----
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

## Compare diag() & fa(,3) ----
mod0_MT_gebvs |>
  left_join(modfa3_MT_gebvs) |>
  ggplot(aes(x=mod0_MT_gebvs,y=modfa3_MT_gebvs)) +
  geom_point()
cor(mod0_MT_gebvs$mod0_MT_gebvs,modfa3_MT_gebvs$modfa3_MT_gebvs)


# Genetic correlations ----

## YLD among trials ----

#ST


## Among groups ----

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

# Change row and column names

# Define replacements
replacements <- c("grain_yield_bu_ac" = "YLD",
                  "heading_time_jd" = "HD",
                  'maturity_jd' = 'MAT',
                  "plant_height_in" = "PH",
                  "test_weight_lb_bu" = "TW")

# Apply replacements to row names
for (pattern in names(replacements)) {
  rownames(CG) <- gsub(pattern, replacements[pattern], rownames(CG))
}

# Apply replacements to column names
for (pattern in names(replacements)) {
  colnames(CG) <- gsub(pattern, replacements[pattern], colnames(CG))
}

CG |> glimpse()

## Plot correlation matrix ----

YLD_MT_CG <- CG[grep("YLD_", rownames(CG)), grep("YLD_", colnames(CG))]

heatmap(YLD_MT_CG)

# Define color palette
my_colors <- colorRampPalette(c("#FF5F0F", "white", "#13294b"))(100)

# Plot
#png('Figures/gen_corr.png', width = 10, height = 10, units = 'in', res = 320)
corrplot(YLD_MT_CG, 
         method = 'number', 
         type = 'lower', 
#         col = my_colors,
         number.cex = 1.2, 
         number.digits = 2,
         diag = T,
         tl.col = 'black', 
         tl.cex = 0.75,
         tl.srt = 10)
#dev.off()


load('Data/ST-GBLUP_fa2.RData')

summary(mod2_YLD)$varcomp |>
  as.data.frame() |>
  rownames_to_column()

# Calculate genetic correlation among groups ----
vc <- summary(mod2_YLD)$varcomp # Variance components
R <- vc$component[2:18] # Variance
L1 <- vc$component[19:35] # Loading 1st order
L2 <- vc$component[36:52] # Loading 2nd order
L <- cbind(L1, L2) # Combine loadings from all three orders

VG <- L%*%t(L) + diag(R) # Variance covariance matrix
#rownames(VG) <- levels(metadata$study) 
#colnames(VG) <- levels(blues_filter$study)
CG <- cov2cor(VG) # Correlation matrix

CG
heatmap(CG)
