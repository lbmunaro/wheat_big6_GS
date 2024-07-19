# Results from Multi-trait GBLUP ----

# Objective ----
# - Predict GEBVs

rm(list=ls()) # clean workspace

# Packages ----

library(tidyverse) # R packages for data science
library(corrplot) # A visual exploratory tool on correlation matrix

# Load data ----
load('Data/Step6_GEBVs_MT-GBLUP_fa3.RData')

# Check variance components ----
summary(mod0_MT.GBLUP)$varcomp |>
  as.data.frame() |>
  rownames_to_column()

# Calculate genetic correlation among groups ----
vc <- summary(mod0_MT.GBLUP)$varcomp # Variance components
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

# Define color palette
my_colors <- colorRampPalette(c("#FF5F0F", "white", "#13294b"))(100)

# Plot
png('Figures/gen_corr.png', width = 10, height = 10, units = 'in', res = 320)
corrplot(CG, 
         method = 'circle', 
         type = 'lower', 
         col = my_colors,
         number.cex = 1.2, 
         number.digits = 2,
         diag = T,
         tl.col = 'black', 
         tl.cex = 0.75,
         tl.srt = 10)
dev.off()

# GEBVs ----

# Metadata
metadata <- blues_group |>
  select(group,study,year,location,trial,trait) |>
  group_by(group) |>
  summarise_at(vars(study:trait),unique) |>
  glimpse()

eval_years <- blues_group |>
  group_by(germplasm) |>
  mutate(eval_years = paste(sort(unique(year)), collapse = "&")) |>
  summarise(eval_years = unique(eval_years)) |>
  glimpse()

# Long format
GEBVs_MT_long <- GEBVs_MT |>
  left_join(metadata) |>
  left_join(eval_years) |>
  glimpse()


summary(mod0_MT.GBLUP)$varcomp

GEBVs_MT_long |> filter(!study%in%c('Big6_Vin_23','Big6_Ith_23','Big6_Neo_24','Big6_Laf_24')) |>
  filter(trait=='heading_time_jd') |>
  arrange(predicted.value) |>
  View()

ggplot(GEBVs_MT_long |> filter(!study%in%c('Big6_Vin_23','Big6_Ith_23','Big6_Neo_24','Big6_Laf_24'))) +
  geom_histogram(aes(x=predicted.value)) +
  facet_wrap(~trait, scales = 'free', ncol = 2) +
  theme_bw()

ggplot(GEBVs_MT_long) +
  geom_histogram(aes(x=predicted.value)) +
  facet_wrap(~trait, scales = 'free', ncol = 2) +
  theme_bw()

ggplot(GEBVs_MT_long) +
  geom_jitter(aes(y=predicted.value, x=year)) +
  facet_wrap(~trait, scales = 'free', ncol = 2) +
  theme_bw()

GEBVs_MT_long |>
  filter(study=='Big6_Ith_23') |>
  ggplot() +
  geom_histogram(aes(x=predicted.value)) +
  facet_wrap(~trait, scales = 'free', ncol = 2) +
  theme_bw()

# Wide format
GEBVs_MT_wide <- GEBVs_MT |>
  select(group, germplasm, predicted.value) |>
  pivot_wider(names_from = group, values_from = predicted.value) |>
  rowwise() |>
  mutate(YLD_mean = mean(c_across(starts_with('grain_yield'))),
         HD_mean = mean(c_across(starts_with('heading_time'))),
         MAT_mean = mean(c_across(starts_with('maturity'))),
         PH_mean = mean(c_across(starts_with('plant_height'))),
         TW_mean = mean(c_across(starts_with('test_weight')))) |>
  arrange(desc(YLD_mean)) |>
  glimpse()

write.csv(GEBVs_MT_wide, file = 'Data/MT-GBLUP_gebvs_wide.csv')


# Comparing 

coef <- summary(mod0_MT.GBLUP, coef = T)$coef.random |>
  as.data.frame() |>
  rownames_to_column() |>
  separate(rowname, into = c('group', 'germplasm'), sep=':') |>
  mutate(group=str_remove(group,'fa\\(group, 3\\)_')) |>
  mutate(germplasm=str_remove(germplasm,'vm\\(germplasm, Ginv\\.sparse\\)_')) |>
  filter(!str_detect(group, c('Comp1|Comp2|Comp3'))) |>
  select(group,germplasm,solution) |>
#  left_join(GEBVs_MT_long) |>
  glimpse()

# Plot
ggplot(coef,aes(x=predicted.value,y=solution)) +
  geom_point(mapping = aes(x = solution, y = predicted.value))


coef |>
  ggplot() +
  geom_histogram(aes(x=solution)) +
  facet_wrap(~trait, scales = 'free', ncol = 2) +
  theme_bw()

hist(coef$predicted.value)

load('Data/Big6_pheno_data.RData')

pheno <- YTpheno_l |>
  group_by(study, trait, germplasm) |>
  summarise(value=mean(value)) |>
  filter(!is.na(value)) |>
  glimpse()

pheno |>
  left_join(coef) |>
  ggplot(aes(x=value,y=predicted.value)) +
  geom_point() +
  facet_wrap(~trait, scales='free')

pheno |>
  left_join(coef) |>
  ggplot(aes(x=value,y=solution)) +
  geom_point()

MT <- coef |>
  filter(str_detect(group,'grain_yield')) |>
  mutate(study=str_remove(group,'grain_yield_bu_ac_')) |>
  arrange(study,germplasm) |>
  glimpse()

ST <- summary(mod1_gy, coef = T)$coef.random |>
  as.data.frame() |>
  rownames_to_column() |>
  mutate(rowname=str_remove(rowname,'vm\\(germplasm, Ginv\\.sparse\\)_')) |>
  mutate(rowname=str_remove(rowname,'fa\\(study, 2\\)_')) |>
  filter(!str_detect(rowname, c('Comp1|Comp2'))) |>
  filter(!str_detect(rowname, 'study_')) |>
  separate(rowname, into = c('germplasm', 'study'), sep=':') |>
  arrange(study,germplasm) |>
  glimpse()

plot(x=ST$solution, y=MT$solution)

