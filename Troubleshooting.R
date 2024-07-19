# Troubleshooting

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

# GEBVs
#GEBVs_mod0_gy <- predict(mod0_gy, classify='germplasm', pworkspace='2gb')$pval

#hist_mod0_gy <- hist(GEBVs_mod0_gy$predicted.value)

save.image('Data/troubleshooting.RData')

### Grain yield diag()----
# model
mod1_gy <- asreml(fixed=predicted.value~1,
                  random=~study+vm(germplasm, Ginv.sparse):diag(study),
                  weights=weight,
                  family=asr_gaussian(dispersion = 1),
                  data=blues_filter |> filter(trait == 'grain_yield_bu_ac'),
                  na.action = na.method(y='include', x='include'),
                  workspace='4gb', maxit= 20)
mod1_gy <- update.asreml(mod1_gy)
mod1_gy <- update.asreml(mod1_gy)


### Grain yield fa2----
# model
mod2_gy <- asreml(fixed=predicted.value~1,
                  random=~study+vm(germplasm, Ginv.sparse):fa(study,2),
                  weights=weight,
                  family=asr_gaussian(dispersion = 1),
                  data=blues_filter |> filter(trait == 'grain_yield_bu_ac'),
                  na.action = na.method(y='include', x='include'),
                  workspace='8gb', maxit= 20)
mod2_gy <- update.asreml(mod2_gy)
mod2_gy <- update.asreml(mod2_gy)

# # GEBVs
# GEBVs_mod1_gy <- predict(mod1_gy, classify='germplasm:study', pworkspace='2gb')$pval
# 
# hist_mod1_gy <- hist(GEBVs_mod1_gy$predicted.value)

save.image('Data/troubleshooting.RData')

load('Data/troubleshooting.RData')


blues_filter <- blues_filter |> droplevels()

# Check variance components ----
summary(mod1_gy)$varcomp |>
  as.data.frame() |>
  rownames_to_column()

# Calculate genetic correlation among groups ----
vc <- summary(mod1_gy)$varcomp # Variance components
R <- vc$component[2:18] # Variance
L1 <- vc$component[19:35] # Loading 1st order
L2 <- vc$component[36:52] # Loading 2nd order
L <- cbind(L1, L2) # Combine loadings from all three orders

VG <- L%*%t(L) + diag(R) # Variance covariance matrix
rownames(VG) <- levels(blues_filter$study) 
colnames(VG) <- levels(blues_filter$study)
CG <- cov2cor(VG) # Correlation matrix

CG

# Define color palette
my_colors <- colorRampPalette(c("#FF5F0F", "white", "#13294b"))(100)
library(corrplot)
# Plot
#png('Figures/gen_corr.png', width = 10, height = 10, units = 'in', res = 320)
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
#dev.off()

corrplot(CG)

heatmap(CG, symm = T)

ASRgenomics::kinship.heatmap(K=CG, dendrogram = T, col.label = T)


# Get the names of the rows and columns
row_names <- rownames(CG)
col_names <- colnames(CG)

# Select the rows and columns that contain '_23'
selected_rows <- grep("_23", row_names, value = TRUE)
selected_cols <- grep("_23", col_names, value = TRUE)

# Subset the data frame or matrix
CG_selected <- CG[selected_rows, selected_cols]

ASRgenomics::kinship.heatmap(K=CG_selected, dendrogram = T, col.label = T)


load('Data/troubleshooting.RData')

