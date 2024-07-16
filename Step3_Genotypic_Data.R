# Genotypic data ----

# Objective ----
# - Subset phenotypic and genotypic data
# - Create relationship matrix


rm(list=objects()) # clean workspace

# Packages ----
library(tidyverse) # R packages for data science
library(gaston) # Genetic Data Handling
#library(rrBLUP)    # R package for Ridge Regression and Best Linear Unbiased Prediction
#library(Matrix)    # R package for sparse and dense matrix classes and methods
library(ASRgenomics) # R package for Genomic Selection Analysis in R

# Load data ----
## BLUES ----
blues <- readRDS('data/Single_Trial_Blues.RDS')

geno <- read.vcf('Data/24NorGrain_23And24Big6_all_regions_Filtered.vcf.gz')

#Correct the line names
geno@ped$id <- sub('^.*:', '', geno@ped$id)
yrseries<- gsub('20', '', as.character(2000:2019))

for(i in 1:length(yrseries)){
  geno@ped$id<- gsub(paste('IL', yrseries[i], sep=''), yrseries[i], geno@ped$id)
}
geno@ped$id<- gsub('IL2020', '2020', geno@ped$id)
geno@ped$id<- gsub('IL21', 'IL2021', geno@ped$id)
geno@ped$id<- gsub('16LCSDH', 'IL16LCSDH', geno@ped$id)
geno@ped$id<-gsub('PIO-25R74', 'Pio25R74', geno@ped$id)
geno@ped$id<-gsub('KASKASKIA', 'Kaskaskia', geno@ped$id)

# Get markers for the genotypes in the Big6 trials
geno <- select.inds(geno, id %in% blues$germplasm)

# Filter snps with maf > 0.05 to include only polymorphic snps
#geno <- select.snps(geno, maf > 0.05)

# Create molecular matrix
genoM <- as.matrix(geno)

# Filtering molecular matrix
genoM_filter <- qc.filtering(M=genoM, maf=0.05, marker.callrate = 0.2,
                             ind.callrate = 0.33, impute = F, plots = T)
dim(genoM_filter$M.clean)
genoM_filter$M.clean[1:5,1:5]

# Calculate the G matrix
G <- G.matrix(M=genoM_filter$M.clean, method = 'VanRaden', na.string = NA)$G

G[1:8,1:8]

# Beding G matrix
G.bend <- G.tuneup(G=G, bend = T, eig.tol = 1e-03)$Gb
G.bend[1:8,1:8]

# Inverse of G matrix
Ginv <- G.inverse(G=G.bend)$Ginv
Ginv[1:8,1:8]
det(Ginv)

# Generate the three-column sparse form of the matrix Ginv
Ginv.sparse <- G.inverse(G=G.bend, sparseform = T)$Ginv
str(Ginv.sparse)

# Filter phenotypic data from individuals with marker data
blues_filter <- blues |>
  filter(germplasm %in% attributes(Ginv.sparse)$colNames) |>
  glimpse()

# Check if Ginv and phenotypic data matches

match.kinship2pheno(K=Ginv, pheno.data = blues,
                    indiv = 'germplasm_name',
                    clean = F, mism = F)

# Save
save(blues_filter, Ginv, Ginv.sparse, file='Data/Pheno_and_Gmat.RData')
