#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 02-genomic_relatedness.R
##
## Version: 0.0.1
##
## Purpose of script: Generate cleaned fam files for LRS and WGS datasets
##
## Author: Robert B Butler III
##
## Date Created: 2024-07-16
##
## Copyright (c) 2024
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##
##   see https://github.com/Longo-Lab/Plasma_Ptau_Proteomics/issues/1
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(SNPRelate)
library(GENESIS)
library(GWASTools)
library(Matrix)
library(BiocParallel)
library(data.table)
cores <- parallel::detectCores()

# directory locations
lr_dir <- file.path('/mnt', 'box', 'Stanford_ADRC', 'LRS_WGS')
sr_dir <- file.path('/mnt', 'box', 'Stanford_ADRC', 'SRS_WGS')


# Functions ---------------------------------------------------------


# Main --------------------------------------------------------------
# load plink file
if (! file.exists(file.path(sr_dir, "srs-lrs_comb.filt.gds"))) {
  gds_fn <- snpgdsBED2GDS(
    bed.fn = file.path(sr_dir, "srs-lrs_comb.filt.bed"), 
    bim.fn = file.path(sr_dir, "srs-lrs_comb.filt.bim"), 
    fam.fn = file.path(sr_dir, "srs-lrs_comb.filt.fam"), 
    out.gdsfn = file.path(sr_dir, "srs-lrs_comb.filt.gds")
  )
} else {
  gds_fn <- file.path(sr_dir, "srs-lrs_comb.filt.gds")
}

# kinship mat
gds <- snpgdsOpen(gds_fn)
KINGmat <- kingToMatrix(snpgdsIBDKING(gds, num.thread = cores))
snpgdsClose(gds)

# run PC-AiR on pruned SNPs
gds <- GenotypeData(GdsGenotypeReader(gds_fn))
mypcair <- pcair(gds, kinobj = KINGmat, divobj = KINGmat, num.cores = cores)

# Make PC plot
png(file = "mypcair.elbow.png")
plot(1:length(mypcair$values), mypcair$values, ylab = "Eigenvalues", xlab = "PCs")
dev.off()

# save PC-AiR rds and csv of 8 PCs
saveRDS(mypcair, file = "mypcair.rds")
dt1 <- data.table(mypcair$vectors[, 1:8], keep.rownames = TRUE)
setnames(dt1, c("PIDN", paste0("PC", 1:8)))
fwrite(dt1, file = "mypcair.PCs.csv")


# PC-Relate -------------------------------
# run PC-Relate
gds <- GenotypeBlockIterator(gds)
mypcrelate <- pcrelate(gds, 
  pcs = mypcair$vectors[,1:8], 
  training.set = mypcair$unrels,
  BPPARAM=BiocParallel::SerialParam()
)

# generate genetic relationship matrix # maybe generate on the fly in the lme4 script.
grm <- pcrelateToMatrix(mypcrelate, thresh = 2^(-11/2), scaleKin = 2)

# write pcrelate RDS and GRM to file
saveRDS(mypcrelate, file = "mypcrelate.rds")

# # to read to Matrix format in lme4 ----------------------------------------
# pcr <- readRDS(file = "mypcrelate.rds")
# grm <- pcrelateToMatrix(pcr, thresh = 2^(-11/2), scaleKin = 2)

