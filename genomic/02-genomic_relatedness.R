#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 02-genomic_relatedness.R
##
## Version: 0.0.3
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
library(plyr)
library(dplyr)
library(future.apply)
cores <- parallel::detectCores()

# directory locations
lr_dir <- file.path('/mnt', 'box', 'Stanford_ADRC', 'LRS_WGS')
sr_dir <- file.path('/mnt', 'box', 'Stanford_ADRC', 'SRS_WGS')


# Functions ---------------------------------------------------------


# Main --------------------------------------------------------------
# load up SomaLogic metadata file and fix misnamed samples
metafile <- '../dataProt_SS-205063.merged.QC.csv'
prot_rename <- fread(
  paste('in2csv', "Somalogic_renamed_samples.xlsx")
)

# set up metadata file for relatedness subsets
a <- fread(metafile)


# # relabel mislabeled values
# a$PIDN <- mapvalues(
#   a$PIDN, 
#   from = prot_rename[keep == TRUE, original_id], 
#   to = prot_rename[keep == TRUE, updated_id]
# )

# checking covariates
covs <- c(
  "Age", "Gender", "Storage_days", "Study", "Subarray", "Fed.Fasted", "APOE", 
  "HybControlNormScale", "PlateId", "Storage"
)
all(covs %in% names(a)) # missing Apoe until added

# checking biomarkers
bios <- c(
  "PlasmaPTau181", "PlasmaAB142P", "PlasmaAB140P", "PlasmaABRatio",
  "Plasma_pTau217", "Plasma_pTau217_AB42Ratio", "Plasma_GFAP", "Plasma_NFL"
)
all(bios %in% names(a)) # missing Apoe until added

qc_cols <- c(
  "PIDN", "SampleID", "ADRC_ID", "Visit", "Visit_Consensus", "ConnectivityZscore",
  "Diagnosis_group", "Date.of.draw", "DS1V"
)
all(qc_cols %in% names(a)) # missing Apoe until added

# protein columns
prots = readLines("../prot_names.txt")

# variables
cols = c(qc_cols, bios, covs)


# # # # exploring the data (Comment out normally) --------------------------------
# # # # generate presence/absence matrix for upset
# # # library(UpSetR)
# # # b <- a[, lapply(.SD, \(i) ifelse(is.na(i), 0, 1)), .SDcols = bios]
# # # upset(b, nsets = 8, nintersects = NA, order.by = "freq")
# # # c <- a[DS1V == 0, lapply(.SD, \(i) ifelse(is.na(i), 0, 1)), .SDcols = bios]
# # # upset(c, nsets = 8, nintersects = NA, order.by = "freq")
# # # d <- a[DS1V == 0 & !is.na(APOE), lapply(.SD, \(i) ifelse(is.na(i), 0, 1)), .SDcols = bios]
# # # upset(d, nsets = 8, nintersects = NA, order.by = "freq")


# load plink file -------------------------------------------------------
if (! file.exists(file.path(sr_dir, "srs-lrs_comb.filt.gds"))) {
  gds_fn <- snpgdsBED2GDS(
    bed.fn = file.path(sr_dir, "srs-lrs_comb.filt.bed"), 
    bim.fn = file.path(sr_dir, "srs-lrs_comb.filt.bim"), 
    fam.fn = file.path(sr_dir, "srs-lrs_comb.filt.fam"), 
    out.gdsfn = file.path(sr_dir, "srs-lrs_comb.filt.gds")
  )
  # # map gds samples to PIDNs
  # gds <- GdsGenotypeReader(gds_fn)
  # sample.sel <- getScanID(gds)
  # a[, gds := PIDN %in% sample.sel ]
  # setcolorder(a, c(names(a)[!names(a) %in% cols], cols))
  # fwrite(a, file = gsub("_plasmaptauadded_Ted", "", metafile))
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
close(gds)

# PC-Relate -------------------------------
plan("multicore", workers = 8)
options(future.globals.maxSize = 4 * 1e9)
# prune to the relevant biomarker overlap
gds <- GenotypeBlockIterator(GenotypeData(GdsGenotypeReader(gds_fn, allow.fork = TRUE)))
sample.sel <- getScanID(gds)

# for each biomarker make grm (second loop includes the APOE covariate)
future_lapply(bios, function(i) {
  # keep only first visit where plasma available
  x <- a[!is.na(get(i))][order(DS1V)][, head(.SD, 1), by = "PIDN"]
  # drop samples with incomplete cases
  cols <- c("SampleID", i, covs)
  x <- x[complete.cases(x[, ..cols])]

  # intersection for variable of interest
  both <- intersect(x[, PIDN], sample.sel)

  # run PC-Relate
  mypcrelate <- pcrelate(gds, 
    pcs = mypcair$vectors[, 1:8], 
    training.set = intersect(mypcair$unrels, both),
    BPPARAM = BiocParallel::SerialParam(),
    sample.include = both
  )

  # write pcrelate RDS and GRM to file
  saveRDS(mypcrelate, file = paste("mypcrelate", i, "rds", sep = "."))
}, future.seed = TRUE)

# # to read to Matrix format in lme4 ----------------------------------------
# pcr <- readRDS(file = "mypcrelate.rds")
# grm <- pcrelateToMatrix(pcr, thresh = 2^(-11/2), scaleKin = 2)




# # # for each biomarker make grm (second loop includes the APOE covariate)
# # lapply(c("all", "APOE"), function(j) {
  # # future_lapply(bios, function(i) {
    # # # keep only first visit where plasma available
    # # x <- a[!is.na(get(i))][order(DS1V)][, head(.SD, 1), by = "PIDN"]
    # # # drop samples with incomplete cases
    # # if (j == "all") {
      # # cols <- c("SampleID", i, covs[-7])
    # # } else {
      # # cols <- c("SampleID", i, covs)
    # # }
    # # x <- x[complete.cases(x[, ..cols])]

    # # # intersection for variable of interest
    # # both <- intersect(x[, PIDN], sample.sel)

    # # # run PC-Relate
    # # mypcrelate <- pcrelate(gds, 
      # # pcs = mypcair$vectors[, 1:8], 
      # # training.set = intersect(mypcair$unrels, both),
      # # BPPARAM = BiocParallel::SerialParam(),
      # # sample.include = both
    # # )

    # # # write pcrelate RDS and GRM to file
    # # saveRDS(mypcrelate, file = paste("mypcrelate", i, j, "rds", sep = "."))
  # # }, future.seed = TRUE)
# # })


