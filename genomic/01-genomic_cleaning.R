#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 01-genomic_cleaning.R
##
## Version: 0.0.1
##
## Purpose of script: Generate cleaned fam files for LRS and WGS datasets
##
## Author: Robert B Butler III
##
## Date Created: 2024-06-27
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
library(data.table)
library(stringi)

# directory locations
rootdir <- r"(/mnt/box/Robert\ Butler\'s\ Files)"
wk_dir <- file.path(rootdir, 'plasma_ptau_proteomics', 'genomic')
lr_dir <- file.path('/mnt', 'box', 'Stanford_ADRC', 'LRS_WGS')
sr_dir <- file.path('/mnt', 'box', 'Stanford_ADRC', 'SRS_WGS')

# variables
nameset <- "combined_seq"
first_prot <- "CRYBB2.10000.28.3"
# ensver <- "111"

# Functions ---------------------------------------------------------
# Main --------------------------------------------------------------

# Read in reference fam files
lrs_fam <- fread(file.path(lr_dir, "LRS_merged_rsids.fam"))
srs_fam <- fread(file.path(sr_dir, "Stanford_ADRC_WGS_QC_step3.fam"))

# Read in sample metadata where available
prot_meta <- fread(
  file.path(
    rootdir, 
    'plasma_ptau_proteomics', 
    'dataProt_SS-205063.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP_ADRC_Feb2021_plasmaptauadded_Ted.csv'
  )
)
# just fetch the metadata to make it easier
cols = c(names(prot_meta)[1:(grep(first_prot, names(prot_meta)) - 1)])
prot_meta <- prot_meta[, ..cols]

# requires system/conda install of csvkit "conda install csvkit"
wgs_meta <- fread(
  paste(
    'in2csv', 
    file.path(sr_dir, "Participant_List_APOE_Race_and_ancestry.xlsx")
  )
)


# Read in Yann's sample renaming info for swaps/duplicates
lrs_rename <- fread(
  paste('in2csv', file.path(wk_dir, "LRS_renamed_samples.xlsx"))
)
lrs_rename <- lrs_rename[!is.na(original_id)]
srs_rename <- fread(
  paste('in2csv', file.path(wk_dir, "WGS_renamed_samples.xlsx"))
)
srs_rename <- lrs_rename[!is.na(original_id)]

# read in kinship for dup checks
lrs_kin <- fread(
  file.path(
    lr_dir, 
    'genesis_processing', 
    'Stanford_ADRC_LRS_QC_for_PCA_relatedpairs.txt'
  )
)
lrs_kin[, ID1 := gsub("ADRC", "", ID1)]
lrs_kin[, ID2 := gsub("ADRC", "", ID2)]

srs_kin <-

# Matching LRS IDs ---------------------------------------------------
# check LRS
# lrs_fam[, .N, by=V1]
# ADRC FIDs are duplicates of IIDs, drop the numbers
lrs_fam$V1 <- gsub("ADRC\\d+", "ADRC", lrs_fam$V1, perl = TRUE)
lrs_fam$V2 <- gsub("ADRC", "", lrs_fam$V2)
# Looks like 35 by ADRC ID, 434 by PBMC, 137 SAMS, three singlet ids

# three singlet ids don't map to any available metadata and will be dropped in plink
drops <- lrs_fam[!V1 %in% c("ADRC", "PBMC", "SAMS"), V2]

# Assume the numbering is PIDN, but ask Yann for confirmation
# renaming based on Yann's rename file for sample 3 swaps/4 duplicates  (plus my additions).
lrs_fam$V2 <- stri_replace_all_fixed(
  lrs_fam$V2, 
  pattern = lrs_rename[keep == TRUE, original_id], 
  replacement = lrs_rename[keep == TRUE, updated_id], 
  vectorize_all = FALSE
)
# the _bis IIDs are all duplicated, so add to drop list
drops <- c(drops, lrs_fam[V2 %like% "_bis", V2])

# also drop the recommended drop samples where keep == FALSE
drops <- c(drops, lrs_rename[keep == FALSE, original_id])

# Checking for duplicates using kinship --------------------------------
# noticed some of the samples might be duplicates by dual enrollment (ADRC_ID + PIDN)
lrs_kin[order(-kin)][kin > 0.25]
# updated in the renaming file, can rerun without finding them