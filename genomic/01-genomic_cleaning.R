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
library(plyr)

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

# Read in sample metadata where available ----
prot_meta <- fread(
  gsub("\\", "", # needs cleaner approach, should be able to handle raw character constants, but just reads file into stdin
    file.path(
      rootdir, 
      'plasma_ptau_proteomics', 
      'dataProt_SS-205063.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP_ADRC_Feb2021_plasmaptauadded_Ted.csv'
    ),
    fixed = TRUE
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

# some PIDNs seem to have an additional suffix checking overlap with Poston
intersect(wgs_meta[Study == "Poston", `ADRC_ID/Linkcode`], wgs_meta[Study == "SAMS", ID])
# confirmed, those and only those are present in the overlap with Poston linking IDs
# those suffixes can be safely dropped in the fam file and liknk to the correct PIDN (below)

# Drop Poston from the metadata when transferring.
wgs_meta <- wgs_meta[Study != "Poston"]

# Read in Yann's sample renaming info for swaps/duplicates ----
# long read
lrs_rename <- fread(
  paste('in2csv', file.path(wk_dir, "LRS_renamed_samples.xlsx"))
)
lrs_rename <- lrs_rename[!is.na(original_id)]

# short read
srs_rename <- fread(
  paste('in2csv', file.path(wk_dir, "WGS_renamed_samples.xlsx"))
)
srs_rename <- srs_rename[!is.na(original_id)]

prot_rename <- fread(
  paste('in2csv', file.path(wk_dir, "Somalogic_renamed_samples.xlsx"))
)

# read in kinship for dup checks ----
# long read
lrs_kin <- fread(
  file.path(
    lr_dir, 
    'genesis_processing', 
    'Stanford_ADRC_LRS_QC_for_PCA_relatedpairs.txt'
  )
)
lrs_kin[, ID1 := gsub("ADRC", "", ID1)]
lrs_kin[, ID2 := gsub("ADRC", "", ID2)]

# short read
srs_kin <- fread(
  file.path(
    sr_dir, 
    'genesis_processing', 
    'Stanford_ADRC_WGS_QC_for_PCA_relatedpairs.txt'
  )
)
srs_kin[, ID1 := sub("^0+", "", sub("^\\w+_", "", ID1))]
srs_kin[, ID2 := sub("^0+", "", sub("^\\w+_", "", ID2))]

# Updating Somalogic metadata PIDNs to address joined IDs
prot_meta$PIDN <- mapvalues(
  prot_meta$PIDN, 
  from = prot_rename[keep == TRUE, original_id], 
  to = prot_rename[keep == TRUE, updated_id]
)
id_map <- unique(prot_meta[, .(ADRC_ID, PIDN)])[!is.na(ADRC_ID)]

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
lrs_fam$V2 <- mapvalues(
  lrs_fam$V2, 
  from = lrs_rename[keep == TRUE, original_id], 
  to = lrs_rename[keep == TRUE, updated_id]
)
# the _bis IIDs are all duplicated, so add to drop list
drops <- c(drops, lrs_fam[V2 %like% "_bis", V2])

# also drop the recommended drop samples where keep == FALSE
drops <- c(drops, lrs_rename[keep == FALSE, original_id])

# Checking for duplicates using kinship ----
# noticed some of the samples might be duplicates by dual enrollment (ADRC_ID + PIDN)
lrs_kin[order(-kin)][kin > 0.25]
# updated in the renaming file, can rerun without finding them

# check remaining ADRC_ID samples in kin
lrs_kin[ID1 %in% lrs_fam[V2 %like% 3900, V2]]
# none with notable matches

# check remaining ADRC_ID samples in prot_meta to make sure they aren't duplicated
prot_meta[ADRC_ID %in% lrs_fam[V2 %like% 3900, V2]]
lrs_fam[V2 %in% prot_meta[ADRC_ID %in% lrs_fam[V2 %like% 3900, V2], PIDN]]
# all are new and not duplicated, lift over the new names to the FAM file
lrs_fam$V2 <- mapvalues(
  lrs_fam$V2, 
  from = id_map$ADRC_ID, 
  to = id_map$PIDN,
  warn_missing = FALSE
)

# drop unnecessary FIDs
lrs_fam$V1 <- 0

# write fam names to "lrs_temp.fam"
fwrite(lrs_fam, 
  file = file.path(lr_dir, "lrs_temp.fam"), 
  sep = " ", 
  col.names = FALSE
)
# copy other plink files to lrs_temp prefix
file.copy(
  from = file.path(lr_dir, "LRS_merged_rsids.bed"),
  to = file.path(lr_dir, "lrs_temp.bed"),
  overwrite = FALSE
)
file.copy(
  from = file.path(lr_dir, "LRS_merged_rsids.bim"),
  to = file.path(lr_dir, "lrs_temp.bim"),
  overwrite = FALSE
)

# write drop to lrs_temp_remove.txt for plink to filter them
# fwrite(list(drops), file = "lrs_remove.txt")


# Matching SRS IDs ---------------------------------------------------
# Looks like mostly ADRC_IDs and some PIDNs for SAMS
# suffixes mentioned above can be safely dropped in the fam file and link to the correct PIDN
# move alphabetic prefix to FIDs
srs_fam[, c("V1", "V2") := tstrsplit(V2, "_", keep = c(1:2))]
srs_fam[, V2 := sub("^0+", "", V2)]

# Checking for duplicates using kinship ----
# noticed some of the samples might be duplicates by dual enrollment (ADRC_ID + PIDN)
srs_kin[order(-kin)][kin > 0.25]
# updated in the renaming file, can rerun without finding them
# all three are sample duplicates in the LRS as well, updated the renaming file

# Assume the numbering is PIDN, but ask Yann for confirmation
# renaming based on Yann's rename file for sample 3 swaps/4 duplicates  (plus my additions).
srs_fam$V2 <- mapvalues(
  srs_fam$V2, 
  from = srs_rename[keep == TRUE, original_id], 
  to = srs_rename[keep == TRUE, updated_id]
)
# the _bis IIDs are all duplicated, so add to drop list
srs_drops <- srs_fam[V2 %like% "_bis", V2]

# also drop the recommended drop samples where keep == FALSE
srs_drops <- c(srs_drops, srs_rename[keep == FALSE, original_id])
srs_drops <- srs_drops[srs_drops %in% srs_fam$V2] # some are not in the fam file, may be from previous steps

# check remaining ADRC_ID samples in prot_meta to make sure they aren't duplicated
prot_meta[ADRC_ID %in% srs_fam[V2 %like% 3900, V2]]
srs_fam[V2 %in% prot_meta[ADRC_ID %in% srs_fam[V2 %like% 3900, V2], PIDN]]
# all are new and not duplicated, lift over the new names to the FAM file
srs_fam$V2 <- mapvalues(
  srs_fam$V2, 
  from = id_map$ADRC_ID, 
  to = id_map$PIDN,
  warn_missing = FALSE
)

# drop unnecessary FIDs
srs_fam$V1 <- 0

# write fam names to "srs_temp.fam"
fwrite(srs_fam, 
  file = file.path(sr_dir, "srs_temp.fam"), 
  sep = " ", 
  col.names = FALSE
)
# copy other plink files to lrs_temp prefix
file.copy(
  from = file.path(sr_dir, "Stanford_ADRC_WGS_QC_step3.bed"),
  to = file.path(sr_dir, "srs_temp.bed"),
  overwrite = FALSE
)
file.copy(
  from = file.path(sr_dir, "Stanford_ADRC_WGS_QC_step3.bim"),
  to = file.path(sr_dir, "srs_temp.bim"),
  overwrite = FALSE
)

# prefer SRS over LRS for biallelic SNP accuracy, so add the overlap 
# (without srs_drops) to the drops list
drops <- unique(c(setdiff(intersect(srs_fam$V2, lrs_fam$V2), srs_drops), drops))

# write drop to srs_temp_remove.txt for plink to filter them
fwrite(list(srs_drops), file = file.path(sr_dir, "srs_remove.txt"))

# write drop to lrs_temp_remove.txt for plink to filter them
fwrite(list(drops), file = file.path(lr_dir, "lrs_remove.txt"))



