#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 01-merge_datasets.R
#'
#' Version: 0.0.4
#'
#' Purpose of script: Generate cleaned master data file for proteomic/biomarkers
#'
#' Author: Robert B Butler III, Michael P. Brown
#'
#' Date Created: 2024-08-01
#'
#' Copyright (c) 2024
#' Email: rrbutler@stanford.edu, mpb@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'
#'
#' ---------------------------

#' load up the packages we will need:  (uncomment as required)
library(data.table)
library(plyr)
library(dplyr)
library(irr)
library(ggplot2)
library(ggpubr)
library(patchwork)
# cores <- parallel::detectCores()

# directory locations
geno_dir <- file.path('../genomic')

# Main --------------------------------------------------------------
# load up SomaLogic metadata file and fix misnamed samples
metafile <- '../dataProt_SS-205063.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP_ADRC_Feb2021_plasmaptauadded_Ted.csv'
prot_rename <- fread(
  paste('in2csv', 
    file.path(geno_dir, "Somalogic_renamed_samples.xlsx")
  )
)
a <- fread(metafile)
a0 <- a

# relabel mislabeled values
a$PIDN <- mapvalues(
  a$PIDN, 
  from = prot_rename[keep == TRUE, original_id], 
  to = prot_rename[keep == TRUE, updated_id]
)

# checking covariates
covs <- c(
  "Age", "Gender", "Storage_days", "Study", "Subarray", "Fed.Fasted", "APOE", 
  "HybControlNormScale", "PlateId", "Storage"
)
all(covs %in% names(a)) # missing Apoe until added

# checking biomarkers
bios <- c(
  "PlasmaPTau181", "PlasmaAB142P", "PlasmaAB140P", "PlasmaABRatio",
  "PlasmapTau217", "PlasmapTau217_AB42Ratio", "PlasmaGFAP", "PlasmaNfL"
)
all(bios %in% names(a)) # missing Apoe until added

qc_cols <- c(
  "PIDN", "SampleID", "ADRC_ID", "Visit", "Visit_Consensus", "ConnectivityZscore",
  "Diagnosis_group", "Date.of.draw"
)
all(qc_cols %in% names(a)) # missing Apoe until added

# protein columns
prots <- readLines("../prot_names.txt")

# dealing with duplicated names
drops <- which(duplicated(names(a)[!names(a) %in% prots]))
a[, (drops) := NULL]

# Initial QC filters
a = a[SampleType == "Sample"]
a = a[Study != "Control"]
a = a[RowCheck == "PASS"]

# sapply(a[, 1:20], table)

# filter out missing diagnosis, PD people
pd_folks = c(unique(a$Diagnosis_group[a$Diagnosis_group %like% "PD"]), "")
lbd_folks = c(unique(a$Diagnosis_group[a$Diagnosis_group %like% "LBD"]), "")
a = a[!Diagnosis_group %in% c(pd_folks, lbd_folks)]

# recode Subarray to letters to avoid confusing as quantitative
a[, Subarray := recode(Subarray, `1`="A", `2`="B", `3`="C", `4`="D", `5`="E", 
                       `6`="F", `7`="G", `8`="H")]

# checking duplicate samples ------------------------------------------------
# first duplicates by SampleID
cols = c(qc_cols, covs[-7])
a[duplicated(SampleID) | duplicated(SampleID, fromLast=T), ..cols]
# fwrite(
#   a[duplicated(SampleID) | duplicated(SampleID, fromLast=T), ..cols], 
#   file = "../20240801-duplicates_by_sampleID.csv"
# )

# same 14 duplicates regardless, check biomarker measures 
# mostly appear to be duplicates across studies or across plates, internal controls
# so select a random one for each sampleID
set.seed(123)
a = a[a[ , .I[sample(.N,1)] , by = SampleID]$V1]

# next duplicate paricipants at the same visit 
a[duplicated(a[, .(PIDN, Visit)]) | duplicated(a[, .(PIDN, Visit)], fromLast=T), ..cols]
# fwrite(
#   a[duplicated(a[, .(PIDN, Visit)]) | duplicated(a[, .(PIDN, Visit)], fromLast=T), ..cols], 
#   file = "../20240801-duplicates_by_PIDN-Visit.csv"
# )
# some are from renaming PIDNs for dual enrolls, will be fixed by Visit reordering

# 1 is a genuine duplicate with the same date of draw, Study, PlateID
a[duplicated(a[, .(PIDN, Visit, Date.of.draw)]) | duplicated(a[, .(PIDN, Visit, Date.of.draw)], fromLast=T), ..cols]
# fwrite(
#   a[duplicated(a[, .(PIDN, Visit, Date.of.draw)]) | duplicated(a[, .(PIDN, Visit, Date.of.draw)], fromLast=T), ..cols],
#   file = "../20240801-duplicates_by_PIDN-Visit.csv"
# )
# randomly drop one of the duplicates
set.seed(123)
a = a[a[ , .I[sample(.N,1)] , by =c("PIDN", "Visit", "Date.of.draw")]$V1]

# the remaining 5 have differing draw dates, suggesting their Visit and 
# Visit_Consensus are wrong. Convert to Date.of.draw to date field to sort
a[, Date.of.draw := as.Date(Date.of.draw, format="%m/%d/%Y")]
dups = a[duplicated(a[, .(PIDN, Visit)]) , PIDN]
a[PIDN %in% dups, ..cols][order(PIDN, Date.of.draw)]
# all are part of the "Wave 2" group from Wagner that have what looks like a 
# 6-month follow-up or something similar. We can convert to Days since first 
# visit to account for non-exact visit spacing. 

# First, reorder visits to be correct for all PIDNs
# a[order(PIDN, Date.of.draw), ..cols]
a[, Visit := frank(Date.of.draw, ties.method="first"), by=c("PIDN")]
# definitely misordered by Date.of.draw... 
# should update metadata, but in the mean time reorder by date
all(a[order(PIDN, Date.of.draw), Visit] == a[order(PIDN, Visit), Visit])

# converting timescales to days since baseline visit -------------------------
# convert Date.of.draw to Days since 1st visit
a[, DS1V := Date.of.draw - Date.of.draw[Visit == 1], by = .(PIDN)]
a[, DS1V := as.numeric(DS1V)]
qc_cols <- c(qc_cols, "DS1V")
cols = c(qc_cols, covs[-7])

# Appears there are a number of folks with date issues resulting in NA(n = 18)
a[is.na(DS1V), ..cols]
a[PIDN %in% unique(a[is.na(DS1V), PIDN]), ..cols]

# Of those, there are 4 with Visit_Consensus "Baseline", but NA draw dates
# Those are easily converted to DS1V = 0
a[is.na(Date.of.draw) & Visit_Consensus == "Baseline", DS1V := 0]

# Correct misgendered folks
# Check for folks without gender continuity
mg = a[, .N, by = .(PIDN, Gender)][, .N, by = PIDN][N > 1]$PIDN
lapply(1:length(mg), function(i){a[PIDN %like% mg[i], Gender]})
# Parse in file with gender corrections
look <- fread("../master_duplicate_list.csv")
look <- look[additional_field_to_correct != ""]
look[, Date.of.draw := as.Date(Date.of.draw, format="%m/%d/%Y")]
# Correct gender swaps
lapply(seq(nrow(look)), function(i){
  a[PIDN == look$PIDN[i] & 
    ADRC_ID == look$alias[i] & 
    Date.of.draw == look$Date.of.draw[i], 
    (look$additional_field_to_correct[i]) := look$corrected_value[i]
  ]
})
# Make sure genders corrected successfully
lapply(1:length(mg), function(i){a[PIDN %like% mg[i], Gender]})

# Merge Apoe tables ---------------------------------------------------------
# Add in Apoe 
b <- fread(
  paste('in2csv', 
    file.path(geno_dir, "Participant_List_APOE_Race_and_ancestry.xlsx")
  )
)
# relabel mislabeled values
b$ID <- mapvalues(
  b$ID, 
  from = prot_rename[keep == TRUE, original_id], 
  to = prot_rename[keep == TRUE, updated_id]
)

# check matches in ID and ADRC_ID/Linkcode to make sure they are PIDNs
# conflicts exist between Poston IDs and PIDNs for one case. Since none have genetic
# information or Apoes, can filter out. Future studies with Poston need to
# verify mapping of 6-digit codes to PIDN

# drop those without any APOE status
b <- b[!is.na(NCRAD_APOE) | !is.na(APOE_WGS)]
# verify PIDN number is present and not duplicated in b
b[!ID %in% a$PIDN, .(ID, `ADRC_ID/Linkcode`)]

# 11 duplicated IDs in the Apoe table
# randomly drop one of the duplicates
set.seed(123)
b <- b[b[ , .I[sample(.N,1)] , by = "ID"]$V1]


# change to factor
apo_cols = c("NCRAD_APOE", "APOE_WGS")
b[, (apo_cols) := lapply(.SD, function(i) {
  recode_factor(i,
    `22` = "e22",
    `23` = "e23",
    `33` = "e33",
    `34` = "e34",
    `44` = "e44",
    `24` = "e24"
  )
}), .SDcols = apo_cols]
# fill missing WGS_APOEs with NCRAD_APOE
b[, APOE := APOE_WGS]
b[is.na(APOE), APOE := NCRAD_APOE]
a <- merge(a,
  b[, .(ID, APOE)],
  by.x = "PIDN",
  by.y = "ID", 
  all.x = TRUE, 
  sort = FALSE
)
cols <- c(qc_cols, covs)
# table(b$APOE)
# table(a$APOE)
# table(a[DS1V == 0, APOE])

# Hardcoded add LRS_apoe types, recode existing APOE as factor
b <- fread("../genomic/lrs_apoe.txt", colClasses = c("character", "factor"))
a <- merge(a, b, by.x = "PIDN", by.y = "IID", all.x = TRUE, sort = FALSE)
#### May be useless since previous APOE that was merged was of type factor
a[, APOE := recode_factor(APOE, 
  `22` = "e22",
  `23` = "e23",
  `33` = "e33",
  `34` = "e34",
  `44` = "e44",
  `24` = "e24"
)]

# fill missing APOEs with LRS_APOE
a[is.na(APOE), APOE := LRS_APOE]
# # Remove LRS_APOE after filling APOE blanks
# a[, LRS_APOE := NULL]
# checking LRS call accuracy
table(a$APOE == a$LRS_APOE)
# looks like 384/400 simultaneous calls agree
irr::kappa2(a[, .(APOE, LRS_APOE)])
# Kappa = 0.94, z = 38, p = 0

# Merge AB Plasma ---------------------------
# Add in AB 
b <- fread(
  paste('in2csv', "Plasma_AB_Data.xlsx")
)
b <- b[!Study %like% "Poston"] # Drop Poston
b <- b[!Sample_ID %like% "lood"] # Drop No Blood Yet
b <- b[!is.na(PlasmaABRatio)] # Drop folks without data

# split column into 
b[, c("PIDN", "Date.of.draw") := tstrsplit(Sample_ID, "\\s(?=\\d{1,2}\\/)", perl = TRUE)]
# removes extra spaces around "+"
b[, PIDN := gsub(" +\\+ +", "+", PIDN)]

# trim off Y1/Wave1 etc. from PIDN
b[, PIDN := sub(" .*$", "", PIDN)]

# relabel mislabeled values
b$PIDN <- mapvalues(
  b$PIDN, 
  from = prot_rename[keep == TRUE, original_id], 
  to = prot_rename[keep == TRUE, updated_id]
)

# remove leading zeros from SampleID
b[, PIDN := gsub("^0+", "", PIDN)]

# preliminary cleaning of date entries, one entry with " ?" 
b[, Date.of.draw := sub(" .*$", "", Date.of.draw)]
# convert Date.of.draw in b to date format and remove NAs
b[, Date.of.draw := as.Date(Date.of.draw, format="%m/%d/%y")]
b <- b[!is.na(Date.of.draw)]
# oldest date of draw in proteomics
thresh <- max(a$Date.of.draw, na.rm = TRUE)
# trim dates after latest draw date in proteomics
b <- b[Date.of.draw <= thresh]

# Check matching on SampleID
# 400 did match, check the ones that didn't 
# look at ones that did not (minus Poston/No Blood)
b[, PIDN := trimws(PIDN)]

# check for duplicates
b[duplicated(b[, .(PIDN, Date.of.draw)]) | duplicated(b[, .(PIDN, Date.of.draw)], fromLast = TRUE)]

# check pluses
a[PIDN %like% "\\+", ..cols] 

# if any do match but cannot be caught programmatically, make 3 column list .csv to map values from
a <- merge(a, 
  b[, .(PIDN, Date.of.draw, PlasmaAB142P, PlasmaAB140P, PlasmaABRatio)], 
  by = c("PIDN", "Date.of.draw"), 
  all.x = TRUE, 
  sort = FALSE
)

sum(!is.na(a$PlasmaABRatio))
# 744 patients with matched AB data 

# Merge NFL/GFAP/PTau217 --------------------
b <- fread(
  paste('in2csv', "ADRC_Plasma_pTau217_GFAP_NLF_part2_06272024.xlsx")
)

# check visit ADRC_ID matches across visits
b[, names(b) %like% "_ADRCID", with = FALSE][1:100]

# split and melt to long format
vis_names <- paste0("P", 1:5)
vist_list <- lapply(vis_names, function(i) {
  cols <- names(b)[names(b) %like% i]
  x <- b[, ..cols]
  setnames(x, gsub(paste0(i, "_"), "", cols))
})
setattr(vist_list, 'names', vis_names)
# merge lists for each visit
c <- rbindlist(vist_list, idcol = "Visit")
# drop NA values for missing visits
c <- c[!is.na(ADRCID)]
# adjust to a uniform column type
c[, Date.of.draw := as.Date(Date_of_Draw, format="%m/%d/%y")]
# oldest date of draw in proteomics
thresh <- max(a$Date.of.draw, na.rm = TRUE)
# trim dates after latest draw date in proteomics
c <- c[Date.of.draw <= thresh]

# now repeat the process, match on ADRCID and Date.of.draw, 
# remember Date_of_draw is unformatted, but can reformat and drop later
c[Plasma_NFL_YN == TRUE & Plasma_NFL_within_range_YN == FALSE]
c[Plasma_pTau217_YN == TRUE & Plasma_pTau217_within_range_YN == FALSE]
c[Plasma_GFAP_YN == TRUE & Plasma_GFAP_within_range_YN == FALSE]
# 1 P1 GFAP not within range reported

# finally after checks merge each biomarker separately
a <- merge(a, 
  c[Plasma_NFL_within_range_YN == TRUE, .(ADRCID, Date.of.draw, Plasma_NFL)],
  by.x = c("ADRC_ID", "Date.of.draw"), 
  by.y = c("ADRCID", "Date.of.draw"), 
  all.x = TRUE, 
  sort = FALSE
)

a <- merge(a, 
  c[Plasma_GFAP_within_range_YN == TRUE, .(ADRCID, Date.of.draw, Plasma_GFAP)],
  by.x = c("ADRC_ID", "Date.of.draw"), 
  by.y = c("ADRCID", "Date.of.draw"), 
  all.x = TRUE, 
  sort = FALSE
)

a <- merge(a, 
  c[Plasma_pTau217_within_range_YN == TRUE, .(ADRCID, Date.of.draw, Plasma_pTau217)],
  by.x = c("ADRC_ID", "Date.of.draw"), 
  by.y = c("ADRCID", "Date.of.draw"),
  all.x = TRUE, 
  sort = FALSE
)

# Adding in Lumipulse pTau217, NfL, GFAP (colnames PlasmapTau217, PlasmaNfL, PlasmaGFAP) --------
# load in each table
b_names <- c(
  "ADRC_Lumipulse_PTau217_nfl_gfap.xlsx", 
  "SAMS_Lumipulse_PTau217_GFAP_NfL.xlsx"
)
b <- lapply(b_names, function(i) fread(paste('in2csv', i)))
setattr(b, 'names', c('ADRC', 'SAMS'))

# join tables
b <- rbindlist(b, idcol = "study", fill = TRUE)
# match character
b$Patient_ID <- as.character(b$Patient_ID)

# checking mismatches
b[!Patient_ID %in% a$PIDN]
b[!Patient_ID %in% a$PIDN & Blood_Draw_Date <= thresh]
# 199 are in the date range to have SomaLogic data, but missing
b[!Patient_ID %in% a0$PIDN & Blood_Draw_Date <= thresh]
# All but 2 are in a0, one of which was renamed and does match in a 
# The other is a sample from 2014 that isn't in Soma at all
# For the rest, since they are in a0 are missing genetic data 

# trim dates after latest draw date in proteomics
b <- b[Blood_Draw_Date <= thresh]

# merge with table
a <- merge(a, 
  b[, .(Patient_ID, Blood_Draw_Date, PlasmapTau217, PlasmaNfL, PlasmaGFAP)],
  by.x = c("PIDN", "Date.of.draw"), 
  by.y = c("Patient_ID", "Blood_Draw_Date"),
  all.x = TRUE, 
  sort = FALSE
)

# Merge Quaternix and Lumipulse ---------------------------------------

pairs <- list(
  c("PlasmaNfL", "Plasma_NFL"),
  c("PlasmaGFAP", "Plasma_GFAP"),
  c("PlasmapTau217", "Plasma_pTau217")
)
# check consistency
p_list <- lapply(pairs, function(i) {
  ggscatter(a, x = i[[2]], y = i[[1]],
     color = "black", size = 3, # Points color, shape and size
     add = "reg.line",  # Add regressin line
     add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
     conf.int = TRUE, # Add confidence interval
     cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
     cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
  ) +
  labs(x = paste("Quanterix Simoa", i[[2]]), y = paste("Lumipulse", i[[1]]))
})

ggsave(
  plot = wrap_plots(p_list),
  file = "platform_comparison_corr.png",
  width = 21
)
# pTau217 does not have a very robust correlation

# checking amount to be gained by adding samples
# > nrow(a[is.na(NfL) & !is.na(Plasma_NFL)][order(DS1V)][, head(.SD, 1), by = "PIDN"])
# [1] 54
# > nrow(a[is.na(GFAP) & !is.na(Plasma_GFAP)][order(DS1V)][, head(.SD, 1), by = "PIDN"])
# [1] 54
# > nrow(a[is.na(pTau217) & !is.na(Plasma_pTau217)][order(DS1V)][, head(.SD, 1), by = "PIDN"])
# [1] 51
# current totals
# > nrow(a[!is.na(NfL)][order(DS1V)][, head(.SD, 1), by = "PIDN"])
# [1] 465
# > nrow(a[!is.na(GFAP)][order(DS1V)][, head(.SD, 1), by = "PIDN"])
# [1] 465
# > nrow(a[!is.na(pTau217)][order(DS1V)][, head(.SD, 1), by = "PIDN"])
# [1] 465

# Create pTau217/AB142P ratio column ----------------------------------
a[, PlasmapTau217_AB42Ratio := PlasmapTau217 / PlasmaAB142P]

# Final checks --------------------------------------------------------
# spot check a with a0 for a few rows
# cols <- cols[c(1:8, 10:15, 17:19)]
# a0[PIDN %in% a$PIDN, ..cols]

# Write to output
out = file.path("../")
new_cols = c(qc_cols, covs, bios, prots)
fwrite(a[, ..new_cols], file = file.path(out, paste("dataProt_SS-205063.merged.QC.csv")))
