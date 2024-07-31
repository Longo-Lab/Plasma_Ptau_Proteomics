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
library(plyr)
library(dplyr)
cores <- parallel::detectCores()

# directory locations
lr_dir <- file.path('/mnt', 'box', 'Stanford_ADRC', 'LRS_WGS')
sr_dir <- file.path('/mnt', 'box', 'Stanford_ADRC', 'SRS_WGS')


# Functions ---------------------------------------------------------


# Main --------------------------------------------------------------
# load up SomaLogic metadata file and fix misnamed samples
metafile <- '../michael_dataProt_SS-205063.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP_ADRC_Feb2021_plasmaptauadded_Ted.csv'
prot_rename <- fread(
  paste('in2csv', "Somalogic_renamed_samples.xlsx")
)

a <- fread(metafile)

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
  "Plasma_pTau217", "Plasma_pTau217_AB42Ratio", "Plasma_GFAP", "Plasma_NFL"
)
all(bios %in% names(a)) # missing Apoe until added

# protein columns
first_prot = "CRYBB2.10000.28.3"
last_prot = "IRF6.9999.1.3"
cols = c(names(a)[grep(first_prot, names(a)):grep(last_prot, names(a))])

# dealing with duplicated names
drops <- which(duplicated(names(a)[!names(a) %in% cols]))
a[, (drops) := NULL]

# # # # exploring the data (Comment out normally) --------------------------------
# # # setcolorder(a, c(names(a)[!names(a) %in% cols], cols))
# # # a[is.na(APOE) & `Genetics?` == 1, c(1:7,66:75)]
# # # length(unique(a[is.na(APOE) & `Genetics?` == 1, PIDN]))

# # # # downsample to just DS1V == 0 (from 01-sample_QC.R)
# # # a = a[SampleType == "Sample"]
# # # a = a[RowCheck == "PASS"]
# # # # filter out missing diagnosis, PD people
# # # pd_folks = c(unique(a$Diagnosis_group[a$Diagnosis_group %like% "PD"]), "")
# # # lbd_folks = c(unique(a$Diagnosis_group[a$Diagnosis_group %like% "LBD"]), "")
# # # a = a[!Diagnosis_group %in% c(pd_folks, lbd_folks)]
# # # a[, Date.of.draw := as.Date(Date.of.draw, format="%m/%d/%Y")]
# # # a[, DS1V := Date.of.draw - Date.of.draw[Visit == 1], .(PIDN)]
# # # a[, DS1V := as.numeric(DS1V)]

# # # # Of those, there are 4 with Visit_Consensus "Baseline", but NA draw dates
# # # # Those are easily converted to DS1V = 0
# # # a[is.na(Date.of.draw) & Visit_Consensus == "Baseline", DS1V := 0]

# # # # 7 remain that apppear to be visits without a baseline. 
# # # # for this analysis, I will calculate their DS1V with Visit 2
# # # a[is.na(DS1V), DS1V := Date.of.draw - Date.of.draw[Visit == 2], .(PIDN)]
# # # a <- a[PIDN != ""]

# # # # generate presence/absence matrix for upset
# # # library(UpSetR)
# # # b <- a[, lapply(.SD, \(i) ifelse(is.na(i), 0, 1)), .SDcols = bios]
# # # upset(b, nsets = 8, nintersects = NA, order.by = "freq")
# # # c <- a[DS1V == 0, lapply(.SD, \(i) ifelse(is.na(i), 0, 1)), .SDcols = bios]
# # # upset(c, nsets = 8, nintersects = NA, order.by = "freq")
# # # d <- a[DS1V == 0 & !is.na(APOE), lapply(.SD, \(i) ifelse(is.na(i), 0, 1)), .SDcols = bios]
# # # upset(d, nsets = 8, nintersects = NA, order.by = "freq")

# Prep Apoe data for analysis --------------------------------------------
apoe <- fread(file.path(lr_dir, "lrs_temp.raw"))

# the _bis IIDs are all duplicated, plus three singlet ids
apoe <- apoe[!IID %like% "[A-Za-z]"]

# also drop the recommended drop samples where keep == FALSE
lrs_rename <- fread(
  paste('in2csv', "LRS_renamed_samples.xlsx")
)
lrs_rename <- lrs_rename[!is.na(original_id)]
apoe <- apoe[!IID %in% lrs_rename[keep == FALSE, original_id]]

# generate apoe factor types
# columns are the concat of Apoe4(19:44908684:T:C) and Apoe2(19:44908822:C:T)
apoe[['LRS_APOE']] <- recode_factor(
  paste0(apoe$`19:44908684:T:C_T`, apoe$`19:44908822:C:T_C`), 
  `20` = "e22",
  `21` = "e23",
  `22` = "e33",
  `12` = "e34",
  `02` = "e44",
  `11` = "e24"
)
fwrite(apoe[, .(IID, LRS_APOE)], file = "lrs_apoe.txt")

# # # checking apoe types against NCRAD calls
# # b <- merge(
  # # unique(a[, .(PIDN, APOE)]), 
  # # apoe[, .(IID, LRS_APOE)], 
  # # by.x = "PIDN", 
  # # by.y = "IID"
# # )
# # b[['APOE']] <- recode_factor(
  # # b[['APOE']], 
  # # `22` = "e22",
  # # `23` = "e23",
  # # `33` = "e33",
  # # `34` = "e34",
  # # `44` = "e44",
  # # `24` = "e24"
# # )
# # table(b$APOE == b$LRS_APOE)
# # # looks like 384/400 simultaneous calls agree
# # irr::kappa2(b[, .(APOE, LRS_APOE)])
# # # Kappa = 0.94, z = 27.3, p = 0

# load plink file -------------------------------------------------------
if (! file.exists(file.path(sr_dir, "srs-lrs_comb.filt.gds"))) {
  gds_fn <- snpgdsBED2GDS(
    bed.fn = file.path(sr_dir, "srs-lrs_comb.filt.bed"), 
    bim.fn = file.path(sr_dir, "srs-lrs_comb.filt.bim"), 
    fam.fn = file.path(sr_dir, "srs-lrs_comb.filt.fam"), 
    out.gdsfn = file.path(sr_dir, "srs-lrs_comb.filt.gds")
  )
  # map gds samples to PIDNs
  gds <- GdsGenotypeReader(gds_fn)
  sample.sel <- getScanID(gds)
  a[, gds := PIDN %in% sample.sel ]
  setcolorder(a, c(names(a)[!names(a) %in% cols], cols))
  fwrite(a, file = gsub("_plasmaptauadded_Ted", "", metafile))
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
# prune to the relevant biomarker overlap
subfile <- tempfile()
gds <- GdsGenotypeReader(gds_fn)
####lapply for each biomarker#######
# intersection for variable of interest
gdsSubset(gdsfile, subfile, sample.include=sample.sel, snp.include=snp.sel)
# save as subfile gds

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

