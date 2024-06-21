#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 02-linreg_glm.R
##
## Version: 0.0.2
##
## Purpose of script: calculate correlation of plasma ptau with SomaLogic plasma
## proteomic levels
##
## Author: Edward Wilson, Robert B Butler III
##
## Date Created: 2023-08-25
##
## Copyright (c) 2023
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##   Linear regression glm for ALL VISITS with corrections 
##   for fixed effects: filter covariates based on importance
##   
##   PlateRunDate & ScannerID couldn't be fit to the model, so dropping those
##   
##   This takes a while, so batch unless you like waiting
##   Usage:
##     sbatch -J ptau --mem=12G -c 16 -t 06:00:00 -p interactive \
##       -o %A_02-linreg_lme4.log \
##       --wrap "ml R/4.2.2; Rscript 02-linreg_glm.R"
##
##   REMEMBER, WHEN COMMENTING IN CODE, DO NOT REFERENCE INDIVIDUAL INFORMATION:
##   SAMPLES/PIDNS, DATES OR ANY INDIVIDUAL LEVEL INFO
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(data.table)
# library(lmerTest) # loads lme4
setDTthreads(percent=75)

# directory locations
rootdir = '/mnt/c/Users/rrbutler/Documents/SRCF_shutdown'
wk_dir = paste0(rootdir, '/plasma_ptau_proteomics/all_visit')
setwd(wk_dir)
refdir = paste0(rootdir, '/reference/biodomains')

# variables
nameset = "dataProt_SS-205063.pTau.all_visit"
first_prot = "CRYBB2.10000.28.3"
last_prot = "IRF6.9999.1.3"


# Functions ----------------------------------------------------------------
#' Generalized linear model based on a flexible formula
#'
#' @param dat Data.table with the original data
#' @param outcome String specifying the column in dat containing the outcome variable 
#' @param exposure String specifying the column in dat containing the exposure variable
#' @param covs Vector of column names containing covariates to be included in the formula
#'
#' @returns A named list: 'Estimate' (beta), 'Std. Error' (se), 't value' (test statistic), 'Pr(>|t|)' (p-value), 'nobs' (number of observations)
get_linreg = function(dat, outcome, exposure, covs=c()){
  model = glm(paste(outcome, "~", paste(c(exposure, covs), collapse=" + ")),
              na.action=na.exclude, data=dat)
  a = summary(model)
  return(as.list(c(a$coefficients[2, ], nobs=nobs(model))))
}

#' Recodes dummy categorical variables and can subset by vector of categories
#'
#' @param dat Data.table with the original data
#' @param dcol String representing the column to be split
#' @param keep Vector of strings that are the categories within dcol to extract. default is all c()
#'
#' @returns modified data.table with new dummy variable columns named "dcol.category" after original column
recode_dummies = function(dat, dcol, keep=c()){
  dt = dat
  # check input
  if (length(keep) == 0){
    inds <- unique(dt[[dcol]])
  } else if (length(keep) > 0 & all(keep %in% unique(dt[[dcol]]))){
    inds <- unique(dt[[dcol]])[unique(dt[[dcol]]) %in% keep]
  } else {
    stop("Keep must be a vector of names in the original column")
  }
  # split out dummies
  dt[, c(paste(dcol, inds, sep=".")) := lapply(inds, function(x) ifelse(get(dcol) == x, 1, 0))]
  # reorder behind original variable
  v = grep(paste0("^", dcol, "$"), names(dt))[1]
  inds = paste(dcol, inds, sep=".")
  setcolorder(dt, c(names(dt)[1:v], inds, names(dt)[(v + 1):(ncol(dt) - length(inds))]))
  return(dt)
}


# Main -------------------------
dat = fread(paste0(nameset, '.INTnorm.csv.gz'))
# dat$Diagnosis_group = factor(dat$Diagnosis_group, c("HC", "SCI", "MCI", "MCI-AD", "AD", "LBD"))
names(dat) = make.names(names(dat))

# variables to drop from the model (based on feature importance): 
# top 30 ----------------------------------------------------
# add in single category columns with feature importance

# Study was important, drop ADRC to force it as the intercept
dat = recode_dummies(dat, "Study", keep=c("Poston", "Wagner", "Kerchner"))

# PlateId was important, but not all plates were top30 in either glmcross-validated or plr estimates
dat = recode_dummies(dat, "PlateId", keep=c("P0028449", "P0028439", "P0028451", 
                                            "P0028565", "P0028437", "P0028448",
                                            "P0028463", "P0028430"))

# best 2of3 (cv, pls, rf; either rf-impur or rf-perm stat can count as only one)
drops = c("Diagnosis_group", "Visit", "DS1V", "SampleID", "PlasmaPTau181", 
          "Storage_days", "Study", "PlateId", "Subarray", "HybControlNormScale",
          "Fed.Fasted", "Storage", "Time.To.Ship", "sv1", "sv5", "sv12", "sv13",
          "sv15", "sv18", "sv19", "sv20", "sv23", "sv24", "sv25", "sv29", "sv31",
          "sv32", "sv33", "sv34", "sv35", "sv36", "sv37", "sv38", "sv39", "sv40",
          "sv42", "sv43", "sv44", "sv45", "sv46", "sv48", "sv49", "sv50", "sv51",
          "sv52", "sv53", "sv54", "sv55", "sv56", "sv57", "sv58", "sv59", "sv60",
          "sv61", "sv62", "sv63", "sv64", "sv65")

# top 50 ----------------------------------------------------
# # add in single category columns with feature importance

# # Study was important, drop ADRC to force it as the intercept
# dat = recode_dummies(dat, "Study", keep=c("Poston", "Wagner", "Kerchner"))

# # PlateId was important, but not all plates were top50 in either glmcross-validated or plr estimates
# dat = recode_dummies(dat, "PlateId", keep=c("P0028449", "P0028439", "P0028451", 
                                            # "P0028565", "P0028437", "P0028448",
                                            # "P0028463", "P0028430", "P0028428",
                                            # "P0028447", "P0028452"))

# # Storage did not show up in Rf top50, but "Negative.80C" came up in both GLM CV and PLS
# dat = recode_dummies(dat, "Storage", keep=c("Negative.80C"))

# # best 2of3 (cv, pls, rf; either rf-impur or rf-perm stat can count as only one)
# drops = c("Diagnosis_group", "Visit", "DS1V", "SampleID", "PlasmaPTau181",
          # "Study", "PlateId", "Subarray", "HybControlNormScale", "Fed.Fasted",
          # "Storage", "sv1", "sv5", "sv12", "sv13", "sv19", "sv20", "sv23", 
          # "sv34", "sv36", "sv39", "sv40", "sv42", "sv45", "sv46", "sv49", "sv50",
          # "sv52", "sv53", "sv54", "sv55", "sv56", "sv57", "sv59", "sv60", 
          # "sv61", "sv62")

# all covariates ---------------------------------------------
# drops = c("Diagnosis_group", "Visit", "DS1V", "SampleID", "PlasmaPTau181")

# Run the model ----------------------------------------------
# make covariates
cols = c(names(dat)[2:(grep(first_prot, names(dat)) - 1)])
cols = cols[! cols %in% drops]
print(paste("Running covs:", paste(cols, collapse=", ")))

# make out table
outs = c("PlasmaPTau181") # define dependent variables by column names 
                          # type them out or specify index range, or by regex)
exps = c(names(dat)[grep(first_prot, names(dat))[1]:last(grep(last_prot, names(dat)))]) # same with independent variables
dt1 = CJ(dx=outs, ix=exps)

# fast run function
dt1[, c("beta", "se", "t stat", "pval", "nobs") := get_linreg(dat, dx, ix, cols), by=dt1[,.I]]

# remove Xs needed for formula
dt1[startsWith(ix, "X."), ix := gsub("^X.", ".", ix)]
fwrite(dt1, file=paste(nameset, "linreg_glm", "csv.gz", sep="."))

# cleanup
system("squeue -j $SLURM_JOBID --Format=TimeUsed")


# # model testing ---------------------
# # spot check 10 proteins
# set(123)
# prots = sample(names(dat)[grep(first_prot, names(dat))[1]:last(grep(last_prot, names(dat)))], 10)

# model = glm(paste("PlasmaPTau181", "~", paste(c(prots[2], cols), collapse=" + ")),
             # na.action=na.exclude, data=dat)
# summary(model)
# plot(fitted(model), resid(model, type = "pearson")) 
# abline(0,0, col="red")
# shapiro.test(resid(model, type="pearson"))


