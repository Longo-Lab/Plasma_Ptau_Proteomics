#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 02-linreg_lme4.R
##
## Version: 0.0.1
##
## Purpose of script: calculate correlation of plasma ptau with SomaLogic plasma
## proteomic levels
##
## Author: Edward Wilson, Robert B Butler III
##
## Date Created: 2023-05-17
##
## Copyright (c) 2023
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##   Linear regression using lme4 with corrections for fixed effects: Age, 
##   Gender, Diagnosis, Study, Storage Days; random effect PIDN
##   
##   This takes a while, so batch unless you like waiting
##   Usage:
##     sbatch -J ptau --mem=12G -c 16 -t 06:00:00 -p interactive \
##       -o %A_02-linreg_lme4.log \
##       --wrap "ml R/4.2.2; Rscript 02-linreg_lme4.R"
##
##   REMEMBER, WHEN COMMENTING IN CODE, DO NOT REFERENCE INDIVIDUAL INFORMATION:
##   SAMPLES/PIDNS, DATES OR ANY INDIVIDUAL LEVEL INFO
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(data.table)
library(lmerTest) # loads lme4
setDTthreads(percent=75)

# directory locations
rootdir = '/labs/flongo/'
wk_dir = paste0(rootdir, '/plasma_ptau_proteomics')
setwd(paste0(rootdir, '/plasma_ptau_proteomics'))
refdir = paste0(rootdir, '/reference/biodomains')

# variables
nameset = "dataProt_SS-205063.pTau_added"
first_prot = "CRYBB2.10000.28.3"
last_prot = "IRF6.9999.1.3"

# Functions -------------------------
get_linreg = function(dat, outcome, exposure){
  model = lmer(get(outcome) ~ get(exposure) + Age + Gender + Study + PlateId +
               Fed.Fasted + Time.To.Freeze + HybControlNormScale + Freeze.Thaw + 
               Time.To.Decant + Storage_days + Time.To.Spin + 
               svd1 + svd2 + svd3 + svd4 + svd5 + 
               (1|PIDN),
               na.action=na.exclude, data=dat)
  # using Satterthwaite's approximation of significance
  a = summary(model)
  return(as.list(c(a$coefficients[2, -3], nobs=nobs(model))))
}

get_linreg2 = function(dat, outcome, exposure){
  model = lmer(get(outcome) ~ get(exposure) + Age + Gender + Study + PlateId +
               Fed.Fasted + Time.To.Freeze + HybControlNormScale + Freeze.Thaw + 
               Time.To.Decant + Storage_days + Time.To.Spin + 
               svd1 + svd2 + svd3 + svd4 + svd5 +
               PIDN + (1|Days_since_first_Visit),
               na.action=na.exclude, data=dat)
  # using Satterthwaite's approximation of significance
  a = summary(model)
  return(as.list(c(a$coefficients[2, -3], nobs=nobs(model))))
}


# Main -------------------------
dat = fread('dataProt_SS-205063.pTau_added.INTnorm.csv')
# dat$Diagnosis_group = factor(dat$Diagnosis_group, c("HC", "SCI", "MCI", "MCI-AD", "AD", "LBD"))

# make out table
outs = c("PlasmaPTau181") # define dependent variables by column names 
                          # type them out or specify index range, or by regex)
exps = c(names(dat)[grep(first_prot, names(dat))[1]:last(grep(last_prot, names(dat)))]) # same with independent variables
dt1 = CJ(dx=outs, ix=exps)

# fast run function
dt1[, c("beta", "se", "t stat", "pval", "nobs") := get_linreg(dat, dx, ix), by=dt1[,.I]]

fwrite(dt1, file=paste(nameset, "linreg_lme4", "csv", sep="."))

# cleanup
system("squeue -j $SLURM_JOBID --Format=TimeUsed")


# model testing ---------------------
model = lmer(PlasmaPTau181 ~ GABARAP.22966.20.3 + Age + Gender + Study + PlateId +
               Fed.Fasted + Time.To.Freeze + HybControlNormScale + Freeze.Thaw + 
               Time.To.Decant + Storage_days + Time.To.Spin + 
               svd1 + svd2 + svd3 + svd4 + svd5 + 
               (1|PIDN),
             na.action=na.exclude, data=dat)
summary(model)
plot(fitted(model), resid(model, type = "pearson")) 
abline(0,0, col="red")

model = lmer(PlasmaPTau181 ~ GABARAP.22966.20.3 + Age + Gender + Study + PlateId +
               Fed.Fasted + Time.To.Freeze + HybControlNormScale + Freeze.Thaw + 
               Time.To.Decant + Storage_days + Time.To.Spin + 
               svd1 + svd2 + svd3 + svd4 + svd5 + 
               PIDN (1|Days_since_first_Visit),
             na.action=na.exclude, data=dat)
summary(model)
plot(fitted(model), resid(model, type = "pearson")) 
abline(0,0, col="red")



# # # test data --------------------------
# # dat = fread('https://articledatas3.s3.eu-central-1.amazonaws.com/CanonicalCorrelationAnalysisData.csv')

# # get_linreg = function(dat, outcome, exposure){
  # # model = lmer(reformulate(c(exposure, "(1|IQ)", "HrsTrn", "HrsWrk"), 
                           # # response=outcome),
               # # na.action = na.exclude, data=dat)
  # # # using Satterthwaite's approximation of significance
  # # a = summary(model)
  # # return(as.list(c(a$coefficients[2, -3], nobs=nobs(model))))
# # }

# # # make out table
# # outs = c(names(dat)[7:9]) # dependent vars
# # exps = c(names(dat)[1:3]) # independent vars
# # dt1 = CJ(dx=outs, ix=exps)

# # # fast run function
# # bench = function(x) {
# # setDTthreads(x)
# # dt1[, c("beta", "se", "t stat", "pval", "nobs") := get_linreg(dat, dx, ix), by=dt1[,.I]]
# # }

# # library(microbenchmark)
# # microbenchmark(bench(1), bench(2), bench(3), bench(4), bench(8), bench(16), times=3L)
# # # test
# # model = lmer(ClientSat ~ PsychTest1 + (1|IQ) + HrsTrn + HrsWrk, na.action=na.exclude, data=dat)
# # get_linreg(dat, "ClientSat", "PsychTest1")
