#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 02-linreg_glm.R
##
## Version: 0.0.4
##
## Purpose of script: calculate correlation of plasma ptau with SomaLogic plasma
## proteomic levels
##
## Author: Edward Wilson, Robert B Butler III, Michael P Brown
##
## Date Created: 2023-08-25
##
## Copyright (c) 2023
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##   Linear regression glm for baseline visit only with corrections 
##   for fixed effects: filter covariates based on importance
##   
##   PlateRunDate & ScannerID couldn't be fit to the model, suggesting some 
##   other issue, so dropping those
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
library(tictoc)
tic()
tic()
library(data.table)
library(optparse)
library(tools)
library(stringr)
library(futile.logger)
# library(lmerTest) # loads lme4
setDTthreads(percent=75)

# Nameset helper function
nameset_default <- function(input){
  file_path_sans_ext(basename(input))
}

# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = ".",
    help = "File path to DE input file (*.csv.gz)[default: %default]"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = ".",
    help = "File path to output directory [default: %default]"
  ),
  make_option(
    c("-p", "--prots"),
    type = "character",
    default = "WSS_names.txt",
    help = "File path to prots file [default: %default]"
  ),
  make_option(
    c("-v", "--covs"),
    type = "character",
    default = "",
    help = "Covariates in comma-separated format [default: %default]"
  ),
  make_option(
    c("-s", "--scale-cov"),
    type = "character",
    dest = "scale_cov",
    default = "",
    help = "Which covariates to scale for relmatGlmer in comma-separated format [default: %default]"
  ),
  make_option(
    c("-d", "--recodes"),
    type = "character",
    default = "",
    help = "Recode covariates in comma-separated format and pipe between covs [default: %default]"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project [default: Nameset of input]"
  ),
  make_option(
    c("-b", "--binomial"),
    action = "store_true",
    default = FALSE,
    help = "Use binomial family for logistic regression [default: %default]"
  ),
  make_option(
    c("-f", "--factors"),
    default = "Q1,Q4",
    help = "when binomial == TRUE, what are the levels (ordered) for comparison [default: %default]"
  ),
  make_option(
    c("-c", "--outcome"),
    type = "character",
    help = "String specifying the column in data containing the outcome variable [required]"     
  )
)

args_list <- parse_args(OptionParser(option_list = option_list))

# Nameset default case
if(is.null(args_list$nameset)){
  args_list$nameset <- nameset_default(args_list$input)
}

# check outpath status
if (!dir.exists(args_list$output)) {
  dir.create(args_list$output, recursive = TRUE)
}

# Set up logger
log_fh <- file.path(
  args_list$output, paste(args_list$nameset, "linreg.log", sep = ".")
)
invisible({
  file.create(log_fh)
  flog.layout(layout.format('[~l] ~m'))
  flog.appender(appender.file(log_fh))
  flog.threshold(INFO)
})

# read in fields
fields <- c(
  "input",
  "output",
  "prots",
  "covs",
  "scale_cov",
  "recodes",
  "nameset",
  "binomial",
  "factors",
  "outcome"
)

# check for optparse values
missing_fields <- setdiff(fields, names(args_list))

if (length(missing_fields) > 0) {
  flog.fatal("Missing required arguments: %s", paste(missing_fields, collapse = ", "))
  stop(parse_args(OptionParser(option_list = option_list), args = "--help"),
       call. = FALSE
  )
}

# Log all provided arguments
flog.info("=== Command-line arguments ===")
for (field in fields) {
  assign(field, args_list[[field]])
  flog.info("  %s: %s", field, get(field))
}
flog.info("==============================")


# Split list variables
covs <- strsplit(covs, ",")[[1]]
scale_cov <- strsplit(scale_cov, ",")[[1]]
prots <- readLines(prots)
recodes <- strsplit(unlist(strsplit(recodes, "\\|")), ",")
factors <- strsplit(factors, ",")[[1]]

# Functions ----------------------------------------------------------------
#' Generalized linear model based on a flexible formula
#'
#' @param dat Data.table with the original data
#' @param outcome String specifying the column in dat containing the outcome variable 
#' @param exposure String specifying the column in dat containing the exposure variable
#' @param covs Vector of column names containing covariates to be included in the formula
#' @param logistic booliean indicating logistic regression instead of gaussian (default)
#'
#' @returns A named list: 'Estimate' (beta), 'Std. Error' (se), 't value' (test statistic), 'Pr(>|t|)' (p-value), 'nobs' (number of observations)
get_linreg = function(dat, outcome, exposure, covs=c(), logistic = FALSE){
  model = glm(
    formula = paste(outcome, "~", paste(c(exposure, covs), collapse=" + ")),
    na.action = na.exclude, 
    data = dat,
    family = if (logistic) binomial() else gaussian
  )
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
flog.info("=== Setting up for analysis ===")
dat = fread(input)
names(dat) = make.names(names(dat))
names(prots) <- make.names(prots)
# similar to other name shielding, but reversed for ease of use
orig_outcome <- outcome
outcome <- make.names(outcome)

# add in single category columns with feature importance
if(length(recodes) > 0) {
  new_covs <- c()
  for(recode_item in recodes) {
    if(recode_item[1] %in% covs) {
      dat <- recode_dummies(dat, 
                            dcol = recode_item[1],
                            keep = recode_item[2:length(recode_item)]
      )
      new_covs <- c(new_covs, paste(recode_item[1], recode_item[2:length(recode_item)], sep = "."))
      covs <- covs[covs != recode_item[1]]
    } else {
      flog.info(sprintf("%s not found in input, dropping & skipping recode...", recode_item[1]))
      covs <- covs[covs != recode_item[1]]
    }
  }
  covs <- c(covs, new_covs)
}
# list covariates
flog.info(paste("Running covs:", paste(covs, collapse = ", ")))

# Setting up binomial for glm
if (binomial) {
  # If outcome is character, convert to factor:
  if(is.character(dat[[outcome]]) & length(unique(dat[[outcome]])) == 2) {
    unique_vals <- unique(dat[[outcome]])
    if (
      length(unique_vals) != 2 
      || !all(unique_vals %in% factors) 
      || !all(factors %in% unique_vals)
    ) {
      flog.error("Invalid factors specified for '%s'. Provided: %s, Available: %s", 
                 outcome, 
                 paste(factors, collapse = ", "), 
                 paste(unique_vals, collapse = ", ")
      )
      stop("exactly two factors must be specified for binomial models")
    }
    dat[[outcome]] <- factor(dat[[outcome]], levels = factors)
  }
  
  # scale covariates for logistic regression
  dat[, (scale_cov) := lapply(.SD, function(x) scale(x)[,1]), .SDcols = scale_cov]
  cols <- c("beta", "se", "z stat", "pval", "nobs")
} else {
  cols <- c("beta", "se", "t stat", "pval", "nobs")
}

# make out table
dt1 = CJ(dx = outcome, ix = names(prots))

elapsed <- toc(log = TRUE, quiet = TRUE)
flog.info("=== Analysis started: %s ===", elapsed$callback_msg)
flog.info(
  "Outcome: %s, Proteins: %d, Cores: %d", outcome, nrow(dt1), getDTthreads()
)

# fast run function
dt1[, (cols) := get_linreg(dat, dx, ix, covs, logistic = binomial), by=.I]

# remove Xs needed for formula
dt1$ix <- plyr::mapvalues(x = dt1$ix, from = names(prots), to = prots)
dt1$dx <- plyr::mapvalues(x = dt1$dx, from = outcome, to = orig_outcome)
# dt1[startsWith(ix, "X."), ix := gsub("^X.", ".", ix)]
fwrite(dt1, file=file.path(output, paste(nameset, "linreg_glm", "csv.gz", sep=".")))

# cleanup
end <- toc(log = TRUE, quiet = TRUE)
flog.info("=== Total run time: %s ===", end$callback_msg)
# system("squeue -j $SLURM_JOBID --Format=TimeUsed")


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


