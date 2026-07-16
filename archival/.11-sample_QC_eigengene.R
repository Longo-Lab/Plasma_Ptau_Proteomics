#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 11-sample_QC_eigengene.R
##
## Version: 0.0.1
##
## Purpose of script: preparing eigengenes for quantile analysis
##
## Author: Robert B Butler III, Michael P Brown, Edward N Wilson III
##
## Date Created: 2023-07-04
##
## Copyright (c) 2023
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##   As noted in 01a-Examining_SomaLogic_QC.R, it is really important to use the
##   data that is scale to the external reference. In this case, eigengenes are being 
##   transformed into upper and lower quantiles 
##
##   REMEMBER, WHEN COMMENTING IN CODE, DO NOT REFERENCE INDIVIDUAL INFORMATION:
##   SAMPLES/PIDNS, DATES OR ANY INDIVIDUAL LEVEL INFO
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(data.table)
library(RNOmni)
library(dplyr)     # for data wrangling
library(ggplot2)   # for awesome plotting
library(patchwork)
library(ggrepel)
library(RColorBrewer)
library(matrixStats)
library(sva)
# Modeling packages
library(pls)     # for fitting PLS models
library(caret)    # for tuning models
library(ranger)   # a c++ implementation of random forest 
library(vip)       # for variable importance
library(optparse)
# library(lme4)


# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "dataProt_SS-205063.merged.QC.csv",
    help = "File path to DE input file [default: %default]"
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
    default = "prot_names.txt",
    help = "File path to prots file [default: %default]"
  ),
  make_option(
    c("-v", "--covs"),
    type = "character",
    default = "",
    help = "Terminal input of covariates in comma-separated format [default: %default]"
  ),
  make_option(
    c("-e", "--eigengenes"),
    type = "character",
    default = "module_scores/dataProt_SS-205063.STRING.INT.eigengenes.csv",
    help = "File path to eigengenes csv [default: %default]"
  ),
  make_option(
    c("-m", "--pcs"),
    type = "character",
    default = "genomic/mypcair.PCs.csv",
    help = "File path to PCs csv [default: %default]"
  ),
  make_option(
    c("-c", "--outcome"),
    type = "character",
    help = "String specifying the column in data containing the outcome variable [required]"     
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project [default: Nameset of input]"
  )
)

args_list <- parse_args(OptionParser(option_list = option_list))


# read in fields
fields <- c(
  "input",
  "output",
  "prots",
  "covs",
  "eigengenes",
  "pcs",
  "outcome",
  "nameset"
)

# check for values
for (field in fields) {
  if (!field %in% names(args_list)) {
    cat(str_c('Argument "', field, '" must be provided\n'), fill = TRUE)
    stop(parse_args(OptionParser(option_list = option_list), args = "--help"),
      call. = FALSE
    )
  }

  assign(field, args_list[[field]])
  print(paste0(field, ": ", get(field)))
}

# Nameset default case
if(is.null(args_list$nameset)){
  args_list$nameset <- nameset_default(args_list$input)
}

# check outpath status
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

# variables
covs <- strsplit(covs,",")[[1]]
prots <- readLines(prots)
# first_prot = prots[1]                                             
# last_prot = prots[length(prots)]

qc_cols <- c(
  "PIDN", "SampleID", "ADRC_ID", "Visit", "Visit_Consensus", "ConnectivityZscore",
  "Diagnosis_group", "DS1V"
)


# Functions ------------------
#' Project Surrogate Variables to New Samples
#'
#' This function learns the relationship between proteins and surrogate variables
#' from a training dataset, then applies this relationship to project SVs onto
#' new samples (e.g., held-out repeated visits).
#'
#' @param prot_train A data frame of protein expression data for training samples.
#'   Rows are samples, columns are proteins. Should match the samples used to
#'   estimate the SVs in \code{svobj}.
#' @param prot_all A data frame of protein expression data for all samples
#'   (including training and new samples). Rows are samples, columns are proteins
#'   (same proteins as \code{prot_train}). Sample order does not need to match
#'   \code{prot_train}.
#' @param svobj An object returned by \code{sva::sva()}, containing the
#'   surrogate variables estimated on the training data. Must have a component
#'   \code{sv} which is a matrix of surrogate variables with rownames matching
#'   \code{rownames(prot_train)}.
#'
#' @return A matrix of projected surrogate variables with dimensions
#'   \code{nrow(prot_all)} x \code{ncol(svobj$sv)}. Rows correspond to samples
#'   in \code{prot_all}, columns correspond to surrogate variables. Row names
#'   are sample IDs from \code{rownames(prot_all)}, column names are "sv1",
#'   "sv2", etc.
#'
#' @details
#' The function fits a linear model for each surrogate variable:
#' \code{SV_i ~ protein_1 + protein_2 + ... + protein_p}
#' using the training data, then uses these models to predict SV values for
#' all samples based on their protein expression profiles.
#'
#' Note: This function expects data frames where rows = samples and columns = proteins,
#' which is the opposite of the format required by \code{sva::sva()}. The function
#' handles this internally.
#'
#' Uses \code{lm.fit()} for computational efficiency while maintaining the
#' numerical stability of linear model fitting.
#'
#' @examples
#' \dontrun{
#' # prot and prot_all are data frames: rows = samples, columns = proteins
#' # svobj was created using: sva(t(prot), mod, mod0, n.sv = 6)
#'
#' # Project to all samples
#' all_svs <- project_svs(prot, prot_all, svobj)
#'
#' # Add to phenotype data
#' pheno_all <- cbind(pheno_all, all_svs[rownames(pheno_all), ])
#' }
#'
#' @seealso \code{\link[sva]{sva}}, \code{\link[sva]{fsva}}
#'
#' @export
project_svs <- function(prot_train, prot_all, svobj) {
  # prot_train and prot_all are data frames: rows = samples, columns = proteins
  # svobj$sv is a matrix: rows = samples (matching prot_train), columns = SVs
  
  # Check that proteins match
  if(!all(colnames(prot_train) == colnames(prot_all))) {
    stop("Protein names (colnames) must match between prot_train and prot_all")
  }
  
  # Check that sample names match between prot_train and svobj
  if(!all(rownames(prot_train) == rownames(svobj$sv))) {
    stop("Sample names must match between prot_train and svobj$sv")
  }
  
  # Convert to matrices and add intercept
  X_train <- as.matrix(cbind(Intercept = 1, prot_train))
  X_all <- as.matrix(cbind(Intercept = 1, prot_all))
  SV_train <- svobj$sv
  
  # Use lm.fit() - faster than lm(), handles collinearity like lm()
  fit <- lm.fit(X_train, SV_train)
  

  # Check if coefficients are valid
  if(all(is.na(fit$coefficients))) {
    stop("All coefficients are NA. Check for perfect collinearity in protein data.")
  }
  
  # Replace NA coefficients with 0 (from collinear predictors)
  coef_matrix <- fit$coefficients
  coef_matrix[is.na(coef_matrix)] <- 0

  # Project onto all samples
  sv_all <- X_all %*% coef_matrix
  
  # Set row and column names
  rownames(sv_all) <- rownames(prot_all)
  colnames(sv_all) <- paste0("sv", 1:ncol(svobj$sv))
  
  return(sv_all)
}


# Main -----------------------
# So, the final INTnorm dataset:
a <- fread(input)

# Map in PCs assumes your primary data table exists 
if (pcs != "") {
  print(sprintf("Reading PCs for inclusion from %s", pcs))
  x <- fread(pcs, colClasses = c(PIDN = "character"))
  a <- merge(a, x, by = "PIDN")
  covs <- unique(c(covs, names(a)[names(a) %like% "^PC\\d{1,2}$"])) # flexible for PCs
}

# Load eigengenes
x <- fread(eigengenes)
a <- merge(a, x, by = "SampleID")

# Inverse rank normalization (INT) to z-scale
a[, (prots) := lapply(.SD, RankNorm), .SDcols=prots]

# cols = c(qc_cols, covs)

# Follow-up visits & Dichotomous filters -----------------------
a <- a[DS1V > 0]
# a$Visit <- factor(a$Visit, levels = sort(unique(a$Visit)), ordered = TRUE)

a[, quartile := cut(get(outcome), 
  breaks = quantile(get(outcome), probs = 0:4/4),
  labels = c("Q1", "Q2", "Q3", "Q4"),
  include.lowest = TRUE
)]

# Filter to only Q1 vs Q4
a <- a[quartile %in% c("Q1", "Q4")]
a$quartile <- droplevels(a$quartile)


# generate surrogate variables ---------------------------------------
print("Generating surrogate variables...")
# drop multiple visits per person (for reduced model complexity)
x <- a[order(DS1V)][, .SD[1], by = PIDN]

# matrix of proteins
cols <- c("SampleID", prots)
prot <- data.frame(x[, ..cols], row.names = "SampleID")
prot_all <- data.frame(a[, ..cols], row.names = "SampleID")

# matrix of phenotypes to regress
# too many PIDNs relative to Samples, so including PIDN won't work, and without it sva is computationally singular
# when PIDN > SampleID * 0.5, reintroduce as Visit number fixed effect
cols <- c("SampleID", "quartile", covs)
pheno <- data.frame(x[, ..cols], row.names = "SampleID")
all(rownames(prot) == rownames(pheno))

# defining models with static set of 6 SVs
mod <- model.matrix( # this will include all covs and the outcome
  as.formula(paste("~", paste(names(pheno), collapse=" + "))), 
  data = pheno
) 
mod0 <- model.matrix( # this will drop the outcome
  as.formula(
    paste("~", paste(names(pheno)[names(pheno) != "quartile"], collapse = " + "))
  ), 
  data = pheno
)

svobj <- sva(t(prot), mod, mod0, n.sv = 6)
# rownames(svobj$sv) <- rownames(prot)

# project svs of multiple visits
all_svs <- data.table(
  project_svs(prot, prot_all, svobj), 
  keep.rownames = "SampleID"
)
a <- merge(a, all_svs, by = "SampleID")

# # extract and add to a
# svs <- data.table(svobj$sv)
# setnames(svs, gsub("V", "sv", names(svs)))
# a <- cbind(a, svs)

# add to covs 
covs <- c(covs, names(all_svs)[names(all_svs) != "SampleID"])


# GLM based feature selection -----------------------------------------
print("Testing feature importance for covariates...")
cols <- c("SampleID", "quartile", covs)

# drop missing values for feature importance testing
pheno <- data.frame(a[order(DS1V)][, .SD[1], by = PIDN][, ..cols], row.names = "SampleID")
pheno <- pheno[complete.cases(pheno),] 

# setup for phenos
# # Standardize: mean=0, sd=1
# pheno$Age <- scale(pheno$Age)[,1]
# pheno$Storage_days <- scale(pheno$Storage_days)[,1]

# formula
fmla <- as.formula(quartile ~ .)
# number of features
n_features <- ncol(pheno) - 1

# include all possible main effects
model = glm(fmla, 
  data = pheno, 
  family = binomial 
) 

# simple training
set.seed(123)
cv_model <- train(
  form = fmla, 
  data = pheno, 
  method = "glm",
  family = binomial, 
  trControl = trainControl(method = "cv", number = 5)
)

# summary(cv_model$finalModel)

# Partial Least Squares feature scoring -------------------------------
# perform 10-fold cross validation on a PLS model tuning the 
# number of principal components to use as predictors from 1-30

# # custom model npcs
# npcs = 39
# model_pls = mvr(outcome ~ ., data=b, ncomp=npcs) 

set.seed(123)
cv_model_pls <- train(
  fmla, 
  data = pheno, 
  method = "pls",
  trControl = trainControl(method = "cv", number = 10),
  preProcess = c("zv", "center", "scale"),
  tuneLength = 100
)

# model with lowest RMSE
cv_model_pls$bestTune
## ncomp
##   1

# results for model with lowest RMSE
cv_model_pls$results %>%
  dplyr::filter(ncomp == pull(cv_model_pls$bestTune))
##  ncomp      RMSE Rsquared       MAE     RMSESD RsquaredSD     MAESD
##      1 0.9399164 0.133929 0.7375119 0.09900189 0.07164967 0.0701038
# plot cross-validated RMSE
# ggplot(cv_model_pls)
p1 <- vip::vip(cv_model, geom="point", num_features = 50, bar = FALSE) + 
  ggtitle("GLM - 5-fold Cross Validation") + 
  geom_vline(xintercept = min(25, floor(length(cv_model$coefnames) / 2)))
p2 <- vip::vip(cv_model_pls, geom="point", num_features = 50, bar = FALSE) +
  ggtitle(paste("Partial Least Squares - Best Fit ", 
                cv_model_pls$bestTune$ncomp, "PCs")) + 
  geom_vline(xintercept = min(25, floor(length(cv_model_pls$coefnames) / 2)))
# p5 <- vip::vip(cv_model_pls, num_features = 55, bar = FALSE) +
  # ggtitle(paste("Partial Least Squares - Best Fit ", npcs, "PCs"))


# Random Forest feature scoring ---------------------------------------

# train a default random forest model
rf1 <- ranger(
  fmla, 
  data = pheno,
  mtry = floor(n_features / 3),
  respect.unordered.factors = "order"
)

default_rmse <- sqrt(rf1$prediction.error)

# create hyperparameter grid
hyper_grid <- expand.grid(
  mtry = floor(n_features * c(.05, .15, .25, .333, .4)),
  min.node.size = c(1, 3, 5, 10), 
  replace = c(TRUE, FALSE),                               
  sample.fraction = c(.5, .63, .8),                       
  rmse = NA                                               
)

# execute full cartesian grid search
for(i in seq_len(nrow(hyper_grid))) {
  # fit model for ith hyperparameter combination
  fit.rf <- ranger(
    formula = fmla, 
    data = pheno,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$min.node.size[i],
    replace         = hyper_grid$replace[i],
    sample.fraction = hyper_grid$sample.fraction[i],
    verbose         = FALSE,
    respect.unordered.factors = 'order',
  )
  # export OOB error 
  hyper_grid$rmse[i] <- sqrt(fit.rf$prediction.error)
}

# assess top 10 models
top = hyper_grid %>%
  arrange(rmse) %>%
  mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) %>%
  head(10)
# mtry min.node.size replace sample.fraction      rmse perc_gain
# 22             5   FALSE            0.63 0.9101213 1.2430475

# re-run model with impurity-based variable importance
rf_impurity <- ranger(
  formula = fmla, 
  data = pheno,
  num.trees = 2000,
  mtry = top$mtry[1],
  min.node.size = top$min.node.size[1],
  sample.fraction = top$sample.fraction[1],
  replace = top$replace[1],
  importance = "impurity",
  respect.unordered.factors = "order",
  verbose = FALSE,
)

# re-run model with permutation-based variable importance
rf_permutation <- ranger(
  formula = fmla, 
  data = pheno,
  num.trees = 2000,
  mtry = top$mtry[1],
  min.node.size = top$min.node.size[1],
  sample.fraction = top$sample.fraction[1],
  replace = top$replace[1],
  importance = "permutation",
  respect.unordered.factors = "order",
  verbose = FALSE,
)

p3 <- vip::vip(rf_impurity, geom="point", num_features = 50, bar = FALSE) + 
  ggtitle("Random Forest - Impurity") + 
  geom_vline(xintercept = min(25, floor(n_features/2)))
p4 <- vip::vip(rf_permutation, geom="point", num_features = 50, bar = FALSE) +
  ggtitle("Random Forest - Permutation") + 
  geom_vline(xintercept = min(25, floor(n_features/2)))

# gridExtra::grid.arrange(p1, p2, nrow = 1)

# feature interaction
# layout <- "
# AACCCCCC
# BBCCCCCC
# "
pdf(
  file = file.path(output, paste(nameset, "feature_importance.pdf", sep=".")), 
  width = 11, 
  height = 11
)
p1 + p2 + p3 + p4
dev.off()

# Post-covariate determination ----------------------------------------------
cols = c(qc_cols, covs, outcome, "quartile", prots)
fwrite(
  a[, ..cols], 
  file = file.path(output, paste(nameset, "INTnorm.csv.gz", sep="."))
)



# Testing -------------------------------------------------------------------
# a = data.table(SampleID=seq(6), 
               # PIDN=c(rep("A",2),rep("B",3),rep("C",1)), 
               # Visit=c(seq(2), seq(3), 1), 
               # Date.of.draw=c("3/12/2018", "4/23/2019", "12/28/2016", "2/14/2018",
                      # "7/3/2019", "4/9/2012"))
# a[20:30, .(SampleID, PIDN, Visit, Date.of.draw, Visit_Consensus, DS1V)]
# ---------------------------------------------------------------------------
# p <- a %>% ggplot(aes(outcome, CSRP3.9171.11.3))
# alpha = 0.5
# p + 
  # geom_point(aes(col = PlateId), alpha = alpha) + 
  # geom_smooth(aes(col = PlateId), method = 'glm') +
  # theme(legend.background = element_rect(fill = "transparent"),
        # legend.justification = c(0, 1),
        # legend.position = c(0, 1))

# # # # MARS based feature scoring ---------------------------------------
# # # # examined potential covariates
# # # cols = c(names(a)[2:(grep(first_prot, names(a)) - 1)])
# # # cols = cols[! cols %in% c("Diagnosis_group", "PIDN", "Visit", "DS1V")]
# # # b = a[, ..cols]
# # # b = b[complete.cases(b)] # drop missing values for feature importance testing

# # # # earth quick feature scoring to examine covariates
# # # fit <- earth(outcome ~ ., data=b, degree=2)
# # # evimp(fit, trim=FALSE) # trim=FALSE to show unused variables
# # # p0 = ~plot(evimp(fit), main="Earth (MARS) Model Importance")

# # # # create a tuning grid
# # # hyper_grid <- expand.grid(
  # # # degree = 1:2, 
  # # # nprune = seq(5, 200, length.out = 10) %>% floor()
# # # )

# # # # Cross-validated model
# # # set.seed(123)  # for reproducibility
# # # cv_mars <- train(
  # # # x = subset(b, select = -outcome),
  # # # y = b[, outcome],
  # # # method = "earth",
  # # # metric = "RMSE",
  # # # trControl = trainControl(method = "cv", number = 10),
  # # # tuneGrid = hyper_grid
# # # )

# # # # View results
# # # cv_mars$bestTune
# # # ##    nprune degree
# # # ## 16     56      2

# # # cv_mars$results %>%
  # # # filter(nprune == cv_mars$bestTune$nprune, degree == cv_mars$bestTune$degree)
# # # ##   degree nprune    RMSE  Rsquared      MAE   RMSESD RsquaredSD    MAESD
# # # ## 1      2     56 26817.1 0.8838914 16439.15 11683.73 0.09785945 1678.672

# # # ggplot(cv_mars)

