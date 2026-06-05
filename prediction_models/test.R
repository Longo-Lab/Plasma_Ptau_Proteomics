#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 09-Prediction_model_construction.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Extract STRING network modules and prepar for Scoring with weighted GSEA
#'  
#'
#' Author: Audrey Weber and Robert B Butler III
#'
#' Date Created: 2025-11-05
#'
#' Copyright (c) 2025
#' Email: audweber@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'   v0.0.1 - initial build
#'   
#'   
#'   
#'
#' ---------------------------
#' ## load up the packages we will need:  (uncomment as required)
library(data.table)
library(RNOmni)
library(dplyr)     # for data wrangling
library(ggplot2)   # for awesome plotting
library(patchwork)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(matrixStats)
library(sva)
library(pls)     # for fitting PLS models
library(caret)    # for tuning models
library(ranger)   # a c++ implementation of random forest 
library(vip)       # for variable importance
library(optparse)
library(auditor)
library(joinet)
library(randomForestSRC)
library(DALEX)
library(glmnet)


#' # read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Biomarker results to use for probe filtering, keeps most significant"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project [default: Nameset of input]"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "cell_type_prot",
    help = "File path to output directory [default: %default]"
  ),
  make_option(
    c("-t", "--threads"),
    default = 0,
    help = "num threads [default: parallelly::availableCores()]"
  )
)

args_list <- parse_args(OptionParser(option_list = option_list))

# Nameset default case
if(is.null(args_list$nameset)){
  args_list$nameset <- nameset_default(args_list$input)
}

# read in fields
fields <- c(
  "input",
  "nameset",
  "output",
  "threads"
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

# check the current active plan
if (threads == 0) {
  threads <- as.integer(parallelly::availableCores()) - 2
} else if (is.numeric(threads) & threads != 0) {
  threads <- as.integer(threads)
} else {
  cat("threads must be an integer value")
  stop(parse_args(OptionParser(option_list = option_list), args = "--help"),
       call. = FALSE
  )
}
print(sprintf("Setting up %d threads for samples", threads))

# input variables
input <- strsplit(input, ",")[[1]]

# check outpath status
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

# input processing (if using optparse)
covs <- strsplit(covs,",")[[1]]

# Functions ------------------
#' Recodes dummy categorical variables and can subset by vector of categories
#'
#' @param dat Data.table with the original data
#' @param dcol String representing the column to be split
#' @param keep Vector of strings that are the categories within dcol to extract. default is all c()
#' @returns modified data.table with new dummy variable columns named "dcol.category" after original column
#' 
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

rroc_curve <- function(y, yhat, n_points = 400) {
  stopifnot(length(y) == length(yhat))
  ok <- is.finite(y) & is.finite(yhat)
  y <- y[ok]; yhat <- yhat[ok]
  e <- yhat - y
  s_min <- -max(e); s_max <- -min(e)
  shifts <- seq(s_min - 1e-8, s_max + 1e-8, length.out = max(20, n_points))
  UNDER <- numeric(length(shifts)); OVER <- numeric(length(shifts))
  for (i in seq_along(shifts)) {
    ee <- e + shifts[i]
    OVER[i]  <- sum(ee[ee > 0])     # x-axis
    UNDER[i] <- sum(-ee[ee < 0])    # y-axis
  }
  data.frame(OVER = OVER, UNDER = UNDER, shift = shifts)
}

rroc_aoc <- function(rroc_df) {
  d <- rroc_df[order(rroc_df$OVER), ]
  sum((d$UNDER[-1] + d$UNDER[-nrow(d)]) * diff(d$OVER) / 2)
}

# Drops values that are not finite and calculates SSE and TSS 
# so that we can evaluate RSE
rse_safe <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]
  sse <- sum((obs - pred)^2)
  tss <- sum((obs - mean(obs))^2)
  if (tss == 0) return(NA_real_)
  sse / tss
}


# Main ------------------------------------------------------------------------------------------------------------
#Load data as datatable 
a <- fread(input)

# Map in PCs assumes your primary data table exists 
if (pcs != "") {
  print(sprintf("Reading PCs for inclusion from %s", pcs))
  x <- fread(pcs, colClasses = c(PIDN = "character"))
  a <- merge(a, x, by = "PIDN")
  covs <- unique(c(covs, names(a)[names(a) %like% "^PC\\d{1,2}$"])) # flexible for PCs
}

# Add cell scores
if (prediction_scores != "") {
  print(sprintf("Reading Prediction scores for inclusion from %s", prediction_scores))
  x <- fread(prediction_scores, colClasses = c(SampleID = "character"))
  a <- merge(a, x, by = "SampleID")
}

# filter outliars 
drop = a[abs(ConnectivityZscore) > 10, SampleID]
a = a[!SampleID %in% drop]

#Only keep baseline visits where plasma is available 
a = a[order(DS1V)][, head(.SD, 1), by = "PIDN"]

# checking for missed duplicates or bad visit labeling
a[duplicated(PIDN), ..qc_cols]

# All variables need to be numeric & dummy coded
a[, Gender := as.numeric(Gender == "M")]
recode_dummies(a, "APOE", keep = c("e34", "e44"))
covs <- unique(c(covs[covs != "APOE"], "APOE.e34", "APOE.e44"))

# Handling incomplete cases
cols <- c("SampleID", outcome, covs, feature)

names(a) <- names(name_map <- setNames(names(a), make.names(names(a), unique = TRUE)))
# Setting up data for input 
pheno <- a[complete.cases(a[, ..cols]), ..cols] # at this point everything is in pheno and is ordered, just keep calling this
pheno <- `rownames<-`(setDF(copy(a[complete.cases(a[, ..cols]), 
                                   setdiff(cols, "SampleID"), with = FALSE])), 
                      make.unique(as.character(a[complete.cases(a[, ..cols]), SampleID])))

#Save a copy of pheno for later 
pheno2 <- pheno

# set predictors of interest 
variants <- list(
  ALL = character(0),
  NO_hybridNormControlScale = "HybControlNormScale",
  NO_PCSandAPOE = c(
    "APOE.e34", "APOE.e44", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"
  ),
  NO_PCS = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8")
)

# Join all predictors from pheno2 and get rid of outcomes
predictor_pool <- setdiff(names(pheno), outcome)

# Tuning grid that tests different alphas and lambdas for glmnet
tuning_grid <- expand.grid(
  alpha  = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),  
  lambda = seq(0.001, 0.1, length = 20)
)

#need a different tuning grid for joinet that accounts for base and meta! 
alpha_grid <- seq(0, 1, by = 0.1)
grid <- CJ(alpha.base = alpha_grid, alpha.meta = alpha_grid)
k_inner = 5 

rows <- list()


# --- Start looping over variants and models ---
for (tag in names(variants)) {
  message("\n==== Variant: ", tag, " ====")
  drops <- variants[[tag]]
  x_cols <- setdiff(names(pheno2), unique(c(outcome, drops)))
  
  
  df <- pheno[, c(outcome, x_cols), drop = FALSE]
  df <- df[complete.cases(df), ]
  if (nrow(df) < 10) next
  
  
# --- glmnet ---
  
  set.seed(42)
  glmnet_model <- train(
    as.formula(paste(outcome, "~ .")),
    data = df,
    method = "glmnet",
    trControl = trainControl(method = "cv", number = K, savePredictions = "final"),
    tuneGrid = tuning_grid,
    metric = "RMSE"
  )
  
  pred_glm <- glmnet_model$pred$pred
  obs_glm <- glmnet_model$pred$obs
  rmse_glm <- sqrt(mean((pred_glm - obs_glm)^2))
  mae_glm <- mean(abs(pred_glm - obs_glm))
  r2_glm <- cor(pred_glm, obs_glm)^2
  rse_glm <- rse_safe(obs_glm, pred_glm)
  mape_glm <- mean(abs((pred_glm - obs_glm)/obs_glm), na.rm=TRUE)*100
  aoc_glm <- rroc_aoc(rroc_curve(obs_glm, pred_glm))
  
  
  rows[[length(rows)+1]] <- data.table(
    model = "glmnet", variant = tag,
    rmse = rmse_glm, mae = mae_glm, r2 = r2_glm,
    rse = rse_glm, mape = mape_glm, aoc_rroc = aoc_glm
  )
  
  
  # --- joinet ---
  Y <- as.matrix(df[, outcome, drop=FALSE])
  X <- as.matrix(df[, x_cols, drop=FALSE])
  foldid.ext <- sample(rep(1:K, length.out = nrow(df)))
  foldid.int <- sample(rep(1:K_inner, length.out = nrow(df)))
  
scores <- vector("list", length(tuning_grid)^2)
g <- 1

for (ab in alpha_grid) for (am in alpha_grid) {
    mcv <- cv.joinet(Y, X, family = "gaussian", type.measure = "mse",
                     alpha.base = ab, alpha.meta = am,
                     foldid.ext = foldid.ext, foldid.int = foldid.int,
                     nfolds.ext = K, nfolds.int = K_inner)
    scores[[g]] <- list(alpha.base = ab, alpha.meta = am,
                        meta_mse = mean(mcv["meta",]))
    g <- g + 1
    }
  
df_scores <- rbindlist(scores)
best <- df_scores[which.min(meta_mse)]
  
# True OOF
  oof_pred <- rep(NA_real_, nrow(Y))
  for (f in sort(unique(foldid.ext))) {
    te <- which(foldid.ext == f)
    tr <- setdiff(seq_len(nrow(Y)), te)
    fit <- joinet(Y[tr,,drop=FALSE], X[tr,,drop=FALSE],
                  alpha.base = best$alpha.base, alpha.meta = best$alpha.meta,
                  family = "gaussian", nfolds = K_inner)
    pred <- predict(fit, newx = X[te,,drop=FALSE])
    oof_pred[te] <- as.numeric(pred$meta)
  }
  
  rmse_j <- sqrt(mean((Y - oof_pred)^2))
  mae_j <- mean(abs(Y - oof_pred))
  r2_j <- cor(Y, oof_pred)^2
  rse_j <- rse_safe(Y, oof_pred)
  mape_j <- mean(abs((Y - oof_pred)/Y), na.rm=TRUE)*100
  aoc_j <- rroc_aoc(rroc_curve(Y, oof_pred))
  
  rows[[length(rows)+1]] <- data.table(
    model = "joinet", variant = tag,
    rmse = rmse_j, mae = mae_j, r2 = r2_j,
    rse = rse_j, mape = mape_j, aoc_rroc = aoc_j
  )
}


# Combine and rank
results <- rbindlist(rows)
results <- results[order(rmse)]
fwrite(results, file.path(output, "model_variant_summary.csv"))


print(results)