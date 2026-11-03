#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 09_all_prediction_models.R
##
## 
##
## Purpose of script: test individual glmnet, multivariate glmnet, joinet, individual random forest, multivariate randon forest 
##for prediction performance of plasma biomarkers from plasma proteoic readouts 
##
## Author: Audrey Weber, Robert B Butler III
##
## Date Created: 2026-01-20
##
## 
## Email: audweber@stanford.edu
##
## ---------------------------
##
## Notes:
##   
##
##   REMEMBER, WHEN COMMENTING IN CODE, DO NOT REFERENCE INDIVIDUAL INFORMATION:
##   SAMPLES/PIDNS, DATES OR ANY INDIVIDUAL LEVEL INFO
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(data.table)
library(RNOmni)
library(Hmisc)     # rcorr.cens()/somers2() for C-index -- must load before
                    # the joinet section's C-index block runs
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
library(parallel)
library(doParallel)
library(data.table)



# Inputs ----------------
inpath <- "/mnt/box/Robert Butler's Files/plasma_ptau_proteomics" # don't repeat things you don't have to
# inpath <- "C:/Users/audweber/Box/plasma_ptau_proteomics" # don't repeat things you don't have to
setwd(inpath)
input <- file.path(inpath, "dataProt_SS-205063.merged.QC.csv")
prots <- file.path(inpath, "prot_names.txt") 
pcs <- file.path(inpath, "genomic", "mypcair.PCs.csv")
base_WSS <- file.path(inpath, "module_scores", 
  "dataProt_SS-205063.STRING.INT.baseline.WSS.csv"
)
follow_WSS <- file.path(inpath, "module_scores", 
  "dataProt_SS-205063.STRING.INT.follow-up.WSS.csv"
)
output <- 'prediction_models/STRING_WSS_v2/model_info/joinet'
covs <- 'Age,Gender,APOE,HybControlNormScale'
# output <- 'prediction_models/eigengenes/mrf/INT'
# covs <- 'Age,Gender,APOE,HybControlNormScale'
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
outcome <- make.names(names(fread(base_WSS, nrows = 0))[-1])
feature <- c(
  "PlasmaPTau181", "PlasmaAB142P","PlasmaAB140P","PlasmaABRatio", 
  "PlasmapTau217","PlasmapTau217_AB42Ratio","PlasmaGFAP","PlasmaNfL"
)
qc_cols <- c(
  "PIDN", "SampleID", "ADRC_ID", "Visit", "Visit_Consensus", "ConnectivityZscore",
  "Diagnosis_group", "Date.of.draw", "DS1V"
)


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

#' Out-of-sample projection matching RNOmni::RankNorm (k = 0.375,
#' ties = "average") exactly for in-sample points; for a new point, treats
#' it as a single insertion into the frozen training reference. Used to
#' project follow-up biomarker values onto the baseline-fit INT reference
#' without leaking follow-up information into the baseline normalization.
project_rank_norm <- function(train_vals, new_vals, k = 0.375) {
  n <- length(train_vals)
  sorted_train <- sort(train_vals)
  r_new <- findInterval(new_vals, sorted_train) + 1
  denom <- (n + 1) - 2 * k + 1
  qnorm((r_new - k) / denom)
}

get_xcols <- function(tag) {
  cat(sprintf("Selecting predictors for variant: %s\n", tag))
  cat(sprintf("Total predictors in pool: %d\n", length(predictor_pool)))
  
  if (tag == "ALL") {
    # use everything that isn't an outcome
    selected <- predictor_pool
  } else {
    # Remove the predictors listed in variants[[tag]]
    selected <- setdiff(predictor_pool, variants[[tag]])
  }
  
  cat(sprintf("Selected %d predictors: %s\n", length(selected), paste(selected, collapse=", ")))
  
  if (length(selected) == 0) {
    stop(sprintf("[%s] No predictors left after removing specified ones. Check variant definition.", tag))
  }
  
  return(selected)
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

# Add baseline cohort WSS scores (now from baseline set so no filtering needed)
if (base_WSS != "") {
  print(sprintf("Reading Prediction scores for inclusion from %s", base_WSS))
  x <- fread(base_WSS, colClasses = c(SampleID = "character"))
  a <- merge(a, x, by = "SampleID")
}

# # filter outliars 
# drop = a[abs(ConnectivityZscore) > 10, SampleID]
# a <- a[!SampleID %in% drop][, !"Neddylation", with = FALSE]

# #Only keep baseline visits where plasma is avaliable 
# a = a[order(DS1V)][, head(.SD, 1), by = "PIDN"]

# checking for missed duplicates or bad visit labeling
a[duplicated(PIDN), ..qc_cols]

names(a) <- names(name_map <- setNames(names(a), make.names(names(a), unique = TRUE)))

# All variables need to be numeric & dummy coded
a[, Gender := as.numeric(Gender == "M")]
recode_dummies(a, "APOE", keep = c("e34", "e44"))
covs <- unique(c(covs[covs != "APOE"], "APOE.e34", "APOE.e44"))


# Handling incomplete cases
cols <- c("SampleID", outcome, covs, feature)

# Setting up data for input 
pheno <- a[complete.cases(a[, ..cols]), ..cols] # at this point everything is in pheno and is ordered, just keep calling this

# INT-normalize biomarkers in-sample on baseline (training); the raw
# baseline values are kept as the frozen reference for projecting
# follow-up biomarkers below, so no hold-out information leaks in.
baseline_feature_ref <- copy(pheno[, ..feature])
pheno[, (feature) := lapply(.SD, RankNorm), .SDcols = feature]

pheno <- `rownames<-`(setDF(copy(pheno[, setdiff(cols, "SampleID"), with = FALSE])),
                      make.unique(as.character(pheno$SampleID)))

#Save a copy of pheno for later
pheno2 <- pheno


variants <- list(
  Biomarkers_Only = c("HybControlNormScale","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8", "APOE.e34", "APOE.e44", "Age", "Gender"),
  Biomarkers_Age_Sex = c("HybControlNormScale","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8", "APOE.e34", "APOE.e44"),
  Age_Sex_HybridNormControl = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8", "APOE.e34", "APOE.e44" ),
  Age_Sex_HybridNormControl_PCs = c("APOE.e34", "APOE.e44"),
  ALL = character(0)
)

predictor_pool <- setdiff(names(pheno2), outcome)

#CV settings 
K_outer <- 5
K_inner <- 10
alpha_grid <- seq(0, 1, 0.2)  # Coarser grid for faster computation

# Set up parallel backend
cores <- detectCores() - 1  # Use all but one core
cl <- makeCluster(cores)
registerDoParallel(cl)
####################################################JOINET##############################################
mod <- "joinet"
#########################################################################################################
# Export necessary functions and data to the parallel workers
clusterExport(cl, c("cv.joinet", "joinet", "get_xcols", "pheno2", "outcome", "K_outer", "K_inner"))

final_models <- list()

for (tag in names(variants)) {
  cat("\n==== Variant:", tag, "====\n")
  
  # 1) Build data for this variant
  x_cols <- get_xcols(tag)
  
  cat(sprintf("[%s] Number of selected predictors: %d\n", tag, length(x_cols)))
  
  df <- as.data.frame(pheno2)[, c(outcome, x_cols), drop = FALSE]
  cat(sprintf("[%s] Dimensions of initial dataframe: %d rows, %d columns\n", tag, nrow(df), ncol(df)))
  
  df <- df[complete.cases(df), ]
  cat(sprintf("[%s] Dimensions after removing NA: %d rows, %d columns\n", tag, nrow(df), ncol(df)))
  
  if (nrow(df) < max(K_outer, 10)) {
    cat(sprintf("[%s] Too few rows after NA filtering (%d). Skipping.\n", tag, nrow(df)))
    next
  }
  
  Y <- as.matrix(df[, outcome, drop = FALSE])
  X <- as.matrix(df[, x_cols, drop = FALSE])
  
  cat(sprintf("[%s] Data prepared. Rows: %d, Predictors: %d\n", tag, nrow(X), ncol(X)))
  
  # 2) Find best alphas in parallel
  set.seed(42)
  grid <- CJ(alpha.base = alpha_grid, alpha.meta = alpha_grid)
  
  results <- tryCatch({
    foreach(i = 1:nrow(grid), .combine = rbind, .packages = c("data.table"), .errorhandling = 'pass') %dopar% {
      ab <- grid$alpha.base[i]
      am <- grid$alpha.meta[i]
      mcv <- cv.joinet(
        Y = Y, X = X, family = "gaussian",
        type.measure = "mse",
        alpha.base = ab, alpha.meta = am,
        nfolds.ext = K_outer, nfolds.int = K_inner
      )
      mse <- mean(mcv["meta", ])
      data.table(alpha.base = ab, alpha.meta = am, mse = mse)
    }
  }, error = function(e) {
    cat(sprintf("[%s] Error in parallel execution: %s\n", tag, e$message))
    return(NULL)
  })
  
  if (is.null(results) || nrow(results) == 0) {
    cat(sprintf("[%s] No valid results. Skipping this variant.\n", tag))
    next
  }
  
  # Check for errors in results
  error_rows <- which(sapply(results, inherits, "error"))
  if (length(error_rows) > 0) {
    cat(sprintf("[%s] Errors occurred in %d iterations:\n", tag, length(error_rows)))
    for (row in error_rows) {
      cat(sprintf("  Row %d: %s\n", row, results[[row]]$message))
    }
    results <- results[-error_rows, ]
  }
  
  if (nrow(results) == 0) {
    cat(sprintf("[%s] No valid results after error removal. Skipping this variant.\n", tag))
    next
  }
  
  best_result <- results[which.min(mse)]
  best_ab <- best_result$alpha.base
  best_am <- best_result$alpha.meta
  
  cat(sprintf("[%s] Best alphas: alpha.base=%.2f  alpha.meta=%.2f\n", tag, best_ab, best_am))
  
  # 3) Create and save final model for this variant
  model_id <- paste(paste(outcome, collapse="_"), tag, sep = "_")
  
  final_models[[model_id]] <- list(
    model = joinet(Y = Y, X = X, family = "gaussian",
                   alpha.base = best_ab, alpha.meta = best_am),
    outcomes = outcome,
    variant = tag,
    alpha_base = best_ab,
    alpha_meta = best_am,
    x_cols = x_cols
  )
  
  cat(sprintf("Processed model: %s\n", model_id))
}

# Stop cluster
stopCluster(cl)

# # Save all final models
saveRDS(final_models, file.path(output, "final_models_IN.rds"))
# 
# # Create a summary table
summary_table <- data.table(
  Model_ID = names(final_models),
  Variant = sapply(final_models, function(m) m$variant),
  Alpha_Base = sapply(final_models, function(m) m$alpha_base),
  Alpha_Meta = sapply(final_models, function(m) m$alpha_meta),
  N_Features = sapply(final_models, function(m) length(m$x_cols))
)
# 
# Save summary table
fwrite(summary_table, file.path(output, "final_models_summary.csv"))

print(summary_table)


## setting up test data--------------------------------------------------------------------

a <- fread(input)

# Map in PCs assumes your primary data table exists 
if (pcs != "") {
  print(sprintf("Reading PCs for inclusion from %s", pcs))
  x <- fread(pcs, colClasses = c(PIDN = "character"))
  a <- merge(a, x, by = "PIDN")
  covs <- unique(c(covs, names(a)[names(a) %like% "^PC\\d{1,2}$"])) # flexible for PCs
}

# Add module scores
if (follow_WSS != "") {
  print(sprintf("Reading scores for inclusion from %s", follow_WSS))
  x <- fread(follow_WSS, colClasses = c(SampleID = "character"))
  a <- merge(a, x, by = "SampleID")
}

# recode gender/APOE for models
a[, Gender := as.numeric(Gender == "M")]
recode_dummies(a, "APOE", keep = c("e34", "e44"))
covs <- unique(c(covs[covs != "APOE"], "APOE.e34", "APOE.e44"))

# Handling incomplete cases
columns_to_check <- c(  
  "Age","Gender","HybControlNormScale", 
  "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","APOE.e34","APOE.e44",
  "PlasmaPTau181","PlasmaAB142P","PlasmaAB140P","PlasmaABRatio","PlasmapTau217",
  "PlasmapTau217_AB42Ratio", "PlasmaGFAP","PlasmaNfL"
)
drops <- a[!complete.cases(a[, ..columns_to_check ]), SampleID]
test_data <- a[!SampleID %in% drops]
names(test_data) <- names(name_map <- setNames(names(test_data), make.names(names(test_data), unique = TRUE)))

# Cleanup test data objects
cols <- c("SampleID", outcome, covs, feature)

# Project follow-up biomarkers onto the frozen baseline INT reference
# (fit above) rather than re-fitting RankNorm in-sample on follow-up,
# so held-out values never influence their own normalization.
for (f in feature) {
  test_data[[f]] <- project_rank_norm(baseline_feature_ref[[f]], test_data[[f]])
}

test_data <- data.frame(test_data[, ..cols], row.names = "SampleID")

actual_values <- test_data[, outcome, drop = FALSE]

## at this point go to INT values code 

##PREDICTIONS--------------------------------------------------------------------------------------------------

predictions <- list()

# Loop through models and make predictions
for (model_id in names(final_models)) {
  model_info <- final_models[[model_id]]
  
  cat("\nProcessing model:", model_id, "\n")
  
  # Check if x_cols exist in test_data
  if (!all(model_info$x_cols %in% colnames(test_data))) {
    missing_cols <- setdiff(model_info$x_cols, colnames(test_data))
    cat("Missing columns in test data:", paste(missing_cols, collapse = ", "), "\n")
    next
  }
  
  # Prepare test data for this model
  X_test <- test_data[, model_info$x_cols, drop = FALSE]
  
  # Make predictions
  tryCatch({
    pred <- predict(model_info$model, newx = X_test)
    
    # For joinet, we want the 'meta' predictions
    if (is.list(pred) && !is.null(pred$meta)) {
      pred_matrix <- as.matrix(pred$meta)
    } else {
      stop("Unexpected prediction structure")
    }
    
    # Ensure the prediction matrix has the correct column names
    if (is.null(colnames(pred_matrix))) {
      colnames(pred_matrix) <- model_info$outcomes
    }
    
    # Store predictions
    predictions[[model_id]] <- pred_matrix
    
    cat("Processed predictions for model:", model_id, "\n")
    cat("Prediction dimensions:", paste(dim(pred_matrix), collapse = " x "), "\n")
    cat("Predicted outcomes:", paste(colnames(pred_matrix), collapse = ", "), "\n")
  }, error = function(e) {
    cat("Error in prediction for model", model_id, ":", e$message, "\n")
  })
}

Biomarkers_Only <- as.data.frame(predictions[[1]])
Biomarkers_Only <- cbind(SampleID = rownames(actual_values), Biomarkers_Only)

fwrite(Biomarkers_Only, file.path(output, "Biomarkers_Only.csv"))

Biomarkers_Age_Sex <- as.data.frame(predictions[[2]])
Biomarkers_Age_Sex <- cbind(SampleID = rownames(actual_values), Biomarkers_Age_Sex)

fwrite(Biomarkers_Age_Sex, file.path(output, "Biomarkers_Age_Sex.csv"))

Age_Sex_HybridNormControl <- as.data.frame(predictions[[3]])
Age_Sex_HybridNormControl<- cbind(SampleID = rownames(actual_values), Age_Sex_HybridNormControl)

fwrite(Age_Sex_HybridNormControl, file.path(output, "Age_Sex_HybridNormControl.csv"))

Age_Sex_HybridNormControl_PCs <- as.data.frame(predictions[[4]])
Age_Sex_HybridNormControl_PCs <- cbind(SampleID = rownames(actual_values), Age_Sex_HybridNormControl_PCs)

fwrite(Age_Sex_HybridNormControl_PCs, file.path(output, "Age_Sex_HybridNormControl_PCs.csv"))

ALL <- as.data.frame(predictions[[5]])
ALL<- cbind(SampleID = rownames(actual_values), ALL)

fwrite(ALL, file.path(output, "ALL.csv"))

metrics <- lapply(names(predictions), function(model_id) {
  predicted <- predictions[[model_id]]
  model_info <- final_models[[model_id]]
  
  outcome_metrics <- lapply(seq_along(model_info$outcomes), function(i) {
    outcome <- model_info$outcomes[i]
    pred <- predicted[, i]
    actual <- actual_values[[outcome]]
    
    # Remove NA pairs
    valid_idx <- !is.na(pred) & !is.na(actual)
    pred_valid <- pred[valid_idx]
    actual_valid <- actual[valid_idx]
    
    # Basic metrics
    RMSE <- sqrt(mean((pred_valid - actual_valid)^2))
    MAE <- mean(abs(pred_valid - actual_valid))
    #R2 <- cor(actual_valid, pred_valid)^2
    meta_model <- model_info$model$meta[[i]]
    idx <- which(meta_model$lambda == meta_model$lambda.min)
    R2_train <- meta_model$glmnet.fit$dev.ratio[idx]
    
    # C-index (like AUC for continuous outcomes)
    c_index_result <- tryCatch({
      rcorr.cens(pred_valid, actual_valid)
    }, error = function(e) {
      c("C Index" = NA, Dxy = NA, S.D. = NA, n = NA)
    })
    
    c_index <- c_index_result["C Index"]
    c_index_se <- c_index_result["S.D."]
    
    # Alternative: Somers' D
    somers_result <- tryCatch({
      somers2(pred_valid, actual_valid)
    }, error = function(e) {
      c(C = NA, Dxy = NA, S.D. = NA, n = NA)
    })
    
    c_index_somers <- somers_result["C"]
    dxy <- somers_result["Dxy"]
    
    # Bias
    Bias <- mean(pred_valid - actual_valid)
    
    # Sample size
    N <- length(pred_valid)
    
    list(
      outcome = outcome,
      N = N,
      RMSE = RMSE,
      MAE = MAE,
      R2 = R2_train,
      C_index = c_index,
      C_index_SE = c_index_se
    )
  })
  
  combined_metrics <- do.call(rbind, lapply(outcome_metrics, as.data.frame))
  combined_metrics$model_id <- model_id
  combined_metrics$variant <- model_info$variant
  
  return(combined_metrics)
})

metrics_dt <- rbindlist(metrics)

# Reorder columns
metrics_dt <- metrics_dt[, .(model_id, outcome, variant, N, 
                             RMSE, MAE, R2, C_index, C_index_SE)]

print(metrics_dt)

# Format for publication with confidence intervals
metrics_dt_pub <- metrics_dt[, .(
  model_id, outcome, variant, N,
  RMSE = round(RMSE, 3),
  MAE = round(MAE, 3),
  R2 = round(R2, 3),
  C_index = sprintf("%.3f ± %.3f", C_index, C_index_SE)
)]

print(metrics_dt_pub)

write.csv(metrics_dt_pub, file.path(output, "joinet_model_metrics_c_index_STRING_WSS_INT_test.csv"), row.names = FALSE)

#run correlations-------------------------------------------------------------------------------------------------- 
correlations <- data.frame()

for (model_id in names(predictions)) {
  pred_matrix <- predictions[[model_id]]
  model_info <- final_models[[model_id]]
  
  for (i in seq_along(model_info$outcomes)) {
    outcome_i <- model_info$outcomes[i]  # local name -- avoid clobbering the global `outcome` vector
    pred <- pred_matrix[, i]
    actual <- actual_values[[outcome_i]]

    cor_test <- cor.test(pred, actual, method = "pearson", use = "complete.obs")

    correlations <- rbind(correlations, data.frame(
      model_id = model_id,
      outcome = outcome_i,
      variant = model_info$variant,
      Pearson_Correlation = cor_test$estimate,
      Correlation_P_Value = cor_test$p.value,
      CI_Lower = cor_test$conf.int[1],
      CI_Upper = cor_test$conf.int[2],
      Sample_Size = length(pred)
    ))
  }
}

write.csv(correlations, file.path(output, "joinet_model_correlations_STRING_WSS_INT.csv"), row.names = FALSE)

#macros-------------------------------------------------------------------------------------------------------------

macro_metrics_dt <- metrics_dt[, .(
  n_outcomes = .N,
  
  RMSE_mean = mean(RMSE, na.rm = TRUE),
  RMSE_median = median(RMSE, na.rm = TRUE),
  
  MAE_mean = mean(MAE, na.rm = TRUE),
  MAE_median = median(MAE, na.rm = TRUE),
  
  R2_mean = mean(R2, na.rm = TRUE),
  R2_median = median(R2, na.rm = TRUE),
  
  C_index_mean = mean(C_index, na.rm = TRUE),
  C_index_median = median(C_index, na.rm = TRUE)
  
  
), by = .(model_id, variant)]

# save
write.csv(macro_metrics_dt, file.path(output, "joinet_macro_metrics_INT.csv"), row.names = FALSE)

#macros correlations 

cor_dt <- as.data.table(correlations)

# Fix sample size: your current code uses length(pred) (includes NAs).
# Use Sample_Size if it's already correct, otherwise recompute it (see note below).
# Here, we'll assume Sample_Size is OK. If not, use complete-case N.

cor_dt[, N_eff := Sample_Size]
cor_dt[, r := as.numeric(Pearson_Correlation)]

# Fisher z pooling (weights ~ N-3)
cor_dt[, z := atanh(pmax(pmin(r, 0.999999), -0.999999))]  # avoid Inf

macro_cor_dt <- cor_dt[
  !is.na(r) & !is.na(N_eff) & N_eff >= 4,
  .(
    n_outcomes = .N,
    N_total = sum(N_eff, na.rm = TRUE),
    
    r_mean = mean(r, na.rm = TRUE),
    r_median = median(r, na.rm = TRUE),
    
    # Fisher z pooled r
    r_fisher = {
      w <- pmax(N_eff - 3, 0)
      zbar <- sum(w * z, na.rm = TRUE) / sum(w, na.rm = TRUE)
      tanh(zbar)
    },
    
    # Approx SE/CI for pooled Fisher z
    r_fisher_ci_lower = {
      w <- pmax(N_eff - 3, 0)
      zbar <- sum(w * z, na.rm = TRUE) / sum(w, na.rm = TRUE)
      se <- sqrt(1 / sum(w, na.rm = TRUE))
      tanh(zbar - 1.96 * se)
    },
    r_fisher_ci_upper = {
      w <- pmax(N_eff - 3, 0)
      zbar <- sum(w * z, na.rm = TRUE) / sum(w, na.rm = TRUE)
      se <- sqrt(1 / sum(w, na.rm = TRUE))
      tanh(zbar + 1.96 * se)
    },
    
    # Combine correlation p-values across outcomes (Stouffer Z; direction-aware)
    p_stouffer = {
      # signed z from each correlation's p-value + sign of r
      p <- as.numeric(Correlation_P_Value)
      s <- sign(r)
      z_i <- qnorm(p / 2, lower.tail = FALSE) * s
      w <- sqrt(pmax(N_eff - 3, 0))
      z_sum <- sum(w * z_i, na.rm = TRUE) / sqrt(sum(w^2, na.rm = TRUE))
      2 * pnorm(abs(z_sum), lower.tail = FALSE)
    }
  ),
  by = .(model_id, variant)
]

write.csv(macro_cor_dt, file.path(output, "joinet_macro_correlations_INT.csv"), row.names = FALSE)


##GLMNET----------------------------------------------------------------------------

output <- "prediction_models/STRING_WSS_v2/model_info/glmnet"
if (!dir.exists(output)) dir.create(output, recursive = TRUE)

alpha_grid <- seq(0, 1, by = 0.2)
K_outer <- 5
lambda_choice <- "lambda.min"  # or "lambda.1se"
# -------------------------
# HELPERS
# -------------------------
prep_X <- function(df, x_cols_raw) {
  X_df <- df[, x_cols_raw, drop = FALSE]
  X <- model.matrix(~ . - 1, data = X_df)  # dummy encode, no intercept
  list(X = X, x_cols_mm = colnames(X))
}

prep_y_numeric <- function(v) {
  if (is.factor(v)) v <- as.character(v)
  as.numeric(v)
}

align_X_test <- function(test_data, x_cols_raw, x_cols_mm_train) {
  X_test_mm <- model.matrix(~ . - 1, data = test_data[, x_cols_raw, drop = FALSE])
  
  missing_mm <- setdiff(x_cols_mm_train, colnames(X_test_mm))
  if (length(missing_mm) > 0) {
    X_test_mm <- cbind(
      X_test_mm,
      matrix(0, nrow(X_test_mm), length(missing_mm),
             dimnames = list(NULL, missing_mm))
    )
  }
  X_test_mm[, x_cols_mm_train, drop = FALSE]
}

# -------------------------
# TRAIN FINAL MODELS (per variant, per outcome)
# -------------------------
final_models <- list()
summary_rows <- list()

for (tag in names(variants)) {
  cat("\n==== Variant:", tag, "====\n")
  
  x_cols_raw <- get_xcols(tag)
  cat(sprintf("[%s] raw predictors: %d\n", tag, length(x_cols_raw)))
  
  # Build base df (keep outcomes + predictors; then complete.cases once)
  df0 <- as.data.frame(pheno2)[, c(outcome, x_cols_raw), drop = FALSE]
  df0 <- df0[complete.cases(df0), , drop = FALSE]
  cat(sprintf("[%s] rows after complete.cases: %d\n", tag, nrow(df0)))
  
  if (nrow(df0) < 10) {
    cat(sprintf("[%s] Too few rows; skipping.\n", tag))
    next
  }
  
  # Build X once per variant
  xobj <- prep_X(df0, x_cols_raw)
  X <- xobj$X
  x_cols_mm <- xobj$x_cols_mm
  
  # Filter finite X rows (paranoia)
  okX <- complete.cases(X) & is.finite(rowSums(X))
  X <- X[okX, , drop = FALSE]
  df0 <- df0[okX, , drop = FALSE]
  
  K_use <- min(K_outer, nrow(X))
  if (K_use < 3) {
    cat(sprintf("[%s] Not enough rows for CV folds; skipping.\n", tag))
    next
  }
  
  # Loop outcomes (single-output cv.glmnet)
  for (outj in outcome) {
    y_raw <- df0[[outj]]
    y <- prep_y_numeric(y_raw)
    
    # drop NA pairs
    ok <- is.finite(y) & complete.cases(X)
    y <- y[ok]
    Xj <- X[ok, , drop = FALSE]
    
    # skip constant / too small
    if (length(y) < 10 || sd(y) == 0) next
    
    # Alpha search
    scores <- data.table(alpha = alpha_grid, mse = NA_real_)
    
    for (i in seq_along(alpha_grid)) {
      a <- alpha_grid[i]
      cvfit <- cv.glmnet(
        x = Xj, y = y,
        family = "gaussian",
        alpha = a,
        nfolds = min(K_use, length(y)),
        type.measure = "mse"
        #standardize = TRUE
      )
      scores$mse[i] <- min(cvfit$cvm, na.rm = TRUE)  # cvm is vector for single outcome
    }
    
    scores <- scores[is.finite(mse)]
    if (nrow(scores) == 0) next
    
    best_alpha <- scores$alpha[which.min(scores$mse)]
    best_mse   <- min(scores$mse)
    
    # Fit final CV model on ALL data at best alpha (stores full-data fit + chosen lambda)
    final_cv <- cv.glmnet(
      x = Xj, y = y,
      family = "gaussian",
      alpha = best_alpha,
      nfolds = min(K_use, length(y)),
      type.measure = "mse"
      #standardize = TRUE
    )
    
    best_lambda <- if (lambda_choice == "lambda.1se") final_cv$lambda.1se else final_cv$lambda.min
    
    model_id <- paste("glmnet", tag, outj, sep = "_")
    
    final_models[[model_id]] <- list(
      model      = final_cv,
      variant    = tag,
      outcome    = outj,
      alpha      = best_alpha,
      lambda     = best_lambda,
      x_cols_raw = x_cols_raw,
      x_cols_mm  = x_cols_mm
    )
    
    summary_rows[[length(summary_rows) + 1]] <- data.table(
      model_id    = model_id,
      variant     = tag,
      outcome     = outj,
      n           = nrow(Xj),
      p           = ncol(Xj),
      best_alpha  = best_alpha,
      best_lambda = best_lambda,
      best_mse    = best_mse
    )
    
    cat(sprintf("[%s|%s] alpha=%.2f lambda=%g mse=%.6f\n", tag, outj, best_alpha, best_lambda, best_mse))
  }
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary_dt, file.path(output, "glmnet_summary_by_variant_outcome_INT.csv"))
saveRDS(final_models, file.path(output, "final_models_glmnet_INT.rds"))

# ============================================================
# PREDICTIONS (per model_id)
# ============================================================
predictions <- list()

for (model_id in names(final_models)) {
  m <- final_models[[model_id]]
  cat("\nProcessing model:", model_id, "\n")
  
  if (!all(m$x_cols_raw %in% colnames(test_data))) {
    missing_cols <- setdiff(m$x_cols_raw, colnames(test_data))
    cat("Missing columns in test data:", paste(missing_cols, collapse = ", "), "\n")
    next
  }
  
  X_test <- align_X_test(test_data, m$x_cols_raw, m$x_cols_mm)
  
  tryCatch({
    pred <- as.numeric(predict(m$model, newx = X_test, s = m$lambda))
    predictions[[model_id]] <- pred
    cat("Processed predictions for model:", model_id, "  n=", length(pred), "\n")
  }, error = function(e) {
    cat("Error in prediction for model", model_id, ":", e$message, "\n")
  })
}

# Save all predictions
saveRDS(predictions, file.path(output, "predictions_glmnet_INT.rds"))

# ============================================================
# METRICS (exactly your style, adapted to single-outcome predictions)
# Assumes actual_values has columns named by outcome
# ============================================================
# 
# som <- Hmisc::somers2(
#   as.numeric(pred_valid),
#   as.numeric(actual_valid)
# )
#
# dxy <- as.numeric(som$Dxy)
# c_index <- (dxy + 1) / 2

metrics <- lapply(names(predictions), function(model_id) {
  pred <- predictions[[model_id]]
  model_info <- final_models[[model_id]]
  outcome <- model_info$outcome
  actual <- actual_values[[outcome]]
  
  valid_idx <- !is.na(pred) & !is.na(actual)
  pred_valid <- pred[valid_idx]
  actual_valid <- actual[valid_idx]
  
  RMSE <- sqrt(mean((pred_valid - actual_valid)^2))
  MAE  <- mean(abs(pred_valid - actual_valid))
  #R2   <- cor(actual_valid, pred_valid)^2
  idx <- which(model_info$model$lambda == model_info$model$lambda.min)
  R2 <- model_info$model$glmnet.fit$dev.ratio[idx]
 
   c_index_result <- tryCatch({
    rcorr.cens(pred_valid, actual_valid)
  }, error = function(e) {
    c("C Index" = NA, Dxy = NA, S.D. = NA, n = NA)
  })
  
  c_index <- c_index_result["C Index"]
  c_index_se <- c_index_result["S.D."]
  
  somers_result <- tryCatch({
    somers2(pred_valid, actual_valid)
  }, error = function(e) {
    c(C = NA, Dxy = NA, S.D. = NA, n = NA)
  })
  
  dxy <- somers_result["Dxy"]
  N <- length(pred_valid)
  
  data.frame(
    model_id = model_id,
    outcome = outcome,
    variant = model_info$variant,
    N = N,
    RMSE = RMSE,
    MAE = MAE,
    R2 = R2,
    C_index = c_index,
    C_index_SE = c_index_se,
    Dxy = dxy,
    row.names = NULL
  )
})

metrics_dt <- rbindlist(metrics, fill = TRUE)
metrics_dt <- metrics_dt[, .(model_id, outcome, variant, N,
                             RMSE, MAE, R2, C_index, C_index_SE, Dxy)]
print(metrics_dt)

metrics_dt_pub <- metrics_dt[, .(
  model_id, outcome, variant, N,
  RMSE = round(RMSE, 3),
  MAE = round(MAE, 3),
  R2 = round(R2, 3),
  C_index = sprintf("%.3f ± %.3f", C_index, C_index_SE)
)]

write.csv(metrics_dt_pub,
          file.path(output, "glmnet_model_metrics_c_index_WSS_INT_updated.csv"),
          row.names = FALSE)

# ============================================================
# CORRELATIONS (your style)
# ============================================================
correlations <- data.frame()

for (model_id in names(predictions)) {
  pred <- predictions[[model_id]]
  model_info <- final_models[[model_id]]
  outcome_i <- model_info$outcome  # local name -- avoid clobbering the global `outcome` vector
  actual <- actual_values[[outcome_i]]

  cor_test <- cor.test(pred, actual, method = "pearson", use = "complete.obs")

  correlations <- rbind(correlations, data.frame(
    model_id = model_id,
    outcome = outcome_i,
    variant = model_info$variant,
    Pearson_Correlation = cor_test$estimate,
    Correlation_P_Value = cor_test$p.value,
    CI_Lower = cor_test$conf.int[1],
    CI_Upper = cor_test$conf.int[2],
    Sample_Size = length(pred)
  ))
}

write.csv(correlations,
          file.path(output, "glmnet_model_correlations_eigengenes_INT.csv"),
          row.names = FALSE)

###########Glmnet (mguassian)#############################################################

# prep_Y <- function(df, outcome_cols) {
#   Y_df <- df[, outcome_cols, drop = FALSE]
#   Y_df <- as.data.frame(lapply(Y_df, function(v) {
#     if (is.factor(v)) v <- as.character(v)
#     as.numeric(v)
#   }))
#   good <- sapply(Y_df, function(v) {
#     v <- v[is.finite(v)]
#     length(v) > 5 && sd(v) > 0
#   })
#   Y_df <- Y_df[, good, drop = FALSE]
#   list(Y = as.matrix(Y_df), outcomes = colnames(Y_df))
# }
# 
# prep_X <- function(df, x_cols_raw) {
#   X_df <- df[, x_cols_raw, drop = FALSE]
#   X <- model.matrix(~ . - 1, data = X_df)  # dummy encode, no intercept
#   list(X = X, x_cols_mm = colnames(X))
# }
# 
# # Robust scorer: works if cvm is matrix (multi-outcome) or vector (single outcome)
# score_cvfit <- function(cvfit) {
#   cvm <- cvfit$cvm
#   
#   # cvm can be:
#   # - matrix: lambdas x outcomes  (multi-response)
#   # - numeric vector: length = #lambdas (single response)
#   if (is.null(cvm)) return(NA_real_)
#   
#   if (is.matrix(cvm)) {
#     rm <- rowMeans(cvm, na.rm = TRUE)
#   } else {
#     rm <- as.numeric(cvm)
#   }
#   
#   rm <- rm[is.finite(rm)]
#   if (length(rm) == 0) return(NA_real_)
#   min(rm)
# }
# 
# # -------------------------
# # MAIN: CV + SAVE FINAL MODELS
# # -------------------------
# final_models <- list()
# summary_rows <- list()
# 
# for (tag in names(variants)) {
#   cat("\n==== Variant:", tag, "====\n")
#   
#   # Must exist in your environment
#   x_cols_raw <- get_xcols(tag)
#   cat(sprintf("[%s] raw predictors: %d\n", tag, length(x_cols_raw)))
#   
#   # Build df and drop NA rows across outcomes + predictors
#   df <- as.data.frame(pheno2)[, c(outcome, x_cols_raw), drop = FALSE]
#   df <- df[complete.cases(df), , drop = FALSE]
#   cat(sprintf("[%s] rows after complete.cases: %d\n", tag, nrow(df)))
#   
#   if (nrow(df) < 10) {
#     cat(sprintf("[%s] Too few rows; skipping.\n", tag))
#     next
#   }
#   
#   # Build Y and X robustly
#   yobj <- prep_Y(df, outcome)
#   Y <- yobj$Y
#   kept_outcomes <- yobj$outcomes
#   if (ncol(Y) < 1) {
#     cat(sprintf("[%s] No valid outcomes after numeric/variance filtering; skipping.\n", tag))
#     next
#   }
#   
#   xobj <- prep_X(df, x_cols_raw)
#   X <- xobj$X
#   x_cols_mm <- xobj$x_cols_mm
#   
#   # Guard: keep only finite rows
#   ok <- complete.cases(X) & complete.cases(Y)
#   X <- X[ok, , drop = FALSE]
#   Y <- Y[ok, , drop = FALSE]
#   
#   cat(sprintf("[%s] Final matrices: n=%d, p=%d, q=%d\n", tag, nrow(X), ncol(X), ncol(Y)))
#   
#   K_use <- min(K_outer, nrow(X))
#   if (K_use < 3) {
#     cat(sprintf("[%s] Not enough rows for CV folds; skipping.\n", tag))
#     next
#   }
#   
#   # ---- Alpha grid search ----
#   scores <- data.table(alpha = alpha_grid, mse = NA_real_)
#   
#   for (i in seq_along(alpha_grid)) {
#     a <- alpha_grid[i]
#     cvfit <- cv.glmnet(
#       x = X, y = Y,
#       family = "mgaussian",
#       alpha = a,
#       nfolds = K_use,
#       type.measure = "mse",
#       standardize = TRUE
#     )
#     scores$mse[i] <- score_cvfit(cvfit)
#     cat(sprintf("[%s] alpha=%.2f  mse=%.6f\n", tag, a, scores$mse[i]))
#   }
#   
#   scores <- scores[is.finite(mse)]
#   if (nrow(scores) == 0) {
#     cat(sprintf("[%s] All alpha scores NA; skipping.\n", tag))
#     next
#   }
#   
#   best_alpha <- scores$alpha[which.min(scores$mse)]
#   best_mse   <- min(scores$mse)
#   cat(sprintf("[%s] BEST alpha=%.2f (mse=%.6f)\n", tag, best_alpha, best_mse))
#   
#   # ---- Fit final CV model on ALL data at best alpha ----
#   final_cv <- cv.glmnet(
#     x = X, y = Y,
#     family = "mgaussian",
#     alpha = best_alpha,
#     nfolds = K_use,
#     type.measure = "mse",
#     standardize = TRUE
#   )
#   
#   best_lambda <- if (lambda_choice == "lambda.1se") final_cv$lambda.1se else final_cv$lambda.min
#   
#   model_id <- paste("glmnet_mgaussian", tag, sep = "_")
#   
#   final_models[[model_id]] <- list(
#     model      = final_cv,
#     variant    = tag,
#     outcomes   = kept_outcomes,
#     alpha      = best_alpha,
#     lambda     = best_lambda,
#     x_cols_raw = x_cols_raw,
#     x_cols_mm  = x_cols_mm,
#     n          = nrow(X),
#     p          = ncol(X),
#     q          = ncol(Y)
#   )
#   
#   summary_rows[[length(summary_rows) + 1]] <- data.table(
#     model_id     = model_id,
#     variant      = tag,
#     n            = nrow(X),
#     p            = ncol(X),
#     q            = ncol(Y),
#     best_alpha   = best_alpha,
#     best_lambda  = best_lambda,
#     best_mse     = best_mse
#   )
# }
# 
# summary_dt <- rbindlist(summary_rows, fill = TRUE)
# print(summary_dt)
# 
# saveRDS(final_models, file.path(output, "final_models_glmnet_mgaussian_z_scale.rds"))
# fwrite(summary_dt, file.path(output, "glmnet_mgaussian_summary_by_variant_z_scale.csv"))
# 
# 
# predictions <- list()
# 
# for (model_id in names(final_models)) {
#   m <- final_models[[model_id]]
#   cat("\nProcessing model:", model_id, "\n")
#   
#   if (!all(m$x_cols_raw %in% colnames(test_data))) {
#     missing_cols <- setdiff(m$x_cols_raw, colnames(test_data))
#     cat("Missing columns in test_data:", paste(missing_cols, collapse = ", "), "\n")
#     next
#   }
#   
#   X_test <- align_X_test(test_data, m$x_cols_raw, m$x_cols_mm)
#   
#   tryCatch({
#     # Use the correct lambda reference from your model structure
#     pred_raw <- predict(m$model, newx = X_test, s = m$lambda, type = "response")
#     
#     # For mgaussian, predictions come as a 3D array: [observations, outcomes, lambda]
#     # Extract predictions - remove the lambda dimension
#     if (length(dim(pred_raw)) == 3) {
#       # 3D array case
#       pred_matrix <- pred_raw[, , 1]  # Remove the lambda dimension
#     } else {
#       # 2D matrix case (if only one lambda value was used)
#       pred_matrix <- pred_raw
#     }
#     
#     # Convert to data frame
#     pred_df <- as.data.frame(pred_matrix)
#     
#     predictions[[model_id]] <- pred_df
#     cat("Processed predictions for model:", model_id, "  shape:", dim(pred_df), "\n")
#     
#   }, error = function(e) {
#     cat("Error in prediction for model", model_id, ":", e$message, "\n")
#     # Debug info
#     cat("Lambda value:", m$lambda, "\n")
#     cat("Lambda class:", class(m$lambda), "\n")
#     cat("Available lambdas in model:", head(m$model$lambda), "\n")
#   })
# }
# 
# 
# saveRDS(predictions, file.path(output, "predictions_glmnet_z_scale.rds"))
# #GLmnet (mgaussian metrics) -----------------------------------------------------------
# 
# metrics <- list()
# counter <- 1
# 
# for (model_id in names(predictions)) {
#   pred_df <- predictions[[model_id]]  # This is a data frame with multiple outcome columns
#   model_info <- final_models[[model_id]]
#   
#   # Loop through each outcome column in the predictions
#   for (outcome_col in colnames(pred_df)) {
#     pred <- pred_df[[outcome_col]]  # Extract predictions for this specific outcome
#     
#     # Get actual values for this outcome
#     if (outcome_col %in% names(actual_values)) {
#       actual <- actual_values[[outcome_col]]
#     } else {
#       cat("Warning: No actual values found for outcome:", outcome_col, "\n")
#       next
#     }
#     
#     # Ensure same length
#     min_length <- min(length(pred), length(actual))
#     pred <- pred[1:min_length]
#     actual <- actual[1:min_length]
#     
#     # Remove missing values
#     valid_idx <- !is.na(pred) & !is.na(actual)
#     pred_valid <- pred[valid_idx]
#     actual_valid <- actual[valid_idx]
#     
#     if (length(pred_valid) < 2) {
#       cat("Warning: Not enough valid predictions for", model_id, outcome_col, "\n")
#       next
#     }
#     
#     # Calculate metrics
#     RMSE <- sqrt(mean((pred_valid - actual_valid)^2))
#     MAE  <- mean(abs(pred_valid - actual_valid))
#     #R2   <- cor(actual_valid, pred_valid)^2
#     idx <- which(model_info$model$lambda == model_info$model$lambda.min)
#     R2 <- model_info$model$glmnet.fit$dev.ratio[idx]
#     
#     # C-index and Somers' D
#     c_index_result <- tryCatch({
#       Hmisc::rcorr.cens(pred_valid, actual_valid)
#     }, error = function(e) {
#       c("C Index" = NA, "Dxy" = NA, "S.D." = NA, "n" = NA)
#     })
#     
#     somers_result <- tryCatch({
#       Hmisc::somers2(pred_valid, actual_valid)
#     }, error = function(e) {
#       c("C" = NA, "Dxy" = NA, "S.D." = NA, "n" = NA)
#     })
#     
#     # Extract values
#     c_index <- as.numeric(c_index_result["C Index"])
#     c_index_se <- as.numeric(c_index_result["S.D."])
#     dxy <- as.numeric(somers_result["Dxy"])
#     bias <- mean(pred_valid - actual_valid)
#     n <- length(pred_valid)
#     
#     # Store metrics
#     metrics[[counter]] <- data.frame(
#       model_id = model_id,
#       outcome = outcome_col,
#       variant = model_info$variant,
#       N = n,
#       RMSE = RMSE,
#       MAE = MAE,
#       R2 = R2,
#       C_index = c_index,
#       C_index_SE = c_index_se,
#       Dxy = dxy,
#       Bias = bias,
#       stringsAsFactors = FALSE
#     )
#     
#     counter <- counter + 1
#   }
# }
# 
# # Combine all metrics
# library(data.table)
# metrics_dt <- rbindlist(metrics, fill = TRUE)
# 
# metrics_dt_pub <- metrics_dt[, .(
#   model_id, outcome, variant, N,
#   RMSE = round(RMSE, 3),
#   MAE = round(MAE, 3),
#   R2 = round(R2, 3),
#   C_index = sprintf("%.3f ± %.3f", C_index, C_index_SE)
# )]
# 
# write.csv(metrics_dt_pub,
#           file.path(output, "glmnet_model_metrics_c_index_eigengenes_gmaussian_z_scale_V2.csv"),
#           row.names = FALSE)
# 		  
# 		  
# correlations <- data.frame()
# 
# for (model_id in names(predictions)) {
#   pred_df <- predictions[[model_id]]  # This is a data frame with multiple outcome columns
#   model_info <- final_models[[model_id]]
#   
#   # Loop through each outcome column in the predictions
#   for (outcome_col in colnames(pred_df)) {
#     pred <- pred_df[[outcome_col]]  # Extract predictions for this specific outcome
#     
#     # Get actual values for this outcome
#     if (outcome_col %in% names(actual_values)) {
#       actual <- actual_values[[outcome_col]]
#     } else {
#       cat("Warning: No actual values found for outcome:", outcome_col, "\n")
#       next
#     }
#     
#     # Ensure same length
#     min_length <- min(length(pred), length(actual))
#     pred <- pred[1:min_length]
#     actual <- actual[1:min_length]
#     
#     # Check if we have enough valid pairs for correlation
#     valid_pairs <- complete.cases(pred, actual)
#     if (sum(valid_pairs) < 3) {
#       cat("Warning: Not enough complete cases for correlation -", model_id, outcome_col, "\n")
#       next
#     }
#     
#     # Calculate correlation
#     cor_test <- tryCatch({
#       cor.test(pred, actual, method = "pearson", use = "complete.obs")
#     }, error = function(e) {
#       cat("Error in correlation for", model_id, outcome_col, ":", e$message, "\n")
#       return(NULL)
#     })
#     
#     if (!is.null(cor_test)) {
#       correlations <- rbind(correlations, data.frame(
#         model_id = model_id,
#         outcome = outcome_col,
#         variant = model_info$variant,
#         Pearson_Correlation = cor_test$estimate,
#         Correlation_P_Value = cor_test$p.value,
#         CI_Lower = cor_test$conf.int[1],
#         CI_Upper = cor_test$conf.int[2],
#         Sample_Size = sum(valid_pairs),  # Number of complete cases used
#         stringsAsFactors = FALSE
#       ))
#     }
#   }
# }
# 
# # Clean up row names
# rownames(correlations) <- NULL
# 
# # Print results
# print(correlations)
# 
# 
# write.csv(correlations,
#           file.path(output, "glmnet_model_correlations_eigengenes_mgaussian_z_scale.csv"),
#           row.names = FALSE)
# 

#Individudal Random forest

output <- "prediction_models/STRING_WSS_v2/model_info/rf"
if (!dir.exists(output)) dir.create(output, recursive = TRUE)


variants <- list(
  Biomarkers_Only = c("HybControlNormScale","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8", "APOE.e34", "APOE.e44", "Age", "Gender"),
  Biomarkers_Age_Sex = c("HybControlNormScale","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8", "APOE.e34", "APOE.e44"),
  Age_Sex_HybridNormControl = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8", "APOE.e34", "APOE.e44" ),
  Age_Sex_HybridNormControl_PCs = c("APOE.e34", "APOE.e44"),
  ALL = character(0)
)

predictor_pool <- setdiff(names(pheno2), outcome)

# CV settings 
K_outer <- 5
K_inner <- 10

# Ranger hyperparameter grid
mtry_grid <- c(2, 3, 5, 7, 10)
num_trees_grid <- c(500, 1000, 1500)  # Including your preferred 1500
min_node_size_grid <- c(1, 3, 5)

# Set up parallel backend
cores <- detectCores() - 1
cl <- makeCluster(cores)
registerDoParallel(cl)

# Export necessary functions and data to the parallel workers
clusterExport(cl, c("get_xcols", "pheno2", "outcome", "K_outer", "K_inner", "variants", "predictor_pool"))

# Custom CV function for Ranger
cv_ranger <- function(X, y, mtry, num_trees, min_node_size, nfolds = 5) {
  set.seed(42)
  fold_ids <- sample(rep(1:nfolds, length.out = nrow(X)))
  mse_folds <- numeric(nfolds)
  
  for (k in 1:nfolds) {
    train_idx <- fold_ids != k
    test_idx <- fold_ids == k
    
    train_data <- data.frame(y = y[train_idx], X[train_idx, , drop = FALSE])
    test_data <- data.frame(X[test_idx, , drop = FALSE])
    
    # Fit Ranger model
    rf_model <- ranger(
      y ~ ., 
      data = train_data,
      mtry = mtry,
      num.trees = num_trees,
      min.node.size = min_node_size,
      sample.fraction = 0.50,
      replace = TRUE
    )
    
    # Predict and calculate MSE
    pred <- predict(rf_model, test_data)$predictions
    mse_folds[k] <- mean((pred - y[test_idx])^2)
  }
  
  return(mean(mse_folds))
}

final_models <- list()

for (tag in names(variants)) {
  cat("\n==== Variant:", tag, "====\n")
  
  # 1) Build data for this variant using your function
  x_cols <- get_xcols(tag)
  
  df <- as.data.frame(pheno2)[, c(outcome, x_cols), drop = FALSE]
  cat(sprintf("[%s] Dimensions of initial dataframe: %d rows, %d columns\n", tag, nrow(df), ncol(df)))
  
  df <- df[complete.cases(df), ]
  cat(sprintf("[%s] Dimensions after removing NA: %d rows, %d columns\n", tag, nrow(df), ncol(df)))
  
  if (nrow(df) < max(K_outer, 10)) {
    cat(sprintf("[%s] Too few rows after NA filtering (%d). Skipping.\n", tag, nrow(df)))
    next
  }
  
  X <- as.matrix(df[, x_cols, drop = FALSE])
  cat(sprintf("[%s] Data prepared. Rows: %d, Predictors: %d\n", tag, nrow(X), ncol(X)))
  
  # Adjust mtry_grid based on number of predictors
  max_mtry <- min(ncol(X), max(mtry_grid))
  mtry_use <- mtry_grid[mtry_grid <= max_mtry]
  if (length(mtry_use) == 0) mtry_use <- c(1, max(1, floor(sqrt(ncol(X)))))
  
  # Loop through each outcome
  for (outcome_i in outcome) {
    cat(sprintf("\n[%s - %s] Processing outcome\n", tag, outcome_i))
    
    y <- df[[outcome_i]]
    
    # 2) Find best hyperparameters in parallel
    set.seed(42)
    grid <- CJ(mtry = mtry_use, num_trees = num_trees_grid, min_node_size = min_node_size_grid)
    
    results <- tryCatch({
      foreach(i = 1:nrow(grid), .combine = rbind, .packages = c("data.table", "ranger"), .errorhandling = 'pass') %dopar% {
        mt <- grid$mtry[i]
        nt <- grid$num_trees[i]
        mns <- grid$min_node_size[i]
        
        mse <- cv_ranger(X, y, mtry = mt, num_trees = nt, min_node_size = mns, nfolds = K_outer)
        
        data.table(mtry = mt, num_trees = nt, min_node_size = mns, mse = mse)
      }
    }, error = function(e) {
      cat(sprintf("[%s - %s] Error in parallel execution: %s\n", tag, outcome_i, e$message))
      return(NULL)
    })
    
    if (is.null(results) || nrow(results) == 0) {
      cat(sprintf("[%s - %s] No valid results. Skipping this outcome.\n", tag, outcome_i))
      next
    }
    
    # Check for errors in results
    error_rows <- which(sapply(results, inherits, "error"))
    if (length(error_rows) > 0) {
      cat(sprintf("[%s - %s] Errors occurred in %d iterations\n", tag, outcome_i, length(error_rows)))
      results <- results[-error_rows, ]
    }
    
    if (nrow(results) == 0) {
      cat(sprintf("[%s - %s] No valid results after error removal. Skipping this outcome.\n", tag, outcome_i))
      next
    }
    
    best_result <- results[which.min(mse)]
    best_mtry <- best_result$mtry
    best_ntree <- best_result$num_trees
    best_nodesize <- best_result$min_node_size
    
    cat(sprintf("[%s - %s] Best params: mtry=%d, ntree=%d, nodesize=%d, MSE=%.4f\n", 
                tag, outcome_i, best_mtry, best_ntree, best_nodesize, best_result$mse))
    
    # 3) Create and save final model for this variant-outcome combination
    model_id <- paste("rf", tag, outcome_i, sep = "_")
    
    final_data <- data.frame(y = y, X)
    final_model <- ranger(
      y ~ ., 
      data = final_data,
      mtry = best_mtry,
      num.trees = best_ntree,
      min.node.size = best_nodesize,
      sample.fraction = 0.50,
      replace = TRUE,
      importance = "impurity"
    )
    
    final_models[[model_id]] <- list(
      model = final_model,
      outcome = outcome_i,
      variant = tag,
      mtry = best_mtry,
      num_trees = best_ntree,
      min_node_size = best_nodesize,
      mse_cv = best_result$mse,
      x_cols = x_cols
    )
    
    cat(sprintf("Processed model: %s\n", model_id))
  }
}

# Stop cluster
stopCluster(cl)

# Save all final models
saveRDS(final_models, file.path(output, "final_models_RF_INT.rds"))

# Create a summary table
summary_table <- data.table(
  Model_ID = names(final_models),
  Outcome = sapply(final_models, function(m) m$outcome),
  Variant = sapply(final_models, function(m) m$variant),
  mtry = sapply(final_models, function(m) m$mtry),
  num_trees = sapply(final_models, function(m) m$num_trees),
  min_node_size = sapply(final_models, function(m) m$min_node_size),
  CV_MSE = sapply(final_models, function(m) m$mse_cv),
  N_Features = sapply(final_models, function(m) length(m$x_cols))
)

# Save summary table
fwrite(summary_table, file.path(output, "final_models_RF_summary_WSS_INT.csv"))

print(summary_table)

predictions <- list()

# Loop through models and make predictions
for (model_id in names(final_models)) {
  model_info <- final_models[[model_id]]
  
  cat("\nProcessing model:", model_id, "\n")
  
  # Check if x_cols exist in test_data
  if (!all(model_info$x_cols %in% colnames(test_data))) {
    missing_cols <- setdiff(model_info$x_cols, colnames(test_data))
    cat("Missing columns in test data:", paste(missing_cols, collapse = ", "), "\n")
    next
  }
  
  # Prepare test data for this model
  X_test <- test_data[, model_info$x_cols, drop = FALSE]
  
  # Make predictions
  tryCatch({
    # For Random Forest - no 'newx', just direct prediction
    pred <- predict(model_info$model, X_test)$predictions
    
    # RF returns a vector, not a matrix like joinet
    # Create a matrix with single column for consistency with your downstream code
    pred_matrix <- as.matrix(pred)
    colnames(pred_matrix) <- model_info$outcome  # Single outcome name
    
    # Store predictions
    predictions[[model_id]] <- pred_matrix
    
    cat("Processed predictions for model:", model_id, "\n")
    cat("Prediction dimensions:", paste(dim(pred_matrix), collapse = " x "), "\n")
    cat("Predicted outcome:", model_info$outcome, "\n")
  }, error = function(e) {
    cat("Error in prediction for model", model_id, ":", e$message, "\n")
  })
}

# Since RF models are variant-outcome combinations, organize predictions by variant
# Group predictions by variant
predictions_by_variant <- list()

for (model_id in names(predictions)) {
  variant <- final_models[[model_id]]$variant
  outcome_i <- final_models[[model_id]]$outcome  # local name -- avoid clobbering the global `outcome` vector

  if (!variant %in% names(predictions_by_variant)) {
    predictions_by_variant[[variant]] <- list()
  }

  predictions_by_variant[[variant]][[outcome_i]] <- predictions[[model_id]][,1]  # Extract vector
}

# Convert to data frames and save by variant
for (variant in names(predictions_by_variant)) {
  variant_df <- as.data.frame(predictions_by_variant[[variant]])
  variant_df <- cbind(SampleID = rownames(test_data), variant_df)  # Use test_data rownames
  
  filename <- paste0(variant, ".csv")
  fwrite(variant_df, file.path(output, filename))
  cat("Saved predictions for variant:", variant, "\n")
}

# Calculate metrics - adapted for single-outcome RF models
metrics <- lapply(names(predictions), function(model_id) {
  predicted <- predictions[[model_id]]
  model_info <- final_models[[model_id]]
  
  outcome <- model_info$outcome  # singular, matching how RF models are stored
  pred <- predicted[, 1]
  actual <- actual_values[[outcome]]
  
  valid_idx <- !is.na(pred) & !is.na(actual)
  pred_valid <- pred[valid_idx]
  actual_valid <- actual[valid_idx]
  
  RMSE <- sqrt(mean((pred_valid - actual_valid)^2))
  MAE  <- mean(abs(pred_valid - actual_valid))
  R2   <- cor(actual_valid, pred_valid)^2
  
  c_index_result <- tryCatch({
    rcorr.cens(pred_valid, actual_valid)
  }, error = function(e) c(`C Index` = NA, `S.D.` = NA))
  
  somers_result <- tryCatch({
    somers2(pred_valid, actual_valid)
  }, error = function(e) c(Dxy = NA))
  
  data.frame(
    model_id = model_id,
    outcome = outcome,
    variant = model_info$variant,
    N = length(pred_valid),
    RMSE = RMSE,
    MAE = MAE,
    R2 = R2,
    C_index = as.numeric(c_index_result["C Index"]),
    C_index_SE = as.numeric(c_index_result["S.D."]),
    Dxy = as.numeric(somers_result["Dxy"]),
    Bias = mean(pred_valid - actual_valid),
    row.names = NULL
  )
})
metrics_dt <- rbindlist(metrics, fill = TRUE)
# Reorder columns
metrics_dt <- metrics_dt[, .(model_id, outcome, variant,
                             RMSE, MAE, R2, C_index, C_index_SE, Dxy, Bias)]

print(metrics_dt)

# Format for publication with confidence intervals
metrics_dt_pub <- metrics_dt[, .(
  model_id, outcome, variant,
  RMSE = round(RMSE, 3),
  MAE = round(MAE, 3),
  R2 = round(R2, 3),
  C_index = sprintf("%.3f ± %.3f", C_index, C_index_SE),
  Dxy = round(Dxy, 3),
  Bias = round(Bias, 3)
)]

print(metrics_dt_pub)


write.csv(metrics_dt_pub, file.path(output, "rf_model_metrics_c_index_WSS_INT_2.csv"), row.names = FALSE)

## correlations 

correlations <- data.frame()

for (model_id in names(predictions)) {
  pred_matrix <- predictions[[model_id]]
  model_info <- final_models[[model_id]]
  
  # RF models have single outcome, not multiple outcomes
  outcome_i <- model_info$outcome  # local name -- avoid clobbering the global `outcome` vector
  pred <- pred_matrix[, 1]       # Single column (not multiple columns)

  # Get actual values for this outcome
  if (outcome_i %in% names(actual_values)) {
    actual <- actual_values[[outcome_i]]

    # Ensure same length
    min_length <- min(length(pred), length(actual))
    pred <- pred[1:min_length]
    actual <- actual[1:min_length]

    # Check for sufficient valid pairs
    valid_pairs <- complete.cases(pred, actual)
    if (sum(valid_pairs) >= 3) {

      cor_test <- tryCatch({
        cor.test(pred, actual, method = "pearson", use = "complete.obs")
      }, error = function(e) {
        cat("Error in correlation for", model_id, outcome_i, ":", e$message, "\n")
        return(NULL)
      })

      if (!is.null(cor_test)) {
        correlations <- rbind(correlations, data.frame(
          model_id = model_id,
          outcome = outcome_i,
          variant = model_info$variant,
          Pearson_Correlation = cor_test$estimate,
          Correlation_P_Value = cor_test$p.value,
          CI_Lower = cor_test$conf.int[1],
          CI_Upper = cor_test$conf.int[2],
          Sample_Size = sum(valid_pairs),  # Number of complete cases used
          stringsAsFactors = FALSE
        ))
      }
    } else {
      cat("Warning: Not enough complete cases for correlation -", model_id, outcome_i, "\n")
    }
  } else {
    cat("Warning: No actual values found for outcome:", outcome, "\n")
  }
}

# Clean up row names
rownames(correlations) <- NULL

# Print results
print(correlations)

# Save correlations
write.csv(correlations, file.path(output, "rf_model_correlations_WSS_INT.csv"), row.names = FALSE)


##mrf----------------------------------------------------------------------------------------------------------------------------- 

library(randomForestSRC)
library(parallel)
library(doParallel)
library(foreach)
library(data.table)

output <- "prediction_models/STRING_WSS_v2/model_info/mrf"
if (!dir.exists(output)) dir.create(output, recursive = TRUE)


variants <- list(
  Biomarkers_Only = c("HybControlNormScale","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8", "APOE.e34", "APOE.e44", "Age", "Gender"),
  Biomarkers_Age_Sex = c("HybControlNormScale","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8", "APOE.e34", "APOE.e44"),
  Age_Sex_HybridNormControl = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8", "APOE.e34", "APOE.e44"),
  Age_Sex_HybridNormControl_PCs = c("APOE.e34", "APOE.e44"),
  ALL = character(0)
)

predictor_pool <- setdiff(names(pheno2), outcome)

# CV settings 
K_outer <- 5
K_inner <- 10

# rfsrc hyperparameter grid
mtry_grid <- c(2, 3, 5, 7, 10)
num_trees_grid <- c(500, 1000, 1500)
min_node_size_grid <- c(1, 3, 5)

# Set up parallel backend -- capped lower than other sections: each worker
# fits a 75-outcome multivariate rfsrc model, which is far more memory-heavy
# per-worker than joinet/glmnet/ranger's single/few-outcome fits, and running
# detectCores()-1 workers here was OOM-killing the process.
cores <- min(6, detectCores() - 1)
cl <- makeCluster(cores)
registerDoParallel(cl)

# Export necessary functions and data to the parallel workers
clusterExport(cl, c("get_xcols", "pheno2", "outcome", "K_outer", "K_inner",
                    "variants", "predictor_pool"))

# Custom CV function for Multivariate rfsrc
cv_rfsrc <- function(data, outcome_vars, x_vars, mtry, num_trees, min_node_size, nfolds = 5) {
  set.seed(42)
  fold_ids <- sample(rep(1:nfolds, length.out = nrow(data)))
  mse_folds <- numeric(nfolds)
  
  # Create multivariate formula
  formula_str <- paste("Multivar(", paste(outcome_vars, collapse = ", "), ") ~", 
                       paste(x_vars, collapse = " + "))
  mv_formula <- as.formula(formula_str)
  
  for (k in 1:nfolds) {
    train_idx <- fold_ids != k
    test_idx <- fold_ids == k
    
    train_data <- data[train_idx, ]
    test_data <- data[test_idx, ]
    
    # Fit Multivariate rfsrc model
    rf_model <- tryCatch({
      rfsrc(
        mv_formula,
        data = train_data,
        mtry = mtry,
        ntree = num_trees,
        nodesize = min_node_size,
        sampsize = function(x) ceiling(0.632 * nrow(train_data)),  # ~63.2% bootstrap
        importance = "none",  # Skip importance for CV speed
        forest = TRUE
      )
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(rf_model)) {
      mse_folds[k] <- NA
      next
    }
    
    # Predict on test set
    pred_obj <- tryCatch({
      predict(rf_model, newdata = test_data)
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(pred_obj)) {
      mse_folds[k] <- NA
      next
    }
    
    # Calculate MSE across all outcomes
    # pred_obj$regrOutput is a list with predictions for each outcome
    mse_outcome <- numeric(length(outcome_vars))
    for (i in seq_along(outcome_vars)) {
      y_true <- test_data[[outcome_vars[i]]]
      y_pred <- pred_obj$regrOutput[[outcome_vars[i]]]$predicted
      mse_outcome[i] <- mean((y_true - y_pred)^2, na.rm = TRUE)
    }
    mse_folds[k] <- mean(mse_outcome)
  }
  
  return(mean(mse_folds, na.rm = TRUE))
}

final_models <- list()

for (tag in names(variants)) {
  cat("\n==== Variant:", tag, "====\n")
  
  # 1) Build data for this variant
  x_cols <- get_xcols(tag)
  
  df <- as.data.frame(pheno2)[, c(outcome, x_cols), drop = FALSE]
  cat(sprintf("[%s] Dimensions of initial dataframe: %d rows, %d columns\n", tag, nrow(df), ncol(df)))
  
  df <- df[complete.cases(df), ]
  cat(sprintf("[%s] Dimensions after removing NA: %d rows, %d columns\n", tag, nrow(df), ncol(df)))
  
  if (nrow(df) < max(K_outer, 10)) {
    cat(sprintf("[%s] Too few rows after NA filtering (%d). Skipping.\n", tag, nrow(df)))
    next
  }
  
  cat(sprintf("[%s] Data prepared. Rows: %d, Predictors: %d, Outcomes: %d\n", 
              tag, nrow(df), length(x_cols), length(outcome)))
  
  # Adjust mtry_grid based on number of predictors
  max_mtry <- min(length(x_cols), max(mtry_grid))
  mtry_use <- mtry_grid[mtry_grid <= max_mtry]
  if (length(mtry_use) == 0) mtry_use <- c(1, max(1, floor(sqrt(length(x_cols)))))
  
  cat(sprintf("\n[%s] Tuning multivariate random forest with rfsrc\n", tag))
  
  # 2) Find best hyperparameters in parallel
  set.seed(42)
  grid <- CJ(mtry = mtry_use, num_trees = num_trees_grid, min_node_size = min_node_size_grid)
  
  results <- tryCatch({
    foreach(i = 1:nrow(grid), .combine = rbind, 
            .packages = c("data.table", "randomForestSRC"), 
            .errorhandling = 'pass') %dopar% {
      mt <- grid$mtry[i]
      nt <- grid$num_trees[i]
      mns <- grid$min_node_size[i]
      
      mse <- cv_rfsrc(df, outcome_vars = outcome, x_vars = x_cols,
                      mtry = mt, num_trees = nt, min_node_size = mns, 
                      nfolds = K_outer)
      
      data.table(mtry = mt, num_trees = nt, min_node_size = mns, mse = mse)
    }
  }, error = function(e) {
    cat(sprintf("[%s] Error in parallel execution: %s\n", tag, e$message))
    return(NULL)
  })
  
  if (is.null(results) || nrow(results) == 0) {
    cat(sprintf("[%s] No valid results. Skipping this variant.\n", tag))
    next
  }
  
  # Check for errors in results
  error_rows <- which(sapply(results$mse, function(x) is.na(x) || is.infinite(x)))
  if (length(error_rows) > 0) {
    cat(sprintf("[%s] Errors occurred in %d iterations\n", tag, length(error_rows)))
    results <- results[-error_rows, ]
  }
  
  if (nrow(results) == 0) {
    cat(sprintf("[%s] No valid results after error removal. Skipping this variant.\n", tag))
    next
  }
  
  best_result <- results[which.min(mse)]
  best_mtry <- best_result$mtry
  best_ntree <- best_result$num_trees
  best_nodesize <- best_result$min_node_size
  
  cat(sprintf("[%s] Best params: mtry=%d, ntree=%d, nodesize=%d, MSE=%.4f\n", 
              tag, best_mtry, best_ntree, best_nodesize, best_result$mse))
  
  # 3) Create and save final model for this variant
  model_id <- paste("rfsrc", tag, sep = "_")
  
  # Build final model on full data with importance
  formula_str <- paste("Multivar(", paste(outcome, collapse = ", "), ") ~", 
                       paste(x_cols, collapse = " + "))
  mv_formula <- as.formula(formula_str)
  
  final_model <- tryCatch({
    rfsrc(
      mv_formula,
      data = df,
      mtry = best_mtry,
      ntree = best_ntree,
      nodesize = best_nodesize,
      sampsize = function(x) ceiling(0.632 * nrow(df)),
      importance = "permute",  # Get variable importance
      forest = TRUE
    )
  }, error = function(e) {
    cat(sprintf("[%s] Error building final model: %s\n", tag, e$message))
    return(NULL)
  })
  
  if (!is.null(final_model)) {
    final_models[[model_id]] <- list(
      model = final_model,
      outcomes = outcome,
      variant = tag,
      mtry = best_mtry,
      num_trees = best_ntree,
      min_node_size = best_nodesize,
      mse_cv = best_result$mse,
      x_cols = x_cols,
      importance = final_model$importance,
      err_rate = final_model$err.rate,
      predicted_oob = final_model$predicted.oob
    )
    
    cat(sprintf("Processed multivariate model: %s\n", model_id))
  }
}

# Stop cluster
stopCluster(cl)

# Display summary
cat("\n=== Final Models Summary ===\n")
for (model_name in names(final_models)) {
  m <- final_models[[model_name]]
  cat(sprintf("\nModel: %s\n", model_name))
  cat(sprintf("  Variant: %s\n", m$variant))
  cat(sprintf("  Outcomes: %s\n", paste(m$outcomes, collapse = ", ")))
  cat(sprintf("  Predictors: %d\n", length(m$x_cols)))
  cat(sprintf("  CV MSE: %.4f\n", m$mse_cv))
  cat(sprintf("  Best params: mtry=%d, ntrees=%d, nodesize=%d\n", 
              m$mtry, m$num_trees, m$min_node_size))
  
  # Show top 5 important variables (averaged across outcomes)
  if (!is.null(m$importance)) {
    avg_imp <- rowMeans(m$importance, na.rm = TRUE)
    top5 <- head(sort(avg_imp, decreasing = TRUE), 5)
    cat("  Top 5 important variables:\n")
    for (var in names(top5)) {
      cat(sprintf("    %s: %.4f\n", var, top5[var]))
    }
  }
}

saveRDS(final_models, file = "final_rfsrc_models.rds")

model_summaries <- do.call(rbind, lapply(names(final_models), function(model_name) {
  m <- final_models[[model_name]]
  data.frame(
    Model = model_name,
    Variant = m$variant,
    Outcomes = paste(m$outcomes, collapse = ", "),
    N_Predictors = length(m$x_cols),
    CV_MSE = m$mse_cv,
    mtry = m$mtry,
    ntrees = m$num_trees,
    nodesize = m$min_node_size,
    stringsAsFactors = FALSE
  )
}))

write.csv(model_summaries, "mrf_Wss_model_summaries.csv", row.names = FALSE)

## predictions for mrf 

predictions <- list()

# Loop through models and make predictions
for (model_id in names(final_models)) {
  model_info <- final_models[[model_id]]
  
  cat("\nProcessing model:", model_id, "\n")
  
  # Check if x_cols exist in test_data
  if (!all(model_info$x_cols %in% colnames(test_data))) {
    missing_cols <- setdiff(model_info$x_cols, colnames(test_data))
    cat("Missing columns in test data:", paste(missing_cols, collapse = ", "), "\n")
    next
  }
  
  # Prepare test data for this model
  test_subset <- test_data[, model_info$x_cols, drop = FALSE]
  
  # Make predictions
  tryCatch({
    pred_obj <- predict(model_info$model, newdata = test_subset)
    
    # Extract predictions for each outcome
    # For rfsrc multivariate, predictions are in pred_obj$regrOutput
    pred_matrix <- matrix(NA, nrow = nrow(test_subset), ncol = length(model_info$outcomes))
    colnames(pred_matrix) <- model_info$outcomes
    
    for (i in seq_along(model_info$outcomes)) {
      outcome_name <- model_info$outcomes[i]
      pred_matrix[, i] <- pred_obj$regrOutput[[outcome_name]]$predicted
    }
    
    # Store predictions
    predictions[[model_id]] <- pred_matrix
    
    cat("Processed predictions for model:", model_id, "\n")
    cat("Prediction dimensions:", paste(dim(pred_matrix), collapse = " x "), "\n")
    cat("Predicted outcomes:", paste(colnames(pred_matrix), collapse = ", "), "\n")
  }, error = function(e) {
    cat("Error in prediction for model", model_id, ":", e$message, "\n")
  })
}

# Save individual prediction files
# Get actual sample IDs from test_data
sample_ids <- rownames(test_data)

# Biomarkers_Only
if ("rfsrc_Biomarkers_Only" %in% names(predictions)) {
  Biomarkers_Only <- as.data.frame(predictions[["rfsrc_Biomarkers_Only"]])
  Biomarkers_Only <- cbind(SampleID = sample_ids, Biomarkers_Only)
  fwrite(Biomarkers_Only, file.path(output, "Biomarkers_Only.csv"))
}

# Biomarkers_Age_Sex
if ("rfsrc_Biomarkers_Age_Sex" %in% names(predictions)) {
  Biomarkers_Age_Sex <- as.data.frame(predictions[["rfsrc_Biomarkers_Age_Sex"]])
  Biomarkers_Age_Sex <- cbind(SampleID = sample_ids, Biomarkers_Age_Sex)
  fwrite(Biomarkers_Age_Sex, file.path(output, "Biomarkers_Age_Sex.csv"))
}

# Age_Sex_HybridNormControl
if ("rfsrc_Age_Sex_HybridNormControl" %in% names(predictions)) {
  Age_Sex_HybridNormControl <- as.data.frame(predictions[["rfsrc_Age_Sex_HybridNormControl"]])
  Age_Sex_HybridNormControl <- cbind(SampleID = sample_ids, Age_Sex_HybridNormControl)
  fwrite(Age_Sex_HybridNormControl, file.path(output, "Age_Sex_HybridNormControl.csv"))
}

# Age_Sex_HybridNormControl_PCs
if ("rfsrc_Age_Sex_HybridNormControl_PCs" %in% names(predictions)) {
  Age_Sex_HybridNormControl_PCs <- as.data.frame(predictions[["rfsrc_Age_Sex_HybridNormControl_PCs"]])
  Age_Sex_HybridNormControl_PCs <- cbind(SampleID = sample_ids, Age_Sex_HybridNormControl_PCs)
  fwrite(Age_Sex_HybridNormControl_PCs, file.path(output, "Age_Sex_HybridNormControl_PCs.csv"))
}

# ALL
if ("rfsrc_ALL" %in% names(predictions)) {
  ALL <- as.data.frame(predictions[["rfsrc_ALL"]])
  ALL <- cbind(SampleID = sample_ids, ALL)
  fwrite(ALL, file.path(output, "ALL.csv"))
}

# Calculate metrics
metrics <- lapply(names(predictions), function(model_id) {
  predicted <- predictions[[model_id]]
  model_info <- final_models[[model_id]]
  
  outcome_metrics <- lapply(seq_along(model_info$outcomes), function(i) {
    outcome <- model_info$outcomes[i]
    pred <- predicted[, i]
    actual <- actual_values[[outcome]]
    
    # Remove NA pairs
    valid_idx <- !is.na(pred) & !is.na(actual)
    pred_valid <- pred[valid_idx]
    actual_valid <- actual[valid_idx]
    
    # Basic metrics
    RMSE <- sqrt(mean((pred_valid - actual_valid)^2))
    MAE <- mean(abs(pred_valid - actual_valid))
    
    # R² - calculate from OOB predictions for rfsrc
    R2 <- if (!is.null(model_info$model$regrOutput)) {
      outcome_regr <- model_info$model$regrOutput[[outcome]]
      
      if (!is.null(outcome_regr$predicted.oob)) {
        pred_oob <- outcome_regr$predicted.oob
        actual_train <- model_info$model$yvar[, outcome]
        valid_oob <- !is.na(pred_oob) & !is.na(actual_train)
        
        ss_res <- sum((actual_train[valid_oob] - pred_oob[valid_oob])^2)
        ss_tot <- sum((actual_train[valid_oob] - mean(actual_train[valid_oob]))^2)
        1 - (ss_res / ss_tot)
      } else {
        NA
      }
    } else {
      NA
    }
    
    # C-index
    c_index_result <- tryCatch({
      rcorr.cens(pred_valid, actual_valid)
    }, error = function(e) {
      c("C Index" = NA, Dxy = NA, S.D. = NA, n = NA)
    })
    
    c_index <- c_index_result["C Index"]
    c_index_se <- c_index_result["S.D."]
    
    # Somers' D
    somers_result <- tryCatch({
      somers2(pred_valid, actual_valid)
    }, error = function(e) {
      c(C = NA, Dxy = NA, S.D. = NA, n = NA)
    })
    
    dxy <- somers_result["Dxy"]
    Bias <- mean(pred_valid - actual_valid)
    N <- length(pred_valid)
    
    data.frame(
      outcome = outcome,
      N = N,
      RMSE = RMSE,
      MAE = MAE,
      R2 = R2,
      C_index = as.numeric(c_index),
      C_index_SE = as.numeric(c_index_se),
      Dxy = as.numeric(dxy),
      Bias = Bias,
      stringsAsFactors = FALSE
    )
  })
  
  combined_metrics <- do.call(rbind, outcome_metrics)
  combined_metrics$model_id <- model_id
  combined_metrics$variant <- model_info$variant
  
  return(combined_metrics)
})

metrics_dt <- rbindlist(metrics)


# Reorder columns
metrics_dt <- metrics_dt[, .(model_id, outcome, variant, N, 
                             RMSE, MAE, R2, C_index, C_index_SE)]

print(metrics_dt)

# Format for publication with confidence intervals
metrics_dt_pub <- metrics_dt[, .(
  model_id, outcome, variant, N,
  RMSE = round(RMSE, 3),
  MAE = round(MAE, 3),
  R2 = round(R2, 3),
  C_index = sprintf("%.3f ± %.3f", C_index, C_index_SE)
)]

print(metrics_dt_pub)

write.csv(metrics_dt_pub, file.path(output, "rfsrc_model_WSS_metrics_2.csv"), row.names = FALSE)

# Run correlations
correlations <- data.frame()

for (model_id in names(predictions)) {
  pred_matrix <- predictions[[model_id]]
  model_info <- final_models[[model_id]]
  
  for (i in seq_along(model_info$outcomes)) {
    outcome_i <- model_info$outcomes[i]  # local name -- avoid clobbering the global `outcome` vector
    pred <- pred_matrix[, i]
    actual <- actual_values[[outcome_i]]

    cor_test <- cor.test(pred, actual, method = "pearson", use = "complete.obs")

    correlations <- rbind(correlations, data.frame(
      model_id = model_id,
      outcome = outcome_i,
      variant = model_info$variant,
      Pearson_Correlation = cor_test$estimate,
      Correlation_P_Value = cor_test$p.value,
      CI_Lower = cor_test$conf.int[1],
      CI_Upper = cor_test$conf.int[2],
      Sample_Size = sum(!is.na(pred) & !is.na(actual))
    ))
  }
}

write.csv(correlations, file.path(output, "rfsrc_model_WSS_correlations.csv"), row.names = FALSE)

# Macro metrics
macro_metrics_dt <- metrics_dt[, .(
  n_outcomes = .N,
  
  RMSE_mean = mean(RMSE, na.rm = TRUE),
  RMSE_median = median(RMSE, na.rm = TRUE),
  
  MAE_mean = mean(MAE, na.rm = TRUE),
  MAE_median = median(MAE, na.rm = TRUE),
  
  R2_mean = mean(R2, na.rm = TRUE),
  R2_median = median(R2, na.rm = TRUE),
  
  C_index_mean = mean(C_index, na.rm = TRUE),
  C_index_median = median(C_index, na.rm = TRUE)
  
), by = .(model_id, variant)]

# Save macro metrics
write.csv(macro_metrics_dt, file.path(output, "rfsrc_macro_metrics_2.csv"), row.names = FALSE)

cat("\n=== All outputs saved to:", output, "===\n")
cat("Files created:\n")
cat("  - Individual prediction CSVs\n")
cat("  - rfsrc_model_metrics.csv\n")
cat("  - rfsrc_model_correlations.csv\n")
cat("  - rfsrc_macro_metrics.csv\n")