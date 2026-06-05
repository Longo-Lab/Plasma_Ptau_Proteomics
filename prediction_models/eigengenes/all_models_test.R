#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 02-linreg_glm.R
##
## Version: 0.0.3
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
library(parallel)
library(doParallel)



# Inputs ----------------
setwd("C:/Users/audweber/Box/plasma_ptau_proteomics")
inpath <- "C:/Users/audweber/Box/plasma_ptau_proteomics" # don't repeat things you don't have to
input <- file.path(inpath, "dataProt_SS-205063.merged.QC.csv")
prots <- file.path(inpath, "prot_names.txt") 
pcs <- file.path(inpath, "genomic/mypcair.PCs.csv")
prediction_scores <- file.path(inpath, # wrap at 80 characters
  "celltype_scores/dataProt_SS-205063.STRING.INT.eigengenes.csv")
output <- 'prediction_models/eigengenes/glmnet/raw'
covs <- 'Age,Gender,APOE,HybControlNormScale'   
output <- 'prediction_models/eigengenes/joinet/INT'
covs <- 'Age,Gender,APOE,HybControlNormScale'   
outcome <- make.names(setdiff(names(fread(prediction_scores, nrows = 0)), c("SampleID", "Neddylation")), unique = TRUE)
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

# Add cell scores
if (prediction_scores != "") {
  print(sprintf("Reading Prediction scores for inclusion from %s", prediction_scores))
  x <- fread(prediction_scores, colClasses = c(SampleID = "character"))
  a <- merge(a, x, by = "SampleID")
}

# filter outliars 
drop = a[abs(ConnectivityZscore) > 10, SampleID]
a <- a[!SampleID %in% drop][, !"Neddylation", with = FALSE]

#Only keep baseline visits where plasma is avaliable 
a = a[order(DS1V)][, head(.SD, 1), by = "PIDN"]

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
pheno <- `rownames<-`(setDF(copy(a[complete.cases(a[, ..cols]), 
                      setdiff(cols, "SampleID"), with = FALSE])), 
                      make.unique(as.character(a[complete.cases(a[, ..cols]), SampleID])))

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
# # Save summary table
# fwrite(summary_table, file.path(output, "final_models_summary.csv"))
# 
print(summary_table)


## setting up test data--------------------------------------------------------------------

a <- fread(input)

# Map in PCs assumes your primary data table exists 
if (pcs != "") {
  print(sprintf("Reading PCs for inclusion from %s", pcs))
  x <- fread(pcs, colClasses = c(PIDN = "character"))
  a <- merge(a, x, by = "PIDN")
  covs <- c(covs, names(a)[names(a) %like% "^PC\\d{1,2}$"]) # flexible for PCs
}

# Add cell enrichments

if (prediction_scores != "") {
  print(sprintf("Reading Prediction scores for inclusion from %s", prediction_scores))
  x <- fread(prediction_scores, colClasses = c(SampleID = "character"))
  a <- merge(a, x, by = "SampleID")
}


drop = a[abs(ConnectivityZscore) > 10, SampleID]
a = a[!SampleID %in% drop]
a[, Gender := as.numeric(Gender == "M")]
recode_dummies(a, "APOE", keep = c("e34", "e44"))
covs <- unique(c(covs[covs != "APOE"], "APOE.e34", "APOE.e44"))


test_data <- a

ids <- test_data [order(DS1V)][, head(.SD, 1), by = "PIDN"][, SampleID]
f_ups <- test_data[!SampleID %in% ids]
columns_to_check <- c(  "Age","Gender","HybControlNormScale", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","APOE.e34","APOE.e44",
                        "PlasmaPTau181","PlasmaAB142P","PlasmaAB140P","PlasmaABRatio","PlasmapTau217","PlasmapTau217_AB42Ratio",
                        "PlasmaGFAP","PlasmaNfL")
drops <- f_ups[!complete.cases(f_ups[, ..columns_to_check ]), SampleID]

test_data <- f_ups[!SampleID %in% drops]

names(test_data) <- names(name_map <- setNames(names(test_data), make.names(names(test_data), unique = TRUE)))

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
    R2 <- cor(actual_valid, pred_valid)^2
    
    # C-index (like AUC for continuous outcomes)
    c_index_result <- tryCatch({
      rcorr.cens(pred_valid, actual_valid)
    }, error = function(e) {
      c(C = NA, Dxy = NA, S.D. = NA, n = NA)
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
      R2 = R2,
      C_index = c_index,
      C_index_SE = c_index_se,
      Dxy = dxy,
      Bias = Bias
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
                             RMSE, MAE, R2, C_index, C_index_SE, Dxy, Bias)]

print(metrics_dt)

# Format for publication with confidence intervals
metrics_dt_pub <- metrics_dt[, .(
  model_id, outcome, variant, N,
  RMSE = round(RMSE, 3),
  MAE = round(MAE, 3),
  R2 = round(R2, 3),
  C_index = sprintf("%.3f ± %.3f", C_index, C_index_SE),
  Dxy = round(Dxy, 3),
  Bias = round(Bias, 3)
)]

print(metrics_dt_pub)

write.csv(metrics_dt_pub, file.path(output, "joinet_model_metrics_c_index_eigengenes_z_score.csv"), row.names = FALSE)

#run correlations-------------------------------------------------------------------------------------------------- 
correlations <- data.frame()

for (model_id in names(predictions)) {
  pred_matrix <- predictions[[model_id]]
  model_info <- final_models[[model_id]]
  
  for (i in seq_along(model_info$outcomes)) {
    outcome <- model_info$outcomes[i]
    pred <- pred_matrix[, i]
    actual <- actual_values[[outcome]]
    
    cor_test <- cor.test(pred, actual, method = "pearson", use = "complete.obs")
    
    correlations <- rbind(correlations, data.frame(
      model_id = model_id,
      outcome = outcome,
      variant = model_info$variant,
      Pearson_Correlation = cor_test$estimate,
      Correlation_P_Value = cor_test$p.value,
      CI_Lower = cor_test$conf.int[1],
      CI_Upper = cor_test$conf.int[2],
      Sample_Size = length(pred)
    ))
  }
}

write.csv(correlations, file.path(output, "joinet_model_correlations_eigengenes.csv"), row.names = FALSE)

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
dir.create(output, showWarnings = FALSE, recursive = TRUE)
write.csv(macro_metrics_dt, file.path(output, "joinet_macro_metrics.csv"), row.names = FALSE)

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

write.csv(macro_cor_dt, file.path(output, "joinet_macro_correlations.csv"), row.names = FALSE)


##GLMNET----------------------------------------------------------------------------

output <- "prediction_models/eigengenes/glmnet/INT"

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

som <- Hmisc::somers2(
  as.numeric(pred_valid),
  as.numeric(actual_valid)
)

dxy <- as.numeric(som$Dxy)
c_index <- (dxy + 1) / 2

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
  R2   <- cor(actual_valid, pred_valid)^2
  
  c_index_result <- tryCatch({
    rcorr.cens(pred_valid, actual_valid)
  }, error = function(e) {
    c(C = NA, Dxy = NA, S.D. = NA, n = NA)
  })
  
  c_index <- c_index_result["C Index"]
  c_index_se <- c_index_result["S.D."]
  
  somers_result <- tryCatch({
    somers2(pred_valid, actual_valid)
  }, error = function(e) {
    c(C = NA, Dxy = NA, S.D. = NA, n = NA)
  })
  
  dxy <- somers_result["Dxy"]
  Bias <- mean(pred_valid - actual_valid)
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
    Bias = Bias,
    row.names = NULL
  )
})

metrics_dt <- rbindlist(metrics, fill = TRUE)
metrics_dt <- metrics_dt[, .(model_id, outcome, variant, N,
                             RMSE, MAE, R2, C_index, C_index_SE, Dxy, Bias)]
print(metrics_dt)

metrics_dt_pub <- metrics_dt[, .(
  model_id, outcome, variant, N,
  RMSE = round(RMSE, 3),
  MAE = round(MAE, 3),
  R2 = round(R2, 3),
  C_index = sprintf("%.3f ± %.3f", C_index, C_index_SE),
  Dxy = round(Dxy, 3),
  Bias = round(Bias, 3)
)]

write.csv(metrics_dt_pub,
          file.path(output, "glmnet_model_metrics_c_index_eigengenes_INT.csv"),
          row.names = FALSE)

# ============================================================
# CORRELATIONS (your style)
# ============================================================
correlations <- data.frame()

for (model_id in names(predictions)) {
  pred <- predictions[[model_id]]
  model_info <- final_models[[model_id]]
  outcome <- model_info$outcome
  actual <- actual_values[[outcome]]
  
  cor_test <- cor.test(pred, actual, method = "pearson", use = "complete.obs")
  
  correlations <- rbind(correlations, data.frame(
    model_id = model_id,
    outcome = outcome,
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
