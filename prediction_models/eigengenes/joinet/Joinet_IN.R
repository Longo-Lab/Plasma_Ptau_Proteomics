# library(data.table)
# 
# # Your existing setup code
# variants <- list(
#   ALL = character(0),
#   NO_hybridNormControlScale = "HybControlNormScale",
#   No_APOE = c("APOE.e34", "APOE.e44"),
#   NO_PCS = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"),
#   NO_PCS_APOE = c("APOE.e34", "APOE.e44", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"),
#   Biomarkers = c("HybControlNormControlScale","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8", "APOE.e34", "APOE.e44")
# )
# 
# # CV settings
# K_outer <- 5
# K_inner <- 10
# alpha_grid <- seq(0, 1, 0.1)
# 
# final_models <- list()
# 
# for (tag in names(variants)) {
#   message("\n==== Variant: ", tag, " ====")
#   
#   # 1) Build data for this variant
#   x_cols <- get_xcols(tag)
#   if (length(x_cols) == 0) {
#     warning(sprintf("[%s] No predictors selected; skipping", tag))
#     next
#   }
#   
#   df <- as.data.frame(pheno2)[, c(outcome, x_cols), drop = FALSE]
#   df <- df[complete.cases(df), ]
#   
#   if (nrow(df) < max(K_outer, 10)) {
#     warning(sprintf("[%s] Too few rows after NA filtering (%d). Skipping.", tag, nrow(df)))
#     next
#   }
#   
#   Y <- as.matrix(df[, outcome, drop = FALSE])
#   X <- as.matrix(df[, x_cols, drop = FALSE])
#   
#   # 2) Find best alphas
#   set.seed(42)
#   grid <- CJ(alpha.base = alpha_grid, alpha.meta = alpha_grid)
#   best_mse <- Inf
#   best_ab <- best_am <- NA
#   
#   for (i in 1:nrow(grid)) {
#     ab <- grid$alpha.base[i]
#     am <- grid$alpha.meta[i]
#     mcv <- cv.joinet(
#       Y = Y, X = X, family = "gaussian",
#       type.measure = "mse",
#       alpha.base = ab, alpha.meta = am,
#       nfolds.ext = K_outer, nfolds.int = K_inner
#     )
#     mse <- mean(mcv["meta", ])
#     if (mse < best_mse) {
#       best_mse <- mse
#       best_ab <- ab
#       best_am <- am
#     }
#   }
#   
#   message(sprintf("[%s] Best alphas: alpha.base=%.2f  alpha.meta=%.2f", tag, best_ab, best_am))
#   
#   # 3) Create and save final model for this variant
#   model_id <- paste(paste(outcome, collapse="_"), tag, sep = "_")
#   
#   final_models[[model_id]] <- list(
#     model = joinet(Y = Y, X = X, family = "gaussian",
#                    alpha.base = best_ab, alpha.meta = best_am),
#     outcomes = outcome,
#     variant = tag,
#     alpha_base = best_ab,
#     alpha_meta = best_am,
#     x_cols = x_cols
#   )
#   
#   # Print progress
#   cat("Processed model:", model_id, "\n")
# }
# 
# # Save all final models
# saveRDS(final_models, file.path(output, "final_models.rds"))
# 
# # Create a summary table
# summary_table <- data.table(
#   Model_ID = names(final_models),
#   Variant = sapply(final_models, function(m) m$variant),
#   Alpha_Base = sapply(final_models, function(m) m$alpha_base),
#   Alpha_Meta = sapply(final_models, function(m) m$alpha_meta),
#   N_Features = sapply(final_models, function(m) length(m$x_cols))
# )
# 
# # Save summary table
# fwrite(summary_table, file.path(output, "final_models_summary.csv"))
# 
# print(summary_table)
# )


# 
# library(data.table)
# library(foreach)
# library(doParallel)
# 
# # Your existing setup code
# variants <- list(
#   ALL = character(0),
#   NO_hybridNormControlScale = "HybControlNormScale",
#   No_APOE = c("APOE.e34", "APOE.e44"),
#   NO_PCS = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"),
#   NO_PCS_APOE = c("APOE.e34", "APOE.e44", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"),
#   Biomarkers = c("HybControlNormScale","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8", "APOE.e34", "APOE.e44", "age,)
# )
# 
# # CV settings
# K_outer <- 5
# K_inner <- 10
# alpha_grid <- seq(0, 1, 0.2)  # Coarser grid for faster computation
# 
# # Set up parallel backend
# cores <- detectCores() - 1  # Use all but one core
# cl <- makeCluster(cores)
# registerDoParallel(cl)
# 
# # Export necessary functions and data to the parallel workers
# clusterExport(cl, c("cv.joinet", "joinet", "get_xcols"))  # Add any other necessary functions
# 
# final_models <- list()
# 
# for (tag in names(variants)) {
#   message("\n==== Variant: ", tag, " ====")
#   
#   # 1) Build data for this variant
#   x_cols <- get_xcols(tag)
#   if (length(x_cols) == 0) {
#     warning(sprintf("[%s] No predictors selected; skipping", tag))
#     next
#   }
#   
#   df <- as.data.frame(pheno2)[, c(outcome, x_cols), drop = FALSE]
#   df <- df[complete.cases(df), ]
#   
#   if (nrow(df) < max(K_outer, 10)) {
#     warning(sprintf("[%s] Too few rows after NA filtering (%d). Skipping.", tag, nrow(df)))
#     next
#   }
#   
#   Y <- as.matrix(df[, outcome, drop = FALSE])
#   X <- as.matrix(df[, x_cols, drop = FALSE])
#   
#   # 2) Find best alphas in parallel
#   set.seed(42)
#   grid <- CJ(alpha.base = alpha_grid, alpha.meta = alpha_grid)
#   
#   results <- foreach(i = 1:nrow(grid), .combine = rbind, .packages = c("data.table")) %dopar% {
#     ab <- grid$alpha.base[i]
#     am <- grid$alpha.meta[i]
#     mcv <- cv.joinet(
#       Y = Y, X = X, family = "gaussian",
#       type.measure = "mse",
#       alpha.base = ab, alpha.meta = am,
#       nfolds.ext = K_outer, nfolds.int = K_inner
#     )
#     mse <- mean(mcv["meta", ])
#     data.table(alpha.base = ab, alpha.meta = am, mse = mse)
#   }
#   
#   best_result <- results[which.min(mse)]
#   best_ab <- best_result$alpha.base
#   best_am <- best_result$alpha.meta
#   
#   message(sprintf("[%s] Best alphas: alpha.base=%.2f  alpha.meta=%.2f", tag, best_ab, best_am))
#   
#   # 3) Create and save final model for this variant
#   model_id <- paste(paste(outcome, collapse="_"), tag, sep = "_")
#   
#   final_models[[model_id]] <- list(
#     model = joinet(Y = Y, X = X, family = "gaussian",
#                    alpha.base = best_ab, alpha.meta = best_am),
#     outcomes = outcome,
#     variant = tag,
#     alpha_base = best_ab,
#     alpha_meta = best_am,
#     x_cols = x_cols
#   )
#   
#   # Print progress
#   cat("Processed model:", model_id, "\n")
# }
# 
# # Stop parallel backend
# stopCluster(cl)
# 
# # Save all final models
# saveRDS(final_models, file.path(output, "final_models.rds"))
# 
# # Create a summary table
# summary_table <- data.table(
#   Model_ID = names(final_models),
#   Variant = sapply(final_models, function(m) m$variant),
#   Alpha_Base = sapply(final_models, function(m) m$alpha_base),
#   Alpha_Meta = sapply(final_models, function(m) m$alpha_meta),
#   N_Features = sapply(final_models, function(m) length(m$x_cols))
# )
# 
# # Save summary table
# fwrite(summary_table, file.path(output, "final_models_summary.csv"))
# 
# print(summary_table)

### UPDATE THIS: 

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


variants <- list(
  Biomarkers_Only = c("HybControlNormScale","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8", "APOE.e34", "APOE.e44", "Age", "Gender"),
  Biomarkers_Age_Sex = c("HybControlNormScale","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8", "APOE.e34", "APOE.e44"),
  Age_Sex_HybridNormControl = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","APOE.e34", "APOE.e44"),
  Age_Sex_HybridNormControl_PCs = c("APOE.e34", "APOE.e44"),
  ALL = character(0)
)

# CV settings
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
saveRDS(final_models, file.path(output, "final_models_z.rds"))
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

fwrite(Age_Sex_HybridNormControl, file.path(output, "Age_Sex_HybridNormControl_z_predictions.csv"))

Age_Sex_HybridNormControl_PCs <- as.data.frame(predictions[[4]])
Age_Sex_HybridNormControl_PCs <- cbind(SampleID = rownames(actual_values), Age_Sex_HybridNormControl_PCs)

fwrite(Age_Sex_HybridNormControl_PCs, file.path(output, "Age_Sex_HybridNormControl_PCs_z_predictions.csv"))

ALL <- as.data.frame(predictions[[5]])
ALL<- cbind(SampleID = rownames(actual_values), ALL)

fwrite(ALL, file.path(output, "ALL_z_predictions.csv"))

# Biomarkers <- as.data.frame(predictions[[6]])
# Biomarkers<- cbind(SampleID = rownames(actual_values), Biomarkers)
# 
# fwrite(Biomarkers, file.path(output, "Biomarkers_eigengenes_z_predictions.csv"))


# metrics <- lapply(names(predictions), function(model_id) {
#   predicted <- predictions[[model_id]]
#   model_info <- final_models[[model_id]]
#   
#   # Calculate metrics for each outcome
#   outcome_metrics <- lapply(seq_along(model_info$outcomes), function(i) {
#     outcome <- model_info$outcomes[i]
#     pred <- predicted[, i]
#     actual <- actual_values[[outcome]]
#     
#     RMSE <- sqrt(mean((pred - actual)^2, na.rm = TRUE))
#     MAE <- mae(actual, pred)
#     R2 <- cor(actual, pred, use = "complete.obs")^2
#     AUC_RROC <- auc(actual, pred)
#     result <- somers2(pred, actual)
#     c_index <- result["C"]
#     
#     list(outcome = outcome, RMSE = RMSE, MAE = MAE, R2 = R2, AUC_RROC = AUC_RROC, c_index = c_index)
#   })
#   
#   # Combine metrics for all outcomes
#   combined_metrics <- do.call(rbind, lapply(outcome_metrics, as.data.frame))
#   combined_metrics$model_id <- model_id
#   combined_metrics$variant <- model_info$variant
#   
#   return(combined_metrics)
# })
# 
# # Combine all metrics into a single data.table
# metrics_dt <- rbindlist(metrics)
# 
# # Reorder columns to match your desired output
# metrics_dt <- metrics_dt[, .(model_id, outcome, variant, RMSE, MAE, R2, AUC_RROC, c_index)]
# 
# # Print metrics
# print(metrics_dt)
# 
# write.csv(metrics_dt, file.path(output, "joinet_model_metrics_eigenegenes.csv"), row.names = FALSE)
library(Hmisc)
library(data.table)

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


# library(Metrics)
# 
# # calculate_auc_rroc <- function(actual, predicted) {
# #   errors <- predicted - actual
# #   n <- length(errors)
# #   ranked_errors <- rank(errors)
# #   auc_rroc <- sum(ranked_errors * errors) / (n * sum(errors))
# #   return(auc_rroc)
# # }
# 
# predictions <- vector("list", length(final_models))
# names(predictions) <- names(final_models)
# 
# # Loop through models and make predictions
# for (model_id in names(final_models)) {
#   model_info <- final_models[[model_id]]
#   
#   # Prepare test data for this model
#   X_test <- test_data[, model_info$x_cols, drop = FALSE]
#   
#   # Make predictions
#   pred <- predict(model_info$model, newx = X_test, s = model_info$lambda)
#   
#   # Store predictions
#   predictions[[model_id]] <- pred[,1]  # pred is a matrix, we take the first column
#   
#   # Print progress
#   cat("Processed predictions for model:", model_id, "\n")
# }
# 
# # Convert predictions list to a data.table
# predictions_dt <- as.data.table(predictions)
# 
# # Calculate metrics for each model
# 
# metrics <- lapply(names(predictions), function(model_id) {
#   predicted <- predictions[[model_id]]
#   actual <- actual_values[[final_models[[model_id]]$outcome]]
#   
#   RMSE <- sqrt(mean((predicted - actual)^2, na.rm = TRUE))
#   MAE <- mae(actual, predicted)
#   R2 <- cor(actual, predicted, use = "complete.obs")^2
#   AUC_RROC <- auc(actual, predicted)
#   
#   list(RMSE = RMSE, MAE = MAE, R2 = R2, AUC_RROC = AUC_RROC)
# })
# 
# # Create metrics data.table
# metrics_dt <- data.table(
#   model_id = names(predictions),
#   outcome = sapply(names(predictions), function(id) final_models[[id]]$outcome),
#   variant = sapply(names(predictions), function(id) final_models[[id]]$variant),
#   RMSE = sapply(metrics, `[[`, "RMSE"),
#   MAE = sapply(metrics, `[[`, "MAE"),
#   R2 = sapply(metrics, `[[`, "R2"),
#   AUC_RROC = sapply(metrics, `[[`, "AUC_RROC")
# )
# 
# print(metrics_dt)
# 
# 
# library(data.table)
# library(Metrics)  # for mae function
# library(pROC)     # for auc function
# 
# predictions <- vector("list", length(final_models))
# names(predictions) <- names(final_models)
# 
# # Loop through models and make predictions
# for (model_id in names(final_models)) {
#   model_info <- final_models[[model_id]]
#   
#   # Check if x_cols exist in test_data
#   if (!all(model_info$x_cols %in% names(test_data))) {
#     missing_cols <- setdiff(model_info$x_cols, names(test_data))
#     cat("Missing columns in test data for model", model_id, ":", paste(missing_cols, collapse = ", "), "\n")
#     next
#   }
#   
#   # Prepare test data for this model
#   X_test <- as.matrix(test_data[, model_info$x_cols, drop = FALSE])
#   
#   # Make predictions
#   tryCatch({
#     if ("joinet" %in% class(model_info$model)) {
#       pred <- predict(model_info$model, newx = X_test)
#       if (is.list(pred) && !is.null(pred$meta)) {
#         pred <- pred$meta
#       }
#     } else {
#       pred <- predict(model_info$model, newx = X_test)
#     }
#     
#     # Store predictions
#     predictions[[model_id]] <- as.vector(pred)
#     
#     # Print progress
#     cat("Processed predictions for model:", model_id, "\n")
#   }, error = function(e) {
#     cat("Error in prediction for model", model_id, ":", e$message, "\n")
#   })
# }
# 
# # Convert predictions list to a data.table
# predictions_dt <- as.data.table(predictions)
# 
# # Calculate metrics for each model
# metrics <- lapply(names(predictions), function(model_id) {
#   predicted <- predictions[[model_id]]
#   actual <- actual_values[[final_models[[model_id]]$outcomes[1]]]  # Assuming single outcome
#   
#   if (length(predicted) != length(actual)) {
#     cat("Mismatch in lengths for model", model_id, "\n")
#     return(list(RMSE = NA, MAE = NA, R2 = NA, AUC_RROC = NA))
#   }
#   
#   RMSE <- sqrt(mean((predicted - actual)^2, na.rm = TRUE))
#   MAE <- mae(actual, predicted)
#   R2 <- cor(actual, predicted, use = "complete.obs")^2
#   AUC_RROC <- auc(actual, predicted)
#   
#   list(RMSE = RMSE, MAE = MAE, R2 = R2, AUC_RROC = AUC_RROC)
# })
# 
# # Create metrics data.table
# metrics_dt <- data.table(
#   model_id = names(predictions),
#   outcome = sapply(names(predictions), function(id) final_models[[id]]$outcomes[1]),  # Assuming single outcome
#   variant = sapply(names(predictions), function(id) final_models[[id]]$variant),
#   RMSE = sapply(metrics, `[[`, "RMSE"),
#   MAE = sapply(metrics, `[[`, "MAE"),
#   R2 = sapply(metrics, `[[`, "R2"),
#   AUC_RROC = sapply(metrics, `[[`, "AUC_RROC")
# )
# 
# print(metrics_dt)
