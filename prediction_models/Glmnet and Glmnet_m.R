library(parallel)
library(doParallel)
library(foreach)
library(data.table)
library(glmnet)

output <- 'prediction_models/eigengenes/glmnet/INT'

alpha_grid <- seq(0, 1, by = 0.2)
K_outer <- 5
lambda_choice <- "lambda.min"  # or "lambda.1se"

# -------------------------
# VARIANTS
# -------------------------
variants <- list(
  ALL = character(0),
  NO_hybridNormControlScale = "HybControlNormScale",
  No_APOE = c("APOE.e34", "APOE.e44"),
  NO_PCS = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"),
  NO_PCS_APOE = c("APOE.e34", "APOE.e44", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"),
  Biomarkers = c("HybControlNormControlScale","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8",
                 "APOE.e34", "APOE.e44", "Age", "Gender")
)

# -------------------------
# HELPERS
# -------------------------
prep_Y <- function(df, outcome_cols) {
  Y_df <- df[, outcome_cols, drop = FALSE]
  Y_df <- as.data.frame(lapply(Y_df, function(v) {
    if (is.factor(v)) v <- as.character(v)
    as.numeric(v)
  }))
  good <- sapply(Y_df, function(v) {
    v <- v[is.finite(v)]
    length(v) > 5 && sd(v) > 0
  })
  Y_df <- Y_df[, good, drop = FALSE]
  list(Y = as.matrix(Y_df), outcomes = colnames(Y_df))
}

prep_X <- function(df, x_cols_raw) {
  X_df <- df[, x_cols_raw, drop = FALSE]
  X <- model.matrix(~ . - 1, data = X_df)  # dummy encode, no intercept
  list(X = X, x_cols_mm = colnames(X))
}

# Robust scorer: works if cvm is matrix (multi-outcome) or vector (single outcome)
score_cvfit <- function(cvfit) {
  cvm <- cvfit$cvm
  
  # cvm can be:
  # - matrix: lambdas x outcomes  (multi-response)
  # - numeric vector: length = #lambdas (single response)
  if (is.null(cvm)) return(NA_real_)
  
  if (is.matrix(cvm)) {
    rm <- rowMeans(cvm, na.rm = TRUE)
  } else {
    rm <- as.numeric(cvm)
  }
  
  rm <- rm[is.finite(rm)]
  if (length(rm) == 0) return(NA_real_)
  min(rm)
}

# -------------------------
# MAIN: CV + SAVE FINAL MODELS
# -------------------------
final_models <- list()
summary_rows <- list()

for (tag in names(variants)) {
  cat("\n==== Variant:", tag, "====\n")
  
  # Must exist in your environment
  x_cols_raw <- get_xcols(tag)
  cat(sprintf("[%s] raw predictors: %d\n", tag, length(x_cols_raw)))
  
  # Build df and drop NA rows across outcomes + predictors
  df <- as.data.frame(pheno2)[, c(outcome, x_cols_raw), drop = FALSE]
  df <- df[complete.cases(df), , drop = FALSE]
  cat(sprintf("[%s] rows after complete.cases: %d\n", tag, nrow(df)))
  
  if (nrow(df) < 10) {
    cat(sprintf("[%s] Too few rows; skipping.\n", tag))
    next
  }
  
  # Build Y and X robustly
  yobj <- prep_Y(df, outcome)
  Y <- yobj$Y
  kept_outcomes <- yobj$outcomes
  if (ncol(Y) < 1) {
    cat(sprintf("[%s] No valid outcomes after numeric/variance filtering; skipping.\n", tag))
    next
  }
  
  xobj <- prep_X(df, x_cols_raw)
  X <- xobj$X
  x_cols_mm <- xobj$x_cols_mm
  
  # Guard: keep only finite rows
  ok <- complete.cases(X) & complete.cases(Y)
  X <- X[ok, , drop = FALSE]
  Y <- Y[ok, , drop = FALSE]
  
  cat(sprintf("[%s] Final matrices: n=%d, p=%d, q=%d\n", tag, nrow(X), ncol(X), ncol(Y)))
  
  K_use <- min(K_outer, nrow(X))
  if (K_use < 3) {
    cat(sprintf("[%s] Not enough rows for CV folds; skipping.\n", tag))
    next
  }
  
  # ---- Alpha grid search ----
  scores <- data.table(alpha = alpha_grid, mse = NA_real_)
  
  for (i in seq_along(alpha_grid)) {
    a <- alpha_grid[i]
    cvfit <- cv.glmnet(
      x = X, y = Y,
      family = "mgaussian",
      alpha = a,
      nfolds = K_use,
      type.measure = "mse",
      standardize = TRUE
    )
    scores$mse[i] <- score_cvfit(cvfit)
    cat(sprintf("[%s] alpha=%.2f  mse=%.6f\n", tag, a, scores$mse[i]))
  }
  
  scores <- scores[is.finite(mse)]
  if (nrow(scores) == 0) {
    cat(sprintf("[%s] All alpha scores NA; skipping.\n", tag))
    next
  }
  
  best_alpha <- scores$alpha[which.min(scores$mse)]
  best_mse   <- min(scores$mse)
  cat(sprintf("[%s] BEST alpha=%.2f (mse=%.6f)\n", tag, best_alpha, best_mse))
  
  # ---- Fit final CV model on ALL data at best alpha ----
  final_cv <- cv.glmnet(
    x = X, y = Y,
    family = "mgaussian",
    alpha = best_alpha,
    nfolds = K_use,
    type.measure = "mse",
    standardize = TRUE
  )
  
  best_lambda <- if (lambda_choice == "lambda.1se") final_cv$lambda.1se else final_cv$lambda.min
  
  model_id <- paste("glmnet_mgaussian", tag, sep = "_")
  
  final_models[[model_id]] <- list(
    model      = final_cv,
    variant    = tag,
    outcomes   = kept_outcomes,
    alpha      = best_alpha,
    lambda     = best_lambda,
    x_cols_raw = x_cols_raw,
    x_cols_mm  = x_cols_mm,
    n          = nrow(X),
    p          = ncol(X),
    q          = ncol(Y)
  )
  
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    model_id     = model_id,
    variant      = tag,
    n            = nrow(X),
    p            = ncol(X),
    q            = ncol(Y),
    best_alpha   = best_alpha,
    best_lambda  = best_lambda,
    best_mse     = best_mse
  )
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
print(summary_dt)

saveRDS(final_models, file.path(output, "final_models_glmnet_mgaussian_INT.rds"))
fwrite(summary_dt, file.path(output, "glmnet_mgaussian_summary_by_variant_INT.csv"))


## predictions---------------------------------------------------------------------------


metrics <- lapply(names(predictions), function(model_id) {
  predicted <- predictions[[model_id]]
  model_info <- final_models[[model_id]]
  
  outcome_metrics <- lapply(seq_along(model_info$outcomes), function(i) {
    outcome <- model_info$outcomes[i]
    pred <- predicted[, i]
    
    # IMPORTANT: keep this EXACTLY as your code intended
    # Use actual_values[[outcome]] like you wrote (not [, outcome])
    actual <- actual_values[[outcome]]
    
    # Remove NA pairs
    valid_idx <- !is.na(pred) & !is.na(actual)
    pred_valid <- pred[valid_idx]
    actual_valid <- actual[valid_idx]
    
    # Basic metrics
    RMSE <- sqrt(mean((pred_valid - actual_valid)^2))
    MAE  <- mean(abs(pred_valid - actual_valid))
    R2   <- cor(actual_valid, pred_valid)^2
    
    # C-index (like AUC for continuous outcomes)
    c_index_result <- tryCatch({
      rcorr.cens(pred_valid, actual_valid)
    }, error = function(e) {
      c(C = NA, Dxy = NA, S.D. = NA, n = NA)
    })
    
    # EXACT extraction you used
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

metrics_dt <- rbindlist(metrics, fill = TRUE)

# Reorder columns (exactly your layout)
metrics_dt <- metrics_dt[, .(model_id, outcome, variant, N,
                             RMSE, MAE, R2, C_index, C_index_SE, Dxy, Bias)]

print(metrics_dt)

# Publication formatting with confidence intervals (exactly your format)
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

write.csv(
  metrics_dt_pub,
  file.path(output, "glmnet_mgaussian_model_metrics_c_index_eigengenes.csv"),
  row.names = FALSE
)

# -------------------------
# Correlations (your exact code structure, just saved under glmnet name)
# -------------------------
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

write.csv(
  correlations,
  file.path(output, "glmnet_mgaussian_model_correlations_eigengenes.csv"),
  row.names = FALSE
)

##GLMNET----------------------------------------------------------------------------

output <- "prediction_models/eigengenes/glmnet/INT"

alpha_grid <- seq(0, 1, by = 0.2)
K_outer <- 5
lambda_choice <- "lambda.min"  # or "lambda.1se"

# -------------------------
# VARIANTS (your list)
# -------------------------
variants <- list(
  ALL = character(0),
  NO_hybridNormControlScale = "HybControlNormScale",
  No_APOE = c("APOE.e34", "APOE.e44"),
  NO_PCS = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"),
  NO_PCS_APOE = c("APOE.e34", "APOE.e44", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"),
  Biomarkers = c("HybControlNormScale","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8",
                 "APOE.e34", "APOE.e44", "Age", "Gender")
)

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
          file.path(output, "glmnet_model_metrics_c_index_eigengenes_raw.csv"),
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
          file.path(output, "glmnet_model_correlations_eigengenes_raw.csv"),
          row.names = FALSE)



#GLmnet (mgaussian metrics) -----------------------------------------------------------

metrics <- list()
counter <- 1

for (model_id in names(predictions)) {
  pred_df <- predictions[[model_id]]  # This is a data frame with multiple outcome columns
  model_info <- final_models[[model_id]]
  
  # Loop through each outcome column in the predictions
  for (outcome_col in colnames(pred_df)) {
    pred <- pred_df[[outcome_col]]  # Extract predictions for this specific outcome
    
    # Get actual values for this outcome
    if (outcome_col %in% names(actual_values)) {
      actual <- actual_values[[outcome_col]]
    } else {
      cat("Warning: No actual values found for outcome:", outcome_col, "\n")
      next
    }
    
    # Ensure same length
    min_length <- min(length(pred), length(actual))
    pred <- pred[1:min_length]
    actual <- actual[1:min_length]
    
    # Remove missing values
    valid_idx <- !is.na(pred) & !is.na(actual)
    pred_valid <- pred[valid_idx]
    actual_valid <- actual[valid_idx]
    
    if (length(pred_valid) < 2) {
      cat("Warning: Not enough valid predictions for", model_id, outcome_col, "\n")
      next
    }
    
    # Calculate metrics
    RMSE <- sqrt(mean((pred_valid - actual_valid)^2))
    MAE  <- mean(abs(pred_valid - actual_valid))
    R2   <- cor(actual_valid, pred_valid)^2
    
    # C-index and Somers' D
    c_index_result <- tryCatch({
      Hmisc::rcorr.cens(pred_valid, actual_valid)
    }, error = function(e) {
      c("C Index" = NA, "Dxy" = NA, "S.D." = NA, "n" = NA)
    })
    
    somers_result <- tryCatch({
      Hmisc::somers2(pred_valid, actual_valid)
    }, error = function(e) {
      c("C" = NA, "Dxy" = NA, "S.D." = NA, "n" = NA)
    })
    
    # Extract values
    c_index <- as.numeric(c_index_result["C Index"])
    c_index_se <- as.numeric(c_index_result["S.D."])
    dxy <- as.numeric(somers_result["Dxy"])
    bias <- mean(pred_valid - actual_valid)
    n <- length(pred_valid)
    
    # Store metrics
    metrics[[counter]] <- data.frame(
      model_id = model_id,
      outcome = outcome_col,
      variant = model_info$variant,
      N = n,
      RMSE = RMSE,
      MAE = MAE,
      R2 = R2,
      C_index = c_index,
      C_index_SE = c_index_se,
      Dxy = dxy,
      Bias = bias,
      stringsAsFactors = FALSE
    )
    
    counter <- counter + 1
  }
}

# Combine all metrics
library(data.table)
metrics_dt <- rbindlist(metrics, fill = TRUE)

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
          file.path(output, "glmnet_model_metrics_c_index_eigengenes_gmaussian_raw.csv"),
          row.names = FALSE)
correlations <- data.frame()

for (model_id in names(predictions)) {
  pred_df <- predictions[[model_id]]  # This is a data frame with multiple outcome columns
  model_info <- final_models[[model_id]]
  
  # Loop through each outcome column in the predictions
  for (outcome_col in colnames(pred_df)) {
    pred <- pred_df[[outcome_col]]  # Extract predictions for this specific outcome
    
    # Get actual values for this outcome
    if (outcome_col %in% names(actual_values)) {
      actual <- actual_values[[outcome_col]]
    } else {
      cat("Warning: No actual values found for outcome:", outcome_col, "\n")
      next
    }
    
    # Ensure same length
    min_length <- min(length(pred), length(actual))
    pred <- pred[1:min_length]
    actual <- actual[1:min_length]
    
    # Check if we have enough valid pairs for correlation
    valid_pairs <- complete.cases(pred, actual)
    if (sum(valid_pairs) < 3) {
      cat("Warning: Not enough complete cases for correlation -", model_id, outcome_col, "\n")
      next
    }
    
    # Calculate correlation
    cor_test <- tryCatch({
      cor.test(pred, actual, method = "pearson", use = "complete.obs")
    }, error = function(e) {
      cat("Error in correlation for", model_id, outcome_col, ":", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(cor_test)) {
      correlations <- rbind(correlations, data.frame(
        model_id = model_id,
        outcome = outcome_col,
        variant = model_info$variant,
        Pearson_Correlation = cor_test$estimate,
        Correlation_P_Value = cor_test$p.value,
        CI_Lower = cor_test$conf.int[1],
        CI_Upper = cor_test$conf.int[2],
        Sample_Size = sum(valid_pairs),  # Number of complete cases used
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Clean up row names
rownames(correlations) <- NULL

# Print results
print(correlations)


write.csv(correlations,
          file.path(output, "glmnet_model_correlations_eigengenes_mgaussian_raw.csv"),
          row.names = FALSE)


