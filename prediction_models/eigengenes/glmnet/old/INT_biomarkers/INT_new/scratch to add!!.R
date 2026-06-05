

results[, model_id := paste(outcome, variant, sep = "_")]

final_models <- vector("list", nrow(results))
names(final_models) <- results$model_id

for (i in seq_len(nrow(results))) {
  row <- results[i]
  
  out_q <- row$outcome
  tag   <- row$variant
  a     <- row$best_alpha
  lam   <- row$best_lambda
  
  # predictors for THIS outcome's variant
  drops  <- variants[[tag]]
  x_cols <- setdiff(predictor_pool, intersect(drops, names(pheno2)))
  
  # training set for this outcome
  if (is.data.table(pheno2)) {
    df <- pheno2[, c(out_q, x_cols), with = FALSE]
    df <- df[!is.na(get(out_q))]
  } else {
    df <- pheno2[, c(out_q, x_cols)]
    df <- df[!is.na(df[[out_q]]), ]
  }
  
  y <- as.numeric(df[[out_q]])
  
  X <- as.matrix(df[, x_cols, drop = FALSE])
  
  final_models[[row$model_id]] <- list(
    model   = glmnet(X, y, alpha = a, lambda = c(lam), family = "gaussian"),
    outcome = out_q,
    variant = tag,
    alpha   = a,
    lambda  = lam,
    x_cols  = x_cols
  )
  
  # Print progress
  cat("Processed model", i, "of", nrow(results), ":", row$model_id, "\n")
}

# Assuming your test data is in a data frame called 'test_data'

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

## at this point go to INT values code 

test_data <- as.matrix(pheno2_norm[, !(names(pheno2_norm) %in% outcome), drop = FALSE])

actual_values <- test_data[, outcome, drop = FALSE]

## loop through to make predictions 

predictions <- vector("list", length(final_models))
names(predictions) <- names(final_models)

# Loop through models and make predictions
for (model_id in names(final_models)) {
  model_info <- final_models[[model_id]]
  
  # Prepare test data for this model
  X_test <- test_data[, model_info$x_cols, drop = FALSE]
  
  # Make predictions
  pred <- predict(model_info$model, newx = X_test, s = model_info$lambda)
  
  # Store predictions
  predictions[[model_id]] <- pred[,1]  # pred is a matrix, we take the first column
  
  # Print progress
  cat("Processed predictions for model:", model_id, "\n")
}

# Convert predictions list to a data.table
predictions_dt <- as.data.table(predictions)

# Calculate RMSE for each model
rmse_values <- sapply(names(predictions), function(model_id) {
  predicted <- predictions[[model_id]]
  actual <- actual_values[[final_models[[model_id]]$outcome]]
  sqrt(mean((predicted - actual)^2, na.rm = TRUE))
})

# Create RMSE data.table
rmse_dt <- data.table(
  model_id = names(rmse_values),
  outcome = sapply(names(rmse_values), function(id) final_models[[id]]$outcome),
  variant = sapply(names(rmse_values), function(id) final_models[[id]]$variant),
  rmse = rmse_values
  
)

# Sort by RMSE
setorder(rmse_dt, rmse)
print(rmse_dt)

# Find the best model for each outcome
best_models <- rmse_dt[, .SD[which.min(rmse)], by = outcome]
print(best_models)



library(Metrics)

# calculate_auc_rroc <- function(actual, predicted) {
#   errors <- predicted - actual
#   n <- length(errors)
#   ranked_errors <- rank(errors)
#   auc_rroc <- sum(ranked_errors * errors) / (n * sum(errors))
#   return(auc_rroc)
# }

predictions <- vector("list", length(final_models))
names(predictions) <- names(final_models)

# Loop through models and make predictions
for (model_id in names(final_models)) {
  model_info <- final_models[[model_id]]
  
  # Prepare test data for this model
  X_test <- test_data[, model_info$x_cols, drop = FALSE]
  
  # Make predictions
  pred <- predict(model_info$model, newx = X_test, s = model_info$lambda)
  
  # Store predictions
  predictions[[model_id]] <- pred[,1]  # pred is a matrix, we take the first column
  
  # Print progress
  cat("Processed predictions for model:", model_id, "\n")
}

# Convert predictions list to a data.table
predictions_dt <- as.data.table(predictions)

# Calculate metrics for each model

metrics <- lapply(names(predictions), function(model_id) {
  predicted <- predictions[[model_id]]
  actual <- actual_values[[final_models[[model_id]]$outcome]]
  
  RMSE <- sqrt(mean((predicted - actual)^2, na.rm = TRUE))
  MAE <- mae(actual, predicted)
  R2 <- cor(actual, predicted, use = "complete.obs")^2
  AUC_RROC <- auc(actual, predicted)
  
  list(RMSE = RMSE, MAE = MAE, R2 = R2, AUC_RROC = AUC_RROC)
})

# Create metrics data.table
metrics_dt <- data.table(
  model_id = names(predictions),
  outcome = sapply(names(predictions), function(id) final_models[[id]]$outcome),
  variant = sapply(names(predictions), function(id) final_models[[id]]$variant),
  RMSE = sapply(metrics, `[[`, "RMSE"),
  MAE = sapply(metrics, `[[`, "MAE"),
  R2 = sapply(metrics, `[[`, "R2"),
  AUC_RROC = sapply(metrics, `[[`, "AUC_RROC")
)

print(metrics_dt)


#save everything
write.csv(metrics_dt, file.path(output, "glmnet_INT_STRING_eigenegene_model_metrics.csv"), row.names = FALSE)
write.csv(predictions, file.path(output, "glmnet_INT_STRING_eigenegene_model_predictions.csv"), row.names = FALSE)
saveRDS(final_models, file.path(output, "all_final_glmnet_INT_STRING_eigenegene_models.rds"))





#calculate correlations and p-values for each one 

correlations <- data.frame(
  model_id = names(predictions),
  outcome = sapply(final_models, function(m) m$outcome),
  variant = sapply(final_models, function(m) m$variant)
)

correlations$Pearson_Correlation <- mapply(function(pred, out) {
  cor(pred, actual_values[[out]], use = "complete.obs")
}, predictions, correlations$outcome)

correlations$Correlation_P_Value <- mapply(function(pred, out) {
  cor.test(pred, actual_values[[out]])$p.value
}, predictions, correlations$outcome)

# Sort by correlation (descending order)
correlations <- correlations[order(-correlations$Pearson_Correlation), ]

# Save correlations as CSV
write.csv(correlations, file.path(output, "model_correlations.csv"), row.names = FALSE)
cat("Correlations saved to", file.path(output, "model_correlations.csv"), "\n")

# Print correlation results for the best model of each outcome
for (outcome in unique(correlations$outcome)) {
  best_model <- correlations[correlations$outcome == outcome, ][1, ]
  cat("\nCorrelation for best model of outcome", outcome, ":\n")
  cat("Model ID:", best_model$model_id, "\n")
  cat("Variant:", best_model$variant, "\n")
  cat("Pearson Correlation:", best_model$Pearson_Correlation, "\n")
  cat("Correlation P-Value:", best_model$Correlation_P_Value, "\n")
}



library(stats)

# Calculate Pearson correlations for each model
correlations <- data.frame(
  model_id = names(predictions),
  outcome = sapply(final_models, function(m) m$outcome),
  variant = sapply(final_models, function(m) m$variant)
)

for (i in 1:nrow(correlations)) {
  pred <- predictions[[correlations$model_id[i]]]
  actual <- actual_values[[correlations$outcome[i]]]
  
  cor_test <- cor.test(pred, actual, method = "pearson", use = "complete.obs")
  
  correlations$Pearson_Correlation[i] <- cor_test$estimate
  correlations$Correlation_P_Value[i] <- cor_test$p.value
  correlations$CI_Lower[i] <- cor_test$conf.int[1]
  correlations$CI_Upper[i] <- cor_test$conf.int[2]
  correlations$Sample_Size[i] <- length(pred)
}


# Save correlations as CSV
write.csv(correlations, file.path(output, "model_correlations.csv"), row.names = FALSE)
cat("Correlations saved to", file.path(output, "model_correlations.csv"), "\n")

# Print correlation results for the best model of each outcome
for (outcome in unique(correlations$outcome)) {
  best_model <- correlations[correlations$outcome == outcome, ][1, ]
  cat("\nCorrelation for best model of outcome", outcome, ":\n")
  cat("Model ID:", best_model$model_id, "\n")
  cat("Variant:", best_model$variant, "\n")
  cat("Pearson Correlation:", best_model$Pearson_Correlation, "\n")
  cat("95% CI: [", best_model$CI_Lower, ",", best_model$CI_Upper, "]\n")
  cat("P-Value:", best_model$Correlation_P_Value, "\n")
  cat("Sample Size:", best_model$Sample_Size, "\n")
}




























