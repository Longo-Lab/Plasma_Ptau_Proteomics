library(dplyr)
col_names <- x[]

target_cols <- colnames(x)[c(2:28)]


outcome <- paste0(outcome, "_q")

a <- a %>%
  mutate(across(all_of(outcome),
                ~ ntile(., 4),
                .names = "{.col}_q"))


pheno2[[outcome]] <- factor(pheno2[[outcome]])

pheno2[outcome] <- lapply(pheno2[outcome], factor)

pheno2 <- pheno

variants <- list(
  ALL = character(0),
  NO_hybridNormControlScale = "HybControlNormScale",
  NO_PCSandAPOE = c(
    "APOE.e34", "APOE.e44", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"
  ),
  NO_PCS = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8")
)

# Join all predictors from pheno2 and get rid of cell types 
predictor_pool <- setdiff(names(pheno2), outcome)

# Tuning grid that tests different alphas and lambdas
tuning_grid <- expand.grid(
  alpha  = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),  
  lambda = seq(0.001, 0.1, length = 20)
)


## ---- CV folds that will normalize the different runs ----
set.seed(42)
K <- 50
folds_train <- createFolds(seq_len(nrow(pheno2)), k = K, returnTrain = TRUE)
index    <- folds_train
indexOut <- lapply(folds_train, function(tr) setdiff(seq_len(nrow(pheno2)), tr))


## ---- loop: outcome × variant with RMSE, RSE, RROC/AOC ------------------------------
rows <- list()
coef_txt <- file.path(output, "coefficients_printouts.txt")
if (file.exists(coef_txt)) file.remove(coef_txt)

# PDF device that will collect all auditor RROCs (one page per (outcome×variant))
# rroc_pdf <- file.path(output, "rroc_curves.pdf")
# pdf(rroc_pdf, width = 7, height = 6)
# on.exit(try(dev.off(), silent = TRUE), add = TRUE)

for (outcome in outcome) {
  for (tag in names(variants)) {
    # Build predictor set for this variant (drop-list only)
    drops <- variants[[tag]]
    missing_drops <- setdiff(drops, names(pheno2))
    if (length(missing_drops)) {
      warning(sprintf("[%s|%s] Drop list contains missing columns: %s",
                      outcome, tag, paste(missing_drops, collapse = ", ")))
    }
    x_cols <- setdiff(predictor_pool, intersect(drops, names(pheno2)))
    
    if (length(x_cols) == 0) {
      warning(sprintf("[%s|%s] No predictors after dropping; skipping.", outcome, tag))
      next
    }
    
    # Restrict data to outcome + selected predictors
    df <- pheno2[, c(outcome, x_cols), drop = FALSE]
    
    
    # Train-----------------------------------------------------------------------------------------------------
    set.seed(42)
    
    cv_model <- train(
      as.formula(paste(outcome, "~ .")),
      data = df,
      method = "glmnet",
      trControl = trainControl(method = "cv", number = K,
                               index = index, indexOut = indexOut,
                               savePredictions = "final", verboseIter = FALSE),
      tuneGrid = tuning_grid,
      metric = "Accuracy"
    )
    
    # Best tune + OOF predictions
    bt  <- cv_model$bestTune
    bundle <- list(
      outcome = outcome,
      variant = tag,
      model   = cv_model$finalModel,   # glmnet object (path of lambdas)
      lambda  = bt$lambda,             # chosen lambda
      alpha   = bt$alpha,
      x_cols  = x_cols                 # exact predictor set/order used
    ) 
    
    #saveRDS(bundle, file = file.path(output, sprintf("glmnet_%s_%s.rds", outcome, tag)))
    
    oof <- subset(cv_model$pred, alpha == bt$alpha & lambda == bt$lambda)
    
    # Guard against empty OOF (can happen if factor levels / NA mismatch)
    # if (nrow(oof) == 0) {
    #     warning(sprintf("[%s|%s] No OOF rows after filtering bestTune; skipping.", outcome, tag))
    #     next
    # }
    #     
    # # Metrics: RMSE, RSE (safe), R2, MAE, MAPE
    acc_val <- mean(oof$pred == oof$obs, na.rm = TRUE)
    # rmse_val <- RMSE(oof$pred, oof$obs)
    # rse_val  <- rse_safe(oof$obs, oof$pred)
    # r2_val   <- 1 - rse_val
    # mae_val <- mean(abs(oof$pred - oof$obs))
    # mape_val <- mean(abs((oof$pred - oof$obs) / oof$obs)) * 100
    
    # ---- auditor RROC from OOF predictions (saved PNG + printed into the PDF) ----
    # expl <- DALEX::explain(
    #   model = list(),  # dummy
    #   data  = data.frame(.oof_id = seq_len(nrow(oof))),
    #   y     = oof$obs,
    #   predict_function = function(m, newdata) oof$pred[newdata$.oof_id],
    #   label = sprintf("%s | %s (glmnet OOF)", outcome, tag),
    #   verbose = FALSE
    # )
    # 
    # resid_obj <- auditor::model_residual(expl)
    # p_rroc <- auditor::plot_rroc(resid_obj) +
    #   ggplot2::ggtitle(sprintf("RROC — %s | %s\nRMSE=%.3f  RSE=%.3f",
    #                            outcome, tag, rmse_val, rse_val))
    
    # Save a PNG per (outcome × variant)
    
    # png_path <- file.path(output, "rroc", sprintf("rroc_%s_%s.png", outcome, tag))
    # 
    # ggsave(png_path, p_rroc, width = 7, height = 6, dpi = 300)
    # 
    # # Also print into the open PDF device
    # print(p_rroc)
    # 
    # #Scalar AOC from our simple RROC to keep your table column
    # d_rroc  <- rroc_curve(oof$obs, oof$pred)
    # aoc_val <- rroc_aoc(d_rroc)
    
    
    # Collect row
    rows[[length(rows) + 1]] <- data.table(
      outcome     = outcome,
      variant     = tag,
      best_alpha  = bt$alpha,
      best_lambda = bt$lambda,
      accuracy    = acc_val,
      # r2_oof      = r2_val,
      # rse_oof     = rse_val,
      # aoc_rroc    = aoc_val,
      # mae_val     = mae_val,
      # mape_val    = mape_val,
      n_features  = length(x_cols)
    )
    
    # ---- Coefficients (raw printout + tidy CSV of non-zeros) ----
    # co <- coef(cv_model$finalModel, s = bt$lambda)
    # 
    # cat(sprintf("\n===== %s | %s =====\n", outcome, tag))
    # print(co)
    # 
    # cat(sprintf("\n===== %s | %s =====\n", outcome, tag), file = coef_txt, append = TRUE)
    # capture.output(print(co), file = coef_txt, append = TRUE)
    # 
    # cz <- as.matrix(co)
    # coef_df <- data.table(feature = rownames(cz), coef = as.numeric(cz))
    # coef_df_nz <- coef_df[coef != 0]
    # if (nrow(coef_df_nz)) {
    #   coef_df_nz[, `:=`(outcome = outcome, variant = tag)]
    #   fwrite(
    #     coef_df_nz,
    #     file = file.path(output, "coefficients_nonzero.csv"),
    #     append = file.exists(file.path(output, "coefficients_nonzero.csv"))
    #   )
    }
  }
}

## ---- save summaries ----
results <- rbindlist(rows)
fwrite(results, file.path(output, "glmnet_outcomes_quartiles.csv"))

# Macro CSV now includes mean RSE (macro_rse) and mean AOC (macro_aoc)
# macro <- results[, .(
#   macro_rmse = mean(rmse_oof, na.rm = TRUE),
#   macro_r2   = mean(r2_oof,   na.rm = TRUE),
#   macro_rse  = mean(rse_oof,  na.rm = TRUE),
#   macro_aoc  = mean(aoc_rroc, na.rm = TRUE)
# ), by = variant][order(macro_rmse)]
# fwrite(macro, file.path(output, "macro_by_variant.csv"))

# cat("\nSaved:\n",
#     "- metrics: ", file.path(output, "glmnet_outcomes_by_variant.csv"), "\n",
#     # "- macro:   ", file.path(output, "macro_by_variant.csv"), "\n",
#     "- coefs:   ", file.path(output, "coefficients_printouts.txt"),
#     " (raw coef() printouts)\n",
#     "- coefs:   ", file.path(output, "coefficients_nonzero.csv"),
#     " (tidy non-zero)\n",
#     # "- RROC PNGs in: ", file.path(output, "rroc"), "\n",
#     # "- RROC PDF: ", rroc_pdf, "\n", sep = "")
# 
# # make sure no open devices remain
# dev.list()
# while (!is.null(dev.list())) dev.off()

#Test 


best <- results[order(-accuracy)][, .SD[1], by = outcome]

final_models <- vector("list", nrow(best))
names(final_models) <- best$outcome

for (i in seq_len(nrow(best))) {
  
  out_q <- best$outcome[i]        # should be the *_q column name
  tag   <- best$variant[i]        # best variant for that outcome
  a     <- best$best_alpha[i]
  lam   <- best$best_lambda[i]
  
  # predictors for THIS outcome's best variant
  drops  <- variants[[tag]]
  x_cols <- setdiff(predictor_pool, intersect(drops, names(pheno2)))
  
  # training set for this outcome
  df <- pheno2[, c(out_q, x_cols), drop = FALSE]
  df <- df[!is.na(df[[out_q]]), , drop = FALSE]
  
  y <- factor(df[[out_q]])
  if (nlevels(y) < 2) {
    warning("Skipping ", out_q, " (needs >=2 classes).")
    next
  }
  
  X <- as.matrix(df[, x_cols, drop = FALSE])
  
  final_models[[out_q]] <- list(
    model   = glmnet(X, y, family = "multinomial", alpha = a, lambda = c(lam)),
    outcome = out_q,
    variant = tag,
    alpha   = a,
    lambda  = lam,
    x_cols  = x_cols,
    levels  = levels(y)
  )
}

## Predict


a_test <- a


names(a) <- names(name_map <- setNames(names(a), make.names(names(a), unique = TRUE)))

# All variables need to be numeric & dummy coded
a[, Gender := as.numeric(Gender == "M")]
recode_dummies(a, "APOE", keep = c("e34", "e44"))
covs <- unique(c(covs[covs != "APOE"], "APOE.e34", "APOE.e44"))

ids <- a[order(DS1V)][, head(.SD, 1), by = "PIDN"][, SampleID]
f_ups <- a[!SampleID %in% ids]
columns_to_check <- c( "PlasmaPTau181","PlasmaAB142P","PlasmaAB140P","PlasmaABRatio","PlasmapTau217","PlasmapTau217_AB42Ratio",
                       "PlasmaGFAP","PlasmaNfL")
drops <- f_ups[!complete.cases(f_ups[, ..columns_to_check ]), SampleID]
a_test <- f_ups[!SampleID %in% drops]

outcome <- make.names(setdiff(names(fread(prediction_scores, nrows = 0)), c("SampleID", "Neddylation")), unique = TRUE)


a_test <- a_test %>%
  mutate(across(all_of(outcome),
                ~ ntile(., 4),
                .names = "{.col}_q"))

test_data <- data.frame(a_test[, ..cols], row.names = "SampleID")
test_data <- test_data[complete.cases(test_data),] 

outcome <- paste0(outcome, "_q")

test_data[outcome] <- lapply(test_data[outcome], factor)


outs <- names(final_models)
if (is.null(outs) || any(outs == "")) {
  # fallback: pull outcome names stored inside each model object (if you saved them)
  outs <- sapply(final_models, function(obj) obj$outcome %||% NA_character_)
}

pred_list <- vector("list", length(outs))
names(pred_list) <- outs

for (out_q in outs) {
  obj <- final_models[[out_q]]
  if (is.null(obj)) next
  
  X_test <- as.matrix(test_data[, obj$x_cols, drop = FALSE])
  pred_list[[out_q]] <- as.character(
    predict(obj$model, newx = X_test, s = obj$lambda, type = "class")[, 1]
  )
}

# drop NULL outcomes
pred_list <- pred_list[!sapply(pred_list, is.null)]

pred_class <- as.data.frame(do.call(cbind, pred_list))
colnames(pred_class) <- names(pred_list)   # force correct names
rownames(pred_class) <- rownames(test_data)


accuracy_per_outcome <- sapply(colnames(pred_class), function(out_q) {
  mean(pred_class[[out_q]] == test_data[[out_q]], na.rm = TRUE)
})

accuracy_per_outcome


accuracy_df <- data.frame(
  outcome = colnames(pred_class),
  accuracy = sapply(colnames(pred_class), function(out_q) {
    mean(pred_class[[out_q]] == test_data[[out_q]], na.rm = TRUE)
  }),
  row.names = NULL
)

accuracy_df

write.csv(
  accuracy_df,
  file = file.path(output, "test_accuracy_by_outcome.csv"),
  row.names = FALSE
)

md_lines <- c(
  "| Outcome | Accuracy |",
  "|---------|----------|",
  paste0(
    "| ",
    accuracy_df$outcome,
    " | ",
    sprintf("%.3f", accuracy_df$accuracy),
    " |"
  )
)

md_text <- paste(md_lines, collapse = "\n")
cat(md_text)



