##
inpath <- "C:/Users/audweber/Box/plasma_ptau_proteomics" # don't repeat things you don't have to
input <- file.path(inpath, "dataProt_SS-205063.merged.QC.csv")
prots <- file.path(inpath, "prot_names.txt") 
pcs <- file.path(inpath, "genomic/mypcair.PCs.csv")
prediction_scores <- file.path(inpath, # wrap at 80 characters
                               "celltype_scores/dataProt_SS-205063.SEA-AD.Neighborhood.1kDS.lvl2.normalized_scores.csv"
)
output <- 'prediction_models/cluster/rf'
covs <- 'Age,Gender,APOE,HybControlNormScale'   
outcome <- c("APP.Metabolism", "Apoptosis",	"Autophagy", "Cell.Cycle",	"DNA.Repair",	"Endolysosome",	"Epigenetic",	
             "Immune.Response",	"Lipid.Metabolism",	"Metal.Binding.and.Homeostasis",	"Mitochondrial.Metabolism",	"Myelination",	
             "Oxidative.Stress",	"Proteostasis",	"RNA.Spliceosome",	"Structural.Stabilization",	"Synapse",	"Tau.Homeostasis",	"Vasculature")
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

cell_types <- outcome


##=========================
## Deps
## =========================
req <- c("data.table","ranger","ggplot2","DALEX","auditor")
inst <- req[!req %in% installed.packages()[,"Package"]]
if (length(inst)) install.packages(inst, repos = "https://cloud.r-project.org")
invisible(lapply(req, library, character.only = TRUE))

## =========================
## Config
## =========================
cell_types <- outcome

variants <- list(
  ALL = character(0),
  NO_hybridNormControlScale = "HybControlNormScale",
  NO_PCSandAPOE = c(
    "APOE.e34", "APOE.e44", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"
  ),
  NO_PCS = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8")
)

hp <- list(num.trees = 1500, mtry = 2, min.node.size = 1, sample.fraction = 0.50, replace = TRUE)

# Use your existing 'output' folder; create needed subfolders
if (!exists("output")) stop("Please define `output` as a folder path (e.g., output <- 'prediction_model').")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output, "rroc"), recursive = TRUE, showWarnings = FALSE)

## =========================
## Helpers
## =========================
if (is.character(covs) && length(covs) == 1 && grepl(",", covs)) {
  covs <- trimws(strsplit(covs, ",")[[1]])
}

rroc_curve <- function(y, yhat, n_points = 400) {
  ok <- is.finite(y) & is.finite(yhat)
  y <- y[ok]; yhat <- yhat[ok]
  e <- yhat - y
  s_min <- -max(e); s_max <- -min(e)
  shifts <- seq(s_min - 1e-8, s_max + 1e-8, length.out = max(20, n_points))
  UNDER <- OVER <- numeric(length(shifts))
  for (i in seq_along(shifts)) {
    ee <- e + shifts[i]
    OVER[i]  <- sum(ee[ee > 0])      # x
    UNDER[i] <- sum(-ee[ee < 0])     # y
  }
  data.frame(OVER = OVER, UNDER = UNDER)
}

rroc_aoc <- function(rroc_df) {
  d <- rroc_df[order(rroc_df$OVER), ]   # <- fix here
  if (nrow(d) < 2) return(NA_real_)
  sum((d$UNDER[-1] + d$UNDER[-nrow(d)]) * diff(d$OVER) / 2)
}

rse_safe <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]
  sse <- sum((obs - pred)^2)
  tss <- sum((obs - mean(obs))^2)
  if (tss == 0) return(NA_real_)
  sse / tss
}

.build_df <- function(ph, outcome, drops_or_keep, keep_only = FALSE) {
  pool <- setdiff(names(ph), outcome)
  if (keep_only) {
    keep <- intersect(drops_or_keep, pool)
    if (!length(keep)) stop("Keep-only set left 0 predictors; check names.")
    ph[, c(outcome, keep), drop = FALSE]
  } else {
    drops <- intersect(drops_or_keep, names(ph))
    xcols <- setdiff(pool, drops)
    if (!length(xcols)) stop("Drop set removed all predictors; check names.")
    ph[, c(outcome, xcols), drop = FALSE]
  }
}
cell_types = outcome
## =========================
## Main loop over cell types
## =========================
all_metrics <- list()
# 
# for (outcome in cell_types) {
#   message("== Cell type: ", outcome)
#   
#   # needed <- unique(c("SampleID", outcome, covs, feature))
  # missing <- setdiff(needed, names(pheno2))
  # if (length(missing)) stop("Missing in 'a' for outcome ", outcome, ": ", paste(missing, collapse = ", "))
  # 
  # if ("Gender" %in% names(pheno2) && !is.numeric(pheno2$Gender)) {
  #   pheno2$Gender <- as.numeric(pheno2$Gender == "M")
  # }

  
run_one <- function(tag, drops_or_keep, keep_only = FALSE) {
    df <- .build_df(cluster2, outcome, drops_or_keep, keep_only)
    set.seed(42)
    
    fit <- ranger(
      as.formula(paste(outcome, "~ .")), data = df,
      
      # Core parameters
      num.trees = hp$num.trees,
      mtry = hp$mtry,
      min.node.size = max(1, hp$min.node.size),   # Ensure >= 1 for classification
      replace = hp$replace,
      sample.fraction = hp$sample.fraction,
      splitrule = "gini",                         # Good default for classification
      respect.unordered.factors = "partition",
      importance = "permutation",
      scale.permutation.importance = TRUE
    )
    
    result_row <- process_ranger_results(fit, outcome, tag, x_cols, output)
    
    return(result_row)
  }
  
  # Modified function for direct ranger() objects
  process_ranger_results <- function(fit, outcome, tag, x_cols, output) {
    
    # ---- Save model components (no bestTune for direct ranger) ----
    bundle <- list(
      outcome = outcome,
      variant = tag,
      model   = fit,                    # The ranger object itself
      mtry    = hp$mtry,               # From hyperparameters
      splitrule = "gini",              # From your settings
      min.node.size = hp$min.node.size,
      num.trees = fit$num.trees,
      x_cols  = x_cols
    )
    
    saveRDS(bundle, file = file.path(output, sprintf("rf_%s_%s.rds", outcome, tag)))
    
    # ---- No OOF predictions available with direct ranger ----
    # We'll use OOB (Out-of-Bag) predictions instead
    oob_pred <- fit$predictions
    
    if (is.null(oob_pred)) {
      warning(sprintf("[%s|%s] No OOB predictions available.", outcome, tag))
      return(NULL)
    }
    
    # For ranger, we don't have cross-validation accuracy, so we'll use OOB error
    oob_accuracy <- 1 - fit$prediction.error
    
    # ---- Variable Importance ----
    var_imp <- importance(fit)
    var_imp_df <- data.frame(
      variable = names(var_imp),
      importance = as.numeric(var_imp)
    ) %>%
      arrange(desc(importance))
    
    # Save variable importance
    fwrite(var_imp_df, file.path(output, sprintf("rf_importance_%s_%s.csv", outcome, tag)))
    
    # ---- Print Variable Importance ----
    cat("\n=== Variable Importance ===\n", file = file.path(output, "rf_importance_printouts.txt"), append = TRUE)
    cat(sprintf("Outcome: %s | Variant: %s\n", outcome, tag), 
        file = file.path(output, "rf_importance_printouts.txt"), append = TRUE)
    
    top_vars <- head(var_imp_df, 20)
    for(i in 1:nrow(top_vars)) {
      cat(sprintf("%2d. %-30s %8.4f\n", i, top_vars$variable[i], top_vars$importance[i]), 
          file = file.path(output, "rf_importance_printouts.txt"), append = TRUE)
    }
    cat("\n", file = file.path(output, "rf_importance_printouts.txt"), append = TRUE)
    
    # ---- Collect results row ----
    result_row <- data.table(
      outcome        = outcome,
      variant        = tag,
      mtry           = hp$mtry,
      splitrule      = "gini",
      min_node_size  = hp$min.node.size,
      num_trees      = fit$num.trees,
      oob_accuracy   = oob_accuracy,        # OOB instead of CV accuracy
      oob_error      = fit$prediction.error,
      n_features     = length(x_cols),
      top_var_importance = var_imp_df$importance[1],
      mean_var_importance = mean(var_imp_df$importance)
    )
    
    return(result_row)
  }
  
  # Run for different variants (you need to define your actual variants)
  variants <- c("variant1", "variant2")  # Replace with your actual variant names
  
  for (tag in variants) {
    result <- run_one(tag, drops_or_keep, keep_only = FALSE)  # Define these parameters
    if (!is.null(result)) {
      all_metrics[[length(all_metrics) + 1]] <- result
    }
  }
}

# ---- Save final results ----
if (length(all_metrics) > 0) {
  results <- rbindlist(all_metrics)
  fwrite(results, file.path(output, "rf_outcomes_by_variant.csv"))
  
  cat("\nSaved:\n",
      "- metrics: ", file.path(output, "rf_outcomes_by_variant.csv"), "\n",
      "- models:  ", file.path(output, "rf_*.rds"), "\n", 
      "- importance: ", file.path(output, "rf_importance_printouts.txt"), "\n",
      sep = "")
}

    
      
    # set.seed(42)
    # fit <- ranger(
    #   as.formula(paste(outcome, "~ .")), data = df,
    #   num.trees = hp$num.trees, mtry = hp$mtry, min.node.size = hp$min.node.size,
    #   replace = hp$replace, sample.fraction = hp$sample.fraction,
    #   importance = "permutation", scale.permutation.importance = TRUE,
    #   respect.unordered.factors = "partition", seed = 42
    # )
    
    # y_obs    <- df[[outcome]]
    # yhat_oob <- as.numeric(fit$predictions)  # OOB preds
    # rmse_oob <- sqrt(mean((y_obs - yhat_oob)^2))
    # rse_oob  <- rse_safe(y_obs, yhat_oob)
    # dr      <- rroc_curve(y_obs, yhat_oob)
    # aoc_val <- rroc_aoc(dr)
    # mae_val <- mean(abs(yhat_oob - y_obs))
    # mape_val <- mean(abs((yhat_oob - y_obs) / y_obs)) * 100
    # mape_new <- calculate_mape(yhat_oob,y_obs)
    
    
    # Auditor-styled RROC from OOB preds
#     expl <- DALEX::explain(
#       model = list(),  # dummy; predict returns OOB preds
#       data  = data.frame(.oof_id = seq_along(yhat_oob)),
#       y     = y_obs,
#       predict_function = function(m, newdata) yhat_oob[newdata$.oof_id],
#       label = sprintf("%s | %s (ranger OOB)", outcome, tag),
#       verbose = FALSE
#     )
#     resid_obj <- auditor::model_residual(expl)
#     p_rroc    <- auditor::plot_rroc(resid_obj)
#     
#     # Save PNG per (cell type × variant)
#     ggsave(file.path(output, "rroc", sprintf("rroc_%s_%s.png", outcome, tag)),
#            p_rroc, width = 7, height = 6, dpi = 300)
#     
#     list(tag = tag, rmse_oob = rmse_oob, rse_oob = rse_oob, aoc_rroc = aoc_val, rroc_plot = p_rroc, mae = mae_val, mape = mape_val)
#   }
#   
#   results <- list(
#     run_one("ALL",             variants$ALL,             keep_only = FALSE),
#     run_one("NO_PCSandAPOE",   variants$NO_PCSandAPOE,   keep_only = FALSE),
#     run_one("NO_PCS", variants$NO_PCS, keep_only = FALSE),
#     run_one("NO_hybridNormControlScale", variants$NO_hybridNormControlScale, keep_only = FALSE),
#     run_one("Biomarkers", variants$Biomarkers, keep_only = FALSE)
#   )
#   
#   # Per-outcome metrics CSV
#   tbl <- rbindlist(lapply(results, function(x)
#     data.table(
#       variant  = x$tag,
#       # rmse_oob = x$rmse_oob,
#       # rse_oob  = x$rse_oob,
#       # aoc_rroc = x$aoc_rroc, 
#       # mae      = x$mae,
#       # mape     = x$mape,
#       # mape_new     = x$mape_new
#     )
#   ))
#   data.table::setorder(tbl, rmse_oob)
#   data.table::fwrite(tbl, file.path(output, sprintf("rf_oob_metrics_%s.csv", outcome)))
#   
#   # One PDF per outcome (one page per variant, auditor look)
#   pdf_path <- file.path(output, sprintf("rroc_%s.pdf", outcome))
#   pdf(pdf_path, width = 7, height = 6)
#   on.exit(try(dev.off(), silent = TRUE), add = TRUE)
#   for (res in results) print(res$rroc_plot)
#   dev.off()
#   
#   tbl[, outcome := outcome]
#   all_metrics[[outcome]] <- tbl
# 
# }
# Macro summary (average across variants per outcome)
# macro <- rbindlist(all_metrics)[, .(
#   macro_rmse = mean(rmse_oob, na.rm = TRUE),
#   macro_rse  = mean(rse_oob,  na.rm = TRUE), 
#   macro_aoc  = mean(aoc_rroc,  na.rm = TRUE), 
#   macro_mae  = mean(mae,  na.rm = TRUE),
#   macro_mape  = mean(mape,  na.rm = TRUE)
# ), by = .(outcome)]
# fwrite(macro, file.path(output, "rf_macro_metrics_by_outcome.csv"))
# 
# macro_by_variant <- rbindlist(all_metrics)[, .(
#   mean_rmse = mean(rmse_oob, na.rm = TRUE),
#   mean_rse  = mean(rse_oob,  na.rm = TRUE),
#   mean_aoc  = mean(aoc_rroc, na.rm = TRUE),
#   mean_mae  = mean(mae,      na.rm = TRUE),
#   mean_mape = mean(mape,     na.rm = TRUE)
# ), by = .(variant)]
# fwrite(macro_by_variant, file.path(output, "rf_macro_metrics_by_variant.csv"))
# print(macro_by_variant)
# 
# cat("Done. Outputs in ", normalizePath(output),
#     " (RROC PNGs in rroc/, RROC PDF per outcome, and CSVs).", sep = "")
# 
# }
##Run on ideal model: 
#NOPCorAPOE + INT for Vascular 

#pick one cell type
others <- make.names(setdiff(names(fread(prediction_scores, nrows = 0)), c( "MEMerged.M2")), unique = TRUE)
others <- setdiff(others, "MEMerged.M2")
pc_apoe <- c(
  "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8",
  "APOE","APOE.e34","APOE.e44"
)

other_cells <- unique(c(others, pc_apoe))
pheno3 <- pheno2[, setdiff(names(pheno2), other_cells), drop = FALSE]
pheno3 <- cluster2[, c(3, 15:22)]

#Best model
rf_cluster <- ranger(
  clust ~ ., data = pheno3,
    num.trees = hp$num.trees, mtry = hp$mtry, min.node.size = hp$min.node.size,
    replace = hp$replace, sample.fraction = hp$sample.fraction,
    importance = "permutation", scale.permutation.importance = TRUE,
    respect.unordered.factors = "partition", seed = 42
)

saveRDS(rf_cluster, file.path(output, "final_rf_cluster.rds"))


#Set up test data 
a <- fread(input)
a_test <- a

ids <- a_test[order(DS1V)][, head(.SD, 1), by = "PIDN"][, SampleID]
f_ups <- a_test[!SampleID %in% ids]
columns_to_check <- c(  "Age","Gender","HybControlNormScale", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","APOE.e34","APOE.e44",
                        "PlasmaPTau181","PlasmaAB142P","PlasmaAB140P","PlasmaABRatio","PlasmapTau217","PlasmapTau217_AB42Ratio",
                        "PlasmaGFAP","PlasmaNfL")
drops <- f_ups[!complete.cases(f_ups[, ..columns_to_check ]), SampleID]
a_test <- f_ups[!SampleID %in% drops]

test_data <- data.frame(a_test[, ..cols], row.names = "SampleID")
# test_data2 <- pheno2_norm


actual_values <- as.data.frame(f_ups2[, c("clust")])
test_data2 <- test_data[, c(80:83, 93:100)]
# test_data2 <- test_data[, setdiff(names(test_data), other_cells), drop = FALSE]

#Time to predict! 

rt_cluster_pred <- predict(rf_cluster, data = new_data)


y_pred <- as.numeric(rt_cluster_pred$predictions)

# Grab actuals 
y_true <- as.numeric(f_ups2$clust)  

# make sure lengths match
stopifnot(length(y_true) == length(y_pred))

#Metrics
rmse <- sqrt(mean((y_true - y_pred)^2, na.rm = TRUE))
mae  <- mean(abs(y_true - y_pred), na.rm = TRUE)
den  <- pmax(abs(y_true), 1e-8)                        # safe MAPE denom
mape <- mean(abs(y_true - y_pred) / den, na.rm = TRUE) * 100
sse  <- sum((y_true - y_pred)^2, na.rm = TRUE)
tss  <- sum((y_true - mean(y_true, na.rm = TRUE))^2, na.rm = TRUE)
r2   <- if (tss == 0) NA_real_ else 1 - sse / tss
rse <- sse /tss
MAPE_new <- calculate_mape(y_true, y_pred)

metrics <- data.frame(RMSE = rmse, MAE = mae, MAPE_before = mape, R2 = r2, RSE = rse, MAPE_new= MAPE_new)

write.csv(metrics,
          file.path(output, "metrics_rf_cluster.csv"),
          row.names = FALSE)


print(metrics)

#save predictions 

pred_df <- data.frame(
  Predicted = rt_cluster_pred
)

merged_test <- cbind(test_data2, pred_df)
idx <- match(rownames(rt_MEMerged.M2_pred), rownames(test_data2))
rf_merged_data <- cbind(
  y_pred,
  test_data2[, drop = FALSE]   # will insert NAs where no match
)

merged_test <- cbind(
  test_data2,
  true_MEMerged.M2 = actual_values,
  Predicted_MEMerged.M2 = rt_MEMerged.M2_pred
)

write.csv(
  merged_test,
  file.path(output, "predictions_MEMerged.M2_rf.csv"),
  row.names = FALSE
)

# Plot Predicted vs Observed
# df_plot <- data.frame(Observed = y_true, Predicted = y_pred)
# p <- ggplot(df_plot, aes(x = Observed, y = Predicted)) +
#   geom_point(alpha = 0.7) +
#   geom_abline(intercept = 0, slope = 1, linetype = 2) +
#   coord_equal() +
#   labs(title = "Predicted vs Observed — RHO.GTPase.Neurotrophin.signaling",
#        x = "Observed", y = "Predicted") +
#   theme_minimal()
# 
# print(p)
# 
# ggsave(file.path(output, paste0("predicted_vs_observed_RHO.GTPase.Neurotrophin.signaling", ".pdf")), 
#        plot = p, width = 7, height = 7, dpi = 300)




