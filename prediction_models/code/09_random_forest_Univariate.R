## =========================
## Deps
## =========================
req <- c("data.table","ranger","ggplot2","DALEX","auditor")
inst <- req[!req %in% installed.packages()[,"Package"]]
if (length(inst)) install.packages(inst, repos = "https://cloud.r-project.org")
invisible(lapply(req, library, character.only = TRUE))

## =========================
## Config
## =========================
cell_types <- if (exists("y_cols")) y_cols else c("Microglia","Astrocyte","Neuron","Oligodendrocyte")

variants <- list(
  ALL             = character(0),
  NO_PCSandAPOE   = c("APOE","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"),
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

## =========================
## Main loop over cell types
## =========================
all_metrics <- list()

for (outcome in cell_types) {
  message("== Cell type: ", outcome)
  
  needed <- unique(c("SampleID", outcome, covs, feature))
  missing <- setdiff(needed, names(a))
  if (length(missing)) stop("Missing in 'a' for outcome ", outcome, ": ", paste(missing, collapse = ", "))
  
  pheno2 <- data.frame(a[, ..needed], row.names = "SampleID", check.names = FALSE)
  pheno2 <- pheno2[complete.cases(pheno2), ]
  
  if ("Gender" %in% names(pheno2) && !is.numeric(pheno2$Gender)) {
    pheno2$Gender <- as.numeric(pheno2$Gender == "M")
  }
  if ("APOE" %in% names(pheno2) && !is.numeric(pheno2$APOE)) {
    pheno2$APOE <- ifelse(grepl("4", pheno2$APOE), 1, 0)
  }
  
  run_one <- function(tag, drops_or_keep, keep_only = FALSE) {
    df <- .build_df(pheno2, outcome, drops_or_keep, keep_only)
    
    set.seed(42)
    fit <- ranger(
      as.formula(paste(outcome, "~ .")), data = df,
      num.trees = hp$num.trees, mtry = hp$mtry, min.node.size = hp$min.node.size,
      replace = hp$replace, sample.fraction = hp$sample.fraction,
      importance = "permutation", scale.permutation.importance = TRUE,
      respect.unordered.factors = "partition", seed = 42
    )
    
    y_obs    <- df[[outcome]]
    yhat_oob <- as.numeric(fit$predictions)  # OOB preds
    rmse_oob <- sqrt(mean((y_obs - yhat_oob)^2))
    rse_oob  <- rse_safe(y_obs, yhat_oob)
    
    dr      <- rroc_curve(y_obs, yhat_oob)
    aoc_val <- rroc_aoc(dr)
    
    # Auditor-styled RROC from OOB preds
    expl <- DALEX::explain(
      model = list(),  # dummy; predict returns OOB preds
      data  = data.frame(.oof_id = seq_along(yhat_oob)),
      y     = y_obs,
      predict_function = function(m, newdata) yhat_oob[newdata$.oof_id],
      label = sprintf("%s | %s (ranger OOB)", outcome, tag),
      verbose = FALSE
    )
    resid_obj <- auditor::model_residual(expl)
    p_rroc    <- auditor::plot_rroc(resid_obj)
    
    # Save PNG per (cell type × variant)
    ggsave(file.path(output, "rroc", sprintf("rroc_%s_%s.png", outcome, tag)),
           p_rroc, width = 7, height = 6, dpi = 300)
    
    list(tag = tag, rmse_oob = rmse_oob, rse_oob = rse_oob, aoc_rroc = aoc_val, rroc_plot = p_rroc)
  }
  
  results <- list(
    run_one("ALL",             variants$ALL,             keep_only = FALSE),
    run_one("NO_PCSandAPOE",   variants$NO_PCSandAPOE,   keep_only = FALSE),
    run_one("NO_PCS", variants$NO_PCS, keep_only = FALSE)
  )
  
  # Per-outcome metrics CSV
  tbl <- rbindlist(lapply(results, function(x)
    data.table(
      variant  = x$tag,
      rmse_oob = x$rmse_oob,
      rse_oob  = x$rse_oob,
      aoc_rroc = x$aoc_rroc   
    )
  ))
  data.table::setorder(tbl, rmse_oob)
  data.table::fwrite(tbl, file.path(output, sprintf("rf_oob_metrics_%s.csv", outcome)))
  
  # One PDF per outcome (one page per variant, auditor look)
  pdf_path <- file.path(output, sprintf("rroc_%s.pdf", outcome))
  pdf(pdf_path, width = 7, height = 6)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)
  for (res in results) print(res$rroc_plot)
  dev.off()
  
  tbl[, outcome := outcome]
  all_metrics[[outcome]] <- tbl
}

# Macro summary (average across variants per outcome)
macro <- rbindlist(all_metrics)[, .(
  macro_rmse = mean(rmse_oob, na.rm = TRUE),
  macro_rse  = mean(rse_oob,  na.rm = TRUE), 
  macro_aoc  = mean(aoc_rroc,  na.rm = TRUE)
), by = .(outcome)]
fwrite(macro, file.path(output, "rf_macro_metrics_by_outcome.csv"))


fwrite(macro_by_variant, file.path(output, "rf_macro_metrics_by_variant.csv"))
print(macro_by_variant)

cat("Done. Outputs in ", normalizePath(output),
    " (RROC PNGs in rroc/, RROC PDF per outcome, and CSVs).", sep = "")



#Run actual model

pheno2_filtered <- pheno2[, !(names(pheno2) %in% c("Microglia", "Astrocyte", "Oligodendrocyte"))]

rf_neuron_INT <- ranger(
  Neuron ~ ., data = pheno2_norm,
  num.trees = hp$num.trees, mtry = hp$mtry, min.node.size = hp$min.node.size,
  replace = hp$replace, sample.fraction = hp$sample.fraction,
  importance = "permutation", scale.permutation.importance = TRUE,
  respect.unordered.factors = "partition", seed = 42
)

neuron_predict <- predict(rf_neuron_INT, data = test_data)

#Organize predictions so everything is aligned 
actual_values <- as.data.frame(pheno2["Neuron"])
pred_neuron <- as.data.frame(neuron_predict$predictions)
rownames(pred_neuron) <- if (!is.null(rownames(test_data))) rownames(test_data) else test_data[, "SampleID"]
colnames(pred_neuron) <- "Neuron"
#sanity check 
identical(colnames(pred_neuron), colnames(actual_values))     
identical(rownames(pred_neuron), rownames(test_data)) 

metrics <- data.frame(
  Outcome = colnames(actual_values),
  RMSE = NA,
  MAE  = NA,
  MAPE = NA,
  R2   = NA
)


for (i in seq_along(metrics$Outcome)) {
  y_true <- actual_values[, i]
  y_pred <- pred_neuron[, i]
  
  metrics$RMSE[i] <- sqrt(mean((y_true - y_pred)^2))
  metrics$MAE[i]  <- mean(abs(y_true - y_pred))
  metrics$MAPE[i] <- mean(abs((y_true - y_pred) / y_true)) * 100  # percentage
  metrics$R2[i]   <- cor(y_true, y_pred)^2
}

metrics

#merge actual and predicted values

common_rows <- intersect(rownames(actual_values), rownames(pred_neuron))
common_cols <- intersect(colnames(actual_values), colnames(pred_neuron))



A <- actual_values[common_rows, common_cols, drop = FALSE]
P <- pred_neuron[common_rows,  common_cols, drop = FALSE]

# add predictions next to actuals, with “…_pred” names
actual_plus_pred <- cbind(
  A,
  setNames(P, paste0(common_cols, "_pred"))
)

merged_data <- merge(pred_neuron, actual_values, by = "row.names", all.x = TRUE)


cells <- c("Neuron")  # edit to your set

for (cell in cells) {
  xcol <- paste0(cell, ".x")
  ycol <- paste0(cell, ".y")
  
  p <- ggplot(merged_data, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_point(alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    labs(title = paste("Predicted vs Observed —", cell),
         x = "Observed", y = "Predicted") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Save each plot with explicit dimensions
  ggsave(file.path(output, paste0("rft__INT_predicted_vs_observed_", cell, ".pdf")), 
         plot = p, width = 7, height = 7, dpi = 300)
}







