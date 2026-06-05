y_cols <- outcome

variants <- list(
  ALL = character(0),
  NO_PCSandAPOE = c("Age","Gender","PlasmaPTau181_norm","PlasmaAB142P_norm","PlasmaAB140P_norm","PlasmaABRatio_norm",
                    "PlasmapTau217_norm","PlasmapTau217_AB42Ratio_norm","PlasmaGFAP_norm","PlasmaNfL_norm","HybControlNormScale"),
  NO_PCS = c("Age","Gender","APOE.e34","APOE.e44","PlasmaPTau181_norm","PlasmaAB142P_norm","PlasmaAB140P_norm","PlasmaABRatio_norm",
             "PlasmapTau217_norm","PlasmapTau217_AB42Ratio_norm","PlasmaGFAP_norm","PlasmaNfL_norm","HybControlNormScale"),
  Impt_biomarkers = c("PlasmaPTau181_norm","PlasmaABRatio_norm","PlasmaGFAP_norm","PlasmaNfL_norm","Age","Gender","APOE.e34","APOE.e44",
                      "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","HybControlNormScale")
)

K_outer <- 5
set.seed(42)
options(contrasts = c("contr.treatment", "contr.poly"))
if (!is.data.table(pheno2)) setDT(pheno2)

# output dirs
stopifnot(exists("output"))
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# helpers
predictor_pool <- setdiff(names(pheno2), y_cols)

get_xcols_keeponly <- function(tag) {
  if (tag == "ALL") {
    return(predictor_pool)
  } else {
    keep <- intersect(variants[[tag]], predictor_pool)
    if (length(keep) == 0) {
      stop(sprintf("[%s] Keep-only set is empty after intersecting with data. Check names.", tag))
    }
    return(keep)
  }
}

rse_safe <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]
  tss <- sum((obs - mean(obs))^2)
  if (length(obs) == 0 || tss == 0) return(NA_real_)
  sum((obs - pred)^2) / tss
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
    OVER[i]  <- sum(ee[ee > 0])
    UNDER[i] <- sum(-ee[ee < 0])
  }
  data.frame(OVER = OVER, UNDER = UNDER)
}

rroc_aoc <- function(rroc_df) {
  d <- rroc_df[order(rroc_df$OVER), ]
  if (nrow(d) < 2) return(NA_real_)
  sum((d$UNDER[-1] + d$UNDER[-nrow(d)]) * diff(d$OVER) / 2)
}

clean_for_rfsrc <- function(df, y_cols) {
  df <- as.data.frame(df)
  # outcomes numeric
  for (yn in y_cols) {
    if (!is.numeric(df[[yn]])) suppressWarnings(df[[yn]] <- as.numeric(df[[yn]]))
  }
  # predictors: char/logical -> factor; dates -> numeric
  preds <- setdiff(names(df), y_cols)
  for (nm in preds) {
    x <- df[[nm]]
    if (is.character(x)) df[[nm]] <- factor(x)
    else if (is.logical(x)) df[[nm]] <- factor(x)
    else if (inherits(x, "Date")) df[[nm]] <- as.numeric(x)
    else if (inherits(x, c("POSIXct","POSIXt"))) df[[nm]] <- as.numeric(x)
  }
  # guard against Inf/NaN in numeric
  for (nm in names(df)) {
    v <- df[[nm]]
    if (is.numeric(v)) {
      bad <- !is.finite(v)
      if (any(bad)) df[[nm]][bad] <- NA_real_
    }
  }
  if ("SampleID" %in% names(df)) df$SampleID <- NULL
  df
}

# accumulators
rows_outcome <- list()
rows_variant <- list()

# variant loop (joinet-style OOF for mRF)
for (tag in names(variants)) {
  message("\n==== Variant: ", tag, " ====")
  
  x_cols <- get_xcols_keeponly(tag)
  sel <- c(y_cols, x_cols)
  if (!all(sel %in% names(pheno2))) {
    miss <- setdiff(sel, names(pheno2))
    stop(sprintf("[%s] Missing required columns in pheno2: %s", tag, paste(miss, collapse=", ")))
  }
  
  df <- as.data.frame(pheno2)[, sel, drop = FALSE]
  keep <- stats::complete.cases(df)
  if (sum(!keep) > 0) message(sprintf("[%s] Dropped %d rows with NA in used columns", tag, sum(!keep)))
  df <- df[keep, , drop = FALSE]
  if (nrow(df) < max(K_outer, 10)) stop(sprintf("[%s] Too few rows after NA filtering (%d).", tag, nrow(df)))
  
  # clean for rfsrc and set up
  df <- clean_for_rfsrc(df, y_cols)
  Y <- as.matrix(df[, y_cols, drop = FALSE])
  n <- nrow(df); m <- ncol(Y)
  
  # build formula
  rhs <- paste(x_cols, collapse = " + ")
  fml <- stats::as.formula(paste0("Multivar(", paste(y_cols, collapse = ", "), ") ~ ", rhs))
  
  # outer OOF folds
  set.seed(42)
  foldid <- sample(rep(1:K_outer, length.out = n))
  oof_pred <- matrix(NA_real_, nrow = n, ncol = m, dimnames = list(NULL, y_cols))
  
  for (f in sort(unique(foldid))) {
    te <- which(foldid == f)
    tr <- setdiff(seq_len(n), te)
    
    fit <- randomForestSRC::rfsrc(
      formula    = fml,
      data       = df[tr, , drop = FALSE],
      importance = TRUE,
      forest     = TRUE,
      na.action  = "na.omit"
    )
    
    pr <- predict(fit, newdata = df[te, , drop = FALSE])
    # try to extract multivariate predictions consistently
    if (!is.null(pr$regrOutput)) {
      # regrOutput is a list per outcome with $predicted
      yhat_mat <- sapply(y_cols, function(nn) as.numeric(pr$regrOutput[[nn]]$predicted))
    } else if (is.matrix(pr$predicted) && all(y_cols %in% colnames(pr$predicted))) {
      yhat_mat <- as.matrix(pr$predicted[, y_cols, drop = FALSE])
    } else {
      stop(sprintf("[%s] Could not parse predict(rfsrc) output.", tag))
    }
    
    if (!all(dim(yhat_mat) == c(length(te), m))) {
      stop(sprintf("[%s] Unexpected predict(rfsrc) shape: got %s, expected %d x %d",
                   tag, paste(dim(yhat_mat), collapse="x"), length(te), m))
    }
    oof_pred[te, ] <- yhat_mat
  }
  
  if (any(!is.finite(rowMeans(oof_pred, na.rm = TRUE)))) {
    stop(sprintf("[%s] Missing OOF predictions for some rows.", tag))
  }
  
  # metrics per outcome (OOF)
  rmse_j <- sqrt(colMeans((oof_pred - Y)^2, na.rm = TRUE))
  rse_j  <- vapply(seq_len(m), function(j) rse_safe(Y[, j], oof_pred[, j]), numeric(1))
  r2_j   <- 1 - rse_j
  mae_j  <- colMeans(abs(oof_pred - Y), na.rm = TRUE)
  mape_j <- {
    den  <- abs(Y)
    num  <- abs(oof_pred - Y)
    mask <- is.finite(den) & is.finite(num) & (den > 1e-8)
    colSums((num / den) * mask, na.rm = TRUE) / pmax(colSums(mask), 1) * 100
  }
  aoc_j  <- vapply(seq_len(m), function(j) {
    dr <- rroc_curve(Y[, j], oof_pred[, j])
    rroc_aoc(dr)
  }, numeric(1))
  
  per_out <- data.table(
    variant    = tag,
    outcome    = y_cols,
    rmse_oof   = as.numeric(rmse_j),
    r2_oof     = as.numeric(r2_j),
    rse_oof    = as.numeric(rse_j),
    aoc_rroc   = as.numeric(aoc_j),
    mae        = as.numeric(mae_j),
    mape       = as.numeric(mape_j),
    n_features = length(x_cols),
    n          = n
  )
  fwrite(per_out, file.path(output, sprintf("mrf_outcomes_%s.csv", tag)))
  
  # save OOF predictions (optional)
  oof_dt <- as.data.table(oof_pred); setnames(oof_dt, y_cols); oof_dt[, .row := seq_len(.N)]
  fwrite(oof_dt, file.path(output, sprintf("mrf_oof_predictions_%s.csv", tag)))
  
  rows_outcome[[length(rows_outcome)+1]] <- per_out
  rows_variant[[length(rows_variant)+1]] <- data.table(
    variant    = tag,
    macro_rmse = mean(rmse_j, na.rm = TRUE),
    macro_r2   = mean(r2_j,   na.rm = TRUE),
    macro_rse  = mean(rse_j,  na.rm = TRUE),
    macro_aoc  = mean(aoc_j,  na.rm = TRUE),
    macro_mae  = mean(mae_j,  na.rm = TRUE),
    macro_mape = mean(mape_j, na.rm = TRUE),
    n_features = length(x_cols),
    n          = n
  )
  
  # one RROC PDF per variant
  rroc_pdf <- file.path(output, sprintf("mrf_rroc_%s.pdf", tag))
  grDevices::pdf(rroc_pdf, width = 7, height = 6)
  for (j in seq_along(y_cols)) {
    y  <- as.numeric(Y[, j])
    yh <- as.numeric(oof_pred[, j])
    expl <- DALEX::explain(
      model = list(),
      data  = data.frame(.oof_id = seq_along(yh)),
      y     = y,
      predict_function = function(m, newdata) yh[newdata$.oof_id],
      label = sprintf("%s | %s (mRF OOF)", y_cols[j], tag),
      verbose = FALSE
    )
    resid_obj <- auditor::model_residual(expl)
    p <- auditor::plot_rroc(resid_obj) + ggplot2::ggtitle(sprintf("RROC — %s | %s", y_cols[j], tag))
    print(p)
  }
  grDevices::dev.off()
  message(sprintf("[%s] Saved RROC PDF: %s", tag, normalizePath(rroc_pdf)))
}

# combine & save master tables
by_outcome <- rbindlist(rows_outcome, use.names = TRUE, fill = TRUE)
by_variant <- rbindlist(rows_variant, use.names = TRUE, fill = TRUE)

# sanity: recompute macro from per-outcome
macro_check <- by_outcome[, .(
  macro_rmse = mean(rmse_oof, na.rm = TRUE),
  macro_r2   = mean(r2_oof,   na.rm = TRUE),
  macro_rse  = mean(rse_oof,  na.rm = TRUE),
  macro_aoc  = mean(aoc_rroc, na.rm = TRUE),
  macro_mae  = mean(mae,      na.rm = TRUE),
  macro_mape = mean(mape,     na.rm = TRUE)
), by = variant]

setkey(macro_check, variant); setkey(by_variant, variant)
merged <- merge(by_variant, macro_check, by = "variant", suffixes = c(".x",".y"))
# check equality (tolerant)
ok_macro <- isTRUE(all.equal(merged$macro_rmse.x, merged$macro_rmse.y, tolerance = 1e-12)) &&
  isTRUE(all.equal(merged$macro_r2.x,   merged$macro_r2.y,   tolerance = 1e-12)) &&
  isTRUE(all.equal(merged$macro_rse.x,  merged$macro_rse.y,  tolerance = 1e-12))
if (!ok_macro) warning("Macro metrics differ between assembled and recomputed versions (rounding expected).")

fwrite(by_outcome, file.path(output, "INT_mrf_outcomes_by_variant.csv"))
fwrite(by_variant, file.path(output, "INT_mrf_macro_by_variant.csv"))

message("\n== Macro by variant (sorted by RMSE) ==")
print(by_variant[order(macro_rmse)])

message("\nDone. Outputs saved in: ", normalizePath(output))


### Running model with ideal conditions

rhs <- paste(x_cols, collapse = " + ")
fml <- as.formula(paste0("Multivar(", paste(y_cols, collapse = ", "), ") ~ ", rhs)) 

set.seed(42)
fit_final <- randomForestSRC::rfsrc(
  formula    = fml,
  data       = pheno2,
  importance = TRUE,
  forest     = TRUE,
  na.action  = "na.omit"
)

saveRDS(fit_final, file.path(output, "final_mrf_model.rds"))

#prepping data 

a_test <- a

ids <- a_test[order(DS1V)][, head(.SD, 1), by = "PIDN"][, SampleID]
f_ups <- a_test[!SampleID %in% ids]
columns_to_check <- c(  "Age","Gender","HybControlNormScale", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","APOE.e34","APOE.e44",
                        "PlasmaPTau181","PlasmaAB142P","PlasmaAB140P","PlasmaABRatio","PlasmapTau217","PlasmapTau217_AB42Ratio",
                        "PlasmaGFAP","PlasmaNfL")
drops <- f_ups[!complete.cases(f_ups[, ..columns_to_check ]), SampleID]
a_test <- f_ups[!SampleID %in% drops]

test_data <- data.frame(a_test[, ..cols], row.names = "SampleID")
test_data <- as.matrix(test_data[complete.cases(test_data),]) 

actual_values <- as.data.frame(test_data[, c("APP.Metabolism", "Apoptosis",	"Autophagy", "Cell.Cycle",	"DNA.Repair",	"Endolysosome",	"Epigenetic",	
                                             "Immune.Response",	"Lipid.Metabolism",	"Metal.Binding.and.Homeostasis",	"Mitochondrial.Metabolism",	"Myelination",	
                                             "Oxidative.Stress",	"Proteostasis",	"RNA.Spliceosome",	"Structural.Stabilization",	"Synapse",	"Tau.Homeostasis",	"Vasculature")])
test_data <- pheno2_norm
test_data <- test_data[, !colnames(test_data) %in% c("APP.Metabolism", "Apoptosis",	"Autophagy", "Cell.Cycle",	"DNA.Repair",	"Endolysosome",	"Epigenetic",	
                                                     "Immune.Response",	"Lipid.Metabolism",	"Metal.Binding.and.Homeostasis",	"Mitochondrial.Metabolism",	"Myelination",	
                                                     "Oxidative.Stress",	"Proteostasis",	"RNA.Spliceosome",	"Structural.Stabilization",	"Synapse",	"Tau.Homeostasis",	"Vasculature")]

#prediction- RANDOM FOREST USES NEWDATA!!
rf_pred <- predict(fit_final,newdata =test_data)

#extract predictions 

cell_types <- c("APP.Metabolism", "Apoptosis",	"Autophagy", "Cell.Cycle",	"DNA.Repair",	"Endolysosome",	"Epigenetic",	
                "Immune.Response",	"Lipid.Metabolism",	"Metal.Binding.and.Homeostasis",	"Mitochondrial.Metabolism",	"Myelination",	
                "Oxidative.Stress",	"Proteostasis",	"RNA.Spliceosome",	"Structural.Stabilization",	"Synapse",	"Tau.Homeostasis",	"Vasculature")
if (!is.null(rf_pred$regrOutput)) {
  pred_mat <- sapply(cell_types, function(y) as.numeric(rf_pred$regrOutput[[y]]$predicted))
  pred_meta <- as.data.frame(pred_mat)
  colnames(pred_meta) <- cell_types
} else if (is.matrix(rf_pred$predicted) && all(cell_types %in% colnames(rf_pred$predicted))) {
  pred_meta <- as.data.frame(rf_pred$predicted[, cell_types, drop = FALSE])
} else {
  stop("Could not parse predict(rfsrc) output for multivariate regression.")
if (!is.null(rownames(test_data))) {
    rownames(pred_meta) <- rownames(test_data)
  } else if ("SampleID" %in% colnames(test_data)) {
    rownames(pred_meta) <- as.character(test_data[["SampleID"]])
  } else {
    rownames(pred_meta) <- sprintf("row_%d", seq_len(nrow(pred_meta)))
  }
}

## Row/col names (mirror your joinet block)
if (!is.null(rownames(test_data))) {
  rownames(pred_meta) <- rownames(test_data)
} else if ("SampleID" %in% colnames(test_data)) {
  rownames(pred_meta) <- as.character(test_data[["SampleID"]])
} else {
  rownames(pred_meta) <- sprintf("row_%d", seq_len(nrow(pred_meta)))
  
#Sanity checks against actual_values (must share rows and columns)
stopifnot(all(cell_types %in% colnames(actual_values)))

# enforce same outcome order
actual_values <- actual_values[, cell_types, drop = FALSE]
  
identical(colnames(pred_meta), colnames(actual_values))

#Get stats 
common_rows <- intersect(rownames(actual_values), rownames(pred_meta))
  A <- actual_values[common_rows, , drop = FALSE]
  P <- pred_meta[common_rows,  , drop = FALSE]
  
metrics <- data.frame(
    Outcome = colnames(A),
    RMSE = NA_real_,
    MAE  = NA_real_,
    MAPE = NA_real_,
    R2   = NA_real_
  )
  
for (i in seq_along(metrics$Outcome)) {
    y_true <- as.numeric(A[, i])
    y_pred <- as.numeric(P[, i])
    
    rmse <- sqrt(mean((y_true - y_pred)^2, na.rm = TRUE))
    mae  <- mean(abs(y_true - y_pred), na.rm = TRUE)
    den  <- pmax(abs(y_true), 1e-8)                              # robust MAPE denom
    mape <- mean(abs(y_true - y_pred) / den, na.rm = TRUE) * 100
    # safer R² (can be negative if worse than mean)
    sse <- sum((y_true - y_pred)^2, na.rm = TRUE)
    tss <- sum((y_true - mean(y_true, na.rm = TRUE))^2, na.rm = TRUE)
    r2  <- if (tss == 0) NA_real_ else 1 - sse / tss
    
    metrics$RMSE[i] <- rmse
    metrics$MAE[i]  <- mae
    metrics$MAPE[i] <- mape
    metrics$R2[i]   <- r2
  }
  
print(metrics)

actual_plus_pred <- cbind(
  A,
  setNames(P, paste0(colnames(P), "_pred"))
)

## 7) Predicted vs Observed plots (one PDF per outcome)
cells <- cell_types
merged_data <- cbind(
  data.frame(SampleID = common_rows, stringsAsFactors = FALSE),
  A[common_rows, , drop = FALSE],
  setNames(P[common_rows, , drop = FALSE], paste0(colnames(P), "_pred"))
)


for (cell in cells) {
  xcol <- cell                      # Observed
  ycol <- paste0(cell, "_pred")     # Predicted
  
  p <- ggplot(merged_data, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_point(alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    labs(title = paste("Predicted vs Observed —", cell),
         x = "Observed", y = "Predicted") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(file.path(output, paste0("rf_multi_predicted_vs_observed_", cell, ".pdf")),
         plot = p, width = 7, height = 7, dpi = 300)
}



