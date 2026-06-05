###### Joinet multivariate 
### NOTE: Index is important because it allows us to compare across conditions
######    by using the same folds for each run 

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

inpath <- "C:/Users/audweber/Box/plasma_ptau_proteomics" # don't repeat things you don't have to
input <- file.path(inpath, "dataProt_SS-205063.merged.QC.csv")
prots <- file.path(inpath, "prot_names.txt") 
pcs <- file.path(inpath, "genomic/mypcair.PCs.csv")
prediction_scores <- file.path(inpath, # wrap at 80 characters
                               "celltype_scores/dataProt_SS-205063.SEA-AD.Neighborhood.1kDS.lvl2.normalized_scores.csv"
)
output <- 'prediction_models/SEA_AD/joinet/INT'
covs <- 'Age,Gender,APOE,HybControlNormScale'   
outcome <- c("Astrocyte","CGE","L2/3.IT","L4/5/6.IT-ET.Car3","MGE","Microglia.PVM","NP.CT.L6b","Oligo","Vascular") 
feature <- c(
  "PlasmaPTau181", "PlasmaAB142P","PlasmaAB140P","PlasmaABRatio", 
  "PlasmapTau217","PlasmapTau217_AB42Ratio","PlasmaGFAP","PlasmaNfL"
)
qc_cols <- c(
  "PIDN", "SampleID", "ADRC_ID", "Visit", "Visit_Consensus", "ConnectivityZscore",
  "Diagnosis_group", "Date.of.draw", "DS1V"
)

# Main -----------------------
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


cols = c(qc_cols, covs)

drop = a[abs(ConnectivityZscore) > 10, SampleID]
a = a[!SampleID %in% drop]


# Setup for baseline visit only analysis ---------------------------
# Correction, keep only first visit where plasma available
a = a[order(DS1V)][, head(.SD, 1), by = "PIDN"]

# checking for missed duplicates or bad visit labeling
a[duplicated(PIDN), ..qc_cols]


# handling incomplete cases
cols <- c("SampleID", outcome, covs, feature)
drops <- a[!complete.cases(a[, ..cols]), SampleID]
a <- a[!SampleID %in% drops]

#Setting up data for input 
pheno2 <- data.frame(a[, ..cols], row.names = "SampleID")
pheno2 <- pheno2[complete.cases(pheno2),] 


#All variables need to be numeric 
a[, Gender := as.numeric(Gender == "M")]
recode_dummies(a, "APOE", keep = c("e34", "e44"))
covs <- unique(c(covs[covs != "APOE"], "APOE.e34", "APOE.e44"))

## outcomes- cell types 

y_cols <- outcome 

#predictors

variants <- list(
  ALL             = character(0),              
  NO_PCSandAPOE   = c( "Age", "Gender", "PlasmaPTau181_norm","PlasmaAB142P_norm","PlasmaAB140P_norm","PlasmaABRatio_norm",
                       "PlasmapTau217_norm","PlasmapTau217_AB42Ratio_norm","PlasmaGFAP_norm","PlasmaNfL_norm", "HybControlNormScale"),
  NO_PCS = c("Age", "Gender", "APOE.e34", "APOE.e44", "PlasmaPTau181","PlasmaAB142P","PlasmaAB140P","PlasmaABRatio",
                      "PlasmapTau217_norm","PlasmapTau217_AB42Ratio_norm","PlasmaGFAP_norm","PlasmaNfL_norm", "HybControlNormScale"),
  Impt_biomarkers = c("PlasmaPTau181_norm","PlasmaABRatio_norm","PlasmaGFAP_norm","PlasmaNfL_norm", "Age", "Gender", "APOE.e34", "APOE.e44","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8", "HybControlNormScale")
)

predictor_pool <- setdiff(names(pheno2), y_cols)

tuning_grid <- expand.grid(
  alpha = c(0, 0.1, 0.2, 0.3, 0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
  lambda = seq(0.001, 0.1, length = 20)  # Test 10 lambda values
)


# CV settings
K_outer <- 5   # OOF folds
K_inner <- 10  # inner folds for lambdas in joinet
alpha_grid <- seq(0, 1, 0.1)

# Reproducibility
set.seed(42)
options(contrasts = c("contr.treatment", "contr.poly"))
if (!is.data.table(pheno2)) setDT(pheno2)

## ---- Helpers ----
predictor_pool <- setdiff(names(pheno2), y_cols)

get_xcols <- function(tag) {
  if (tag == "ALL") {
    # use everything that isn't an outcome
    return(predictor_pool)
  } else {
    # KEEP-ONLY: intersect requested keeps with what exists
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

rroc_aoc <- function(rroc_df) {
  d <- rroc_df[order(rroc_df$OVER), ]
  if (nrow(d) < 2) return(NA_real_)
  sum((d$UNDER[-1] + d$UNDER[-nrow(d)]) * diff(d$OVER) / 2)
}


## ---- Accumulators ----
rows_outcome <- list()  # per outcome per variant
rows_variant <- list()  # macro per variant

## ---- Variant loop ----
for (tag in names(variants)) {
  message("\n==== Variant: ", tag, " ====")
  
  # 1) Build data for this variant
  x_cols <- get_xcols(tag)
  if (length(x_cols) == 0) {
    warning(sprintf("[%s] No predictors selected; skipping", tag))
    next
  }
  
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
  
  Y <- as.matrix(df[, y_cols, drop = FALSE])
  X <- as.matrix(df[, x_cols,   drop = FALSE])
  n <- nrow(Y); m <- ncol(Y)
  stopifnot(nrow(X) == n, ncol(Y) == length(y_cols), identical(colnames(Y), y_cols))
  
  # 2) Alpha selection via cv.joinet (saved to CSV)
  set.seed(42)
  foldid.ext <- sample(rep(1:K_outer, length.out = n))
  foldid.int <- sample(rep(1:K_inner, length.out = n))
  grid <- CJ(alpha.base = alpha_grid, alpha.meta = alpha_grid)
  scores <- vector("list", nrow(grid))
  
  for (k in seq_len(nrow(grid))) {
    ab <- grid$alpha.base[k]; am <- grid$alpha.meta[k]
    mcv <- cv.joinet(
      Y = Y, X = X, family = "gaussian",
      type.measure = "mse",
      alpha.base = ab, alpha.meta = am,
      nfolds.ext = K_outer, nfolds.int = K_inner,
      foldid.ext = foldid.ext, foldid.int = foldid.int
    )
    # mcv rows: "none","base","meta"; cols: outcomes
    scores[[k]] <- data.table(alpha.base=ab, alpha.meta=am,
                              method=rownames(mcv), mse=rowMeans(mcv))
  }
  
  sw <- dcast(rbindlist(scores), alpha.base + alpha.meta ~ method, value.var = "mse")
  if (!all(c("none","base","meta") %in% names(sw))) stop(sprintf("[%s] cv.joinet output missing methods.", tag))
  setorderv(sw, "meta")
  best_ab <- sw$alpha.base[1]
  best_am <- sw$alpha.meta[1]
  message(sprintf("[%s] Best alphas: alpha.base=%.2f  alpha.meta=%.2f", tag, best_ab, best_am))
  
  # Save the alpha grid for this variant
  fwrite(sw[, .(alpha.base, alpha.meta, none, base, meta)],
         file.path(output, sprintf("joinet_alpha_grid_%s.csv", tag)))
  
  # 3) True OOF predictions (K_outer folds) at chosen alphas
  set.seed(42)
  foldid <- sample(rep(1:K_outer, length.out = n))
  oof_pred <- matrix(NA_real_, nrow = n, ncol = m, dimnames = list(NULL, colnames(Y)))
  
  for (f in sort(unique(foldid))) {
    te <- which(foldid == f)
    tr <- setdiff(seq_len(n), te)
    
    fit <- joinet(
      Y = Y[tr, , drop = FALSE],
      X = X[tr, , drop = FALSE],
      family = "gaussian",
      nfolds = K_inner,
      type.measure = "mse",
      alpha.base = best_ab,
      alpha.meta = best_am,
      compare = TRUE
    )
    
    pred <- predict(fit, newx = X[te, , drop = FALSE])
    # handle different shapes from joinet versions
    yhat <- if (is.list(pred) && !is.null(pred$meta)) as.matrix(pred$meta) else as.matrix(pred)
    if (!all(dim(yhat) == c(length(te), m))) {
      stop(sprintf("[%s] Unexpected predict(joinet) shape: got %s, expected %d x %d",
                   tag, paste(dim(yhat), collapse="x"), length(te), m))
    }
    oof_pred[te, ] <- yhat
  }
  
  # Sanity: OOF coverage
  if (any(!is.finite(rowMeans(oof_pred, na.rm = TRUE))))
    stop(sprintf("[%s] Missing OOF predictions for some rows.", tag))
  
  # 4) Metrics per outcome
  rmse_j <- sqrt(colMeans((oof_pred - Y)^2, na.rm = TRUE))
  rse_j  <- vapply(seq_len(m), function(j) rse_safe(Y[, j], oof_pred[, j]), numeric(1))
  r2_j   <- 1 - rse_j
  mae_j <- colMeans(abs(oof_pred - Y), na.rm = TRUE)
  mape_j <- {
    den  <- abs(Y)
    num  <- abs(oof_pred - Y)
    mask <- is.finite(den) & is.finite(num) & (den > 1e-8)
    colSums((num / den) * mask, na.rm = TRUE) / pmax(colSums(mask), 1) * 100
  }
  stopifnot(all(is.finite(rmse_j)))
  
  aoc_j  <- vapply(seq_len(m), function(j) {
    dr <- rroc_curve(Y[, j], oof_pred[, j])
    rroc_aoc(dr)
  }, numeric(1))
  
  # Save per-outcome metrics for this variant
  per_out <- data.table(
    variant    = tag,
    outcome    = colnames(Y),
    rmse_oof   = as.numeric(rmse_j),
    r2_oof     = as.numeric(r2_j),
    rse_oof    = as.numeric(rse_j),
    aoc_rroc   = as.numeric(aoc_j), 
    mae        = as.numeric(mae_j),
    mape       = as.numeric(mape_j),
    n_features = ncol(X),
    alpha_base = best_ab,
    alpha_meta = best_am,
    n          = n
  )
  fwrite(per_out, file.path(output, sprintf("joinet_outcomes_%s.csv", tag)))
  
  # Save OOF predictions (optional but useful)
  oof_dt <- as.data.table(oof_pred)
  setnames(oof_dt, colnames(Y))
  oof_dt[, .row := seq_len(.N)]
  fwrite(oof_dt, file.path(output, sprintf("joinet_oof_predictions_%s.csv", tag)))
  
  # Append to accumulators
  rows_outcome[[length(rows_outcome)+1]] <- per_out
  rows_variant[[length(rows_variant)+1]] <- data.table(
    variant    = tag,
    macro_rmse = mean(rmse_j, na.rm = TRUE),
    macro_r2   = mean(r2_j,   na.rm = TRUE),
    macro_rse  = mean(rse_j,  na.rm = TRUE),
    macro_aoc  = mean(aoc_j,  na.rm = TRUE),
    macro_mae  = mean(mae_j, na.rm = TRUE), 
    macro_mape = mean(mape_j, na.rm = TRUE),
    alpha_base = best_ab,
    alpha_meta = best_am,
    n_features = ncol(X),
    n          = n
  )
  
  # 5) RROC PDF (one file per variant; one page per outcome)
  rroc_pdf <- file.path(output, sprintf("joinet_rroc_%s.pdf", tag))
  pdf(rroc_pdf, width = 7, height = 6)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)  # ensure close even on error
  
  for (j in seq_len(ncol(Y))) {
    y  <- as.numeric(Y[, j])
    yh <- as.numeric(oof_pred[, j])
    
    expl <- DALEX::explain(
      model = list(),  # dummy model; we feed OOF preds via predict_function
      data  = data.frame(.oof_id = seq_along(yh)),
      y     = y,
      predict_function = function(m, newdata) yh[newdata$.oof_id],
      label = sprintf("%s | %s (Joinet meta)", colnames(Y)[j], tag),
      verbose = FALSE
    )
    resid_obj <- auditor::model_residual(expl)
    p <- auditor::plot_rroc(resid_obj)
    print(p)
  }
  dev.off()
  message(sprintf("[%s] Saved RROC PDF: %s", tag, normalizePath(rroc_pdf)))
}

## ---- Combine & save master tables ----
by_outcome <- rbindlist(rows_outcome, use.names = TRUE, fill = TRUE)
by_variant <- rbindlist(rows_variant, use.names = TRUE, fill = TRUE)

# Recompute macro to sanity-check
macro_check <- by_outcome[, .(
  macro_rmse = mean(rmse_oof, na.rm = TRUE),
  macro_r2   = mean(r2_oof,   na.rm = TRUE),
  macro_rse  = mean(rse_oof,  na.rm = TRUE)
), by = variant]
setkey(macro_check, variant); setkey(by_variant, variant)
merged <- merge(by_variant[, .(variant, macro_rmse, macro_r2, macro_rse)],
                macro_check, by = "variant", suffixes = c(".x",".y"))
ok_macro <- isTRUE(all.equal(merged$macro_rmse.x, merged$macro_rmse.y, tolerance = 1e-12)) &&
  isTRUE(all.equal(merged$macro_r2.x,   merged$macro_r2.y,   tolerance = 1e-12)) &&
  isTRUE(all.equal(merged$macro_rse.x,  merged$macro_rse.y,  tolerance = 1e-12))
if (!ok_macro) stop("Macro metrics mismatch between assembled and recomputed versions.")

fwrite(by_outcome, file.path(output, "joinet_outcomes_by_variant.csv"))
fwrite(by_variant, file.path(output, "joinet_macro_by_variant.csv"))

message("\n== Saved files ==")
print(list.files(output, full.names = TRUE))

message("\n== Macro by variant (sorted by RMSE) ==")
print(by_variant[order(macro_rmse)])

message("\nDone. All outputs saved to: ", normalizePath(output))


### Running model with ideal conditions- INT and no PC or APOE
remove <- c("Astrocyte","CGE","L2.3.IT","L4.5.6.IT.ET.Car3","MGE","Microglia.PVM","NP.CT.L6b","Oligo","Vascular","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8", "APOE.e34","APOE.e44")
cell <- pheno2[, ..y_cols, drop = FALSE]
variant <- pheno2[, .SD, .SDcols = !names(pheno2) %in% remove]


rownames(cell) == rownames(variant)    

fit <-joinet(Y=cell,X=variant,
                family = "gaussian",
                type.measure = "mse",
                alpha.base = 0.5,
                alpha.meta = 1)

#Rerunning so names are the same for y  

original_names <- colnames(cell)

# 2. Create a mapping between generic and original names
generic_names <- paste0("y", 1:length(original_names))
name_mapping <- setNames(original_names, generic_names)

# 3. Define the renaming function
rename_joinet_outcomes <- function(fit_obj, mapping) {
  # Create a copy to avoid modifying the original
  new_fit <- fit_obj
  
  # Rename in base component
  for (old_name in names(mapping)) {
    if (old_name %in% names(new_fit$base)) {
      new_name <- mapping[old_name]
      new_fit$base[[new_name]] <- new_fit$base[[old_name]]
      new_fit$base[[old_name]] <- NULL
    }
  }
  
  # Rename in meta component
  for (old_name in names(mapping)) {
    if (old_name %in% names(new_fit$meta)) {
      new_name <- mapping[old_name]
      new_fit$meta[[new_name]] <- new_fit$meta[[old_name]]
      new_fit$meta[[old_name]] <- NULL
    }
  }
  
  return(new_fit)
}

# 4. Now apply the renaming
fit_renamed <- rename_joinet_outcomes(fit, name_mapping)

# 5. Check the result
names(fit_renamed$base)

saveRDS(fit_renamed, file.path(output, "final_joinet_model.rds"))
#preparing test data 

a_test <- a

ids <- a_test[order(DS1V)][, head(.SD, 1), by = "PIDN"][, SampleID]
f_ups <- a_test[!SampleID %in% ids]
columns_to_check <- c(  "Age","Gender","HybControlNormScale", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8", "APOE.e34","APOE.e44",
                        "PlasmaPTau181","PlasmaAB142P","PlasmaAB140P","PlasmaABRatio","PlasmapTau217","PlasmapTau217_AB42Ratio",
                       "PlasmaGFAP","PlasmaNfL")
drops <- f_ups[!complete.cases(f_ups[, ..columns_to_check ]), SampleID]
a_test <- f_ups[!SampleID %in% drops]

#R cant read names as they are so need to go back and change them and redefine outcomes 
names(a_test) <- names(name_map <- setNames(names(a_test), make.names(names(a_test), unique = TRUE)))
outcome <- c("Astrocyte","CGE","L2.3.IT","L4.5.6.IT.ET.Car3","MGE",                    
             "Microglia.PVM","NP.CT.L6b","Oligo","Vascular") 
cols <- c("SampleID", outcome, covs, feature)

test_data <- data.frame(a_test[, ..cols], row.names = "SampleID")

test_data <- pheno2_norm[, !(names(pheno2_norm) %in% remove), drop = FALSE]


actual_values <- as.data.frame(pheno2_norm[, c("Astrocyte","CGE","L2.3.IT","L4.5.6.IT.ET.Car3","MGE",                    
                                               "Microglia.PVM","NP.CT.L6b","Oligo","Vascular")])
# test_data <- test_data[, !colnames(test_data) %in% c("Microglia", "Astrocyte", "Neuron", "Oligodendrocyte")]

#prediction 
joinetprediction <- predict(fit_renamed,newx=test_data)
P <- joinetprediction
pred_meta <- as.data.frame(P$meta)
rownames(pred_meta) <- if (!is.null(rownames(test_data))) rownames(test_data) else test_data[, "SampleID"]
colnames(pred_meta) <- y_cols

#sanity check 
identical(colnames(pred_meta), colnames(actual_values))     
identical(rownames(pred_meta), rownames(test_data)) 

# Model accuracy on validation set:  

metrics <- data.frame(
  Outcome = colnames(actual_values),
  RMSE = NA,
  MAE  = NA,
  MAPE = NA,
  R2   = NA
)

for (i in seq_along(metrics$Outcome)) {
  y_true <- actual_values[, i]
  y_pred <- pred_meta[, i]
  
  metrics$RMSE[i] <- sqrt(mean((y_true - y_pred)^2))
  metrics$MAE[i]  <- mean(abs(y_true - y_pred))
  metrics$MAPE[i] <- mean(abs((y_true - y_pred) / y_true)) * 100  # percentage
  metrics$R2[i]   <- cor(y_true, y_pred)^2
}

metrics

common_rows <- intersect(rownames(actual_values), rownames(pred_meta))
common_cols <- intersect(colnames(actual_values), colnames(pred_meta))

A <- actual_values[common_rows, common_cols, drop = FALSE]
P <- pred_meta[common_rows,  common_cols, drop = FALSE]

# add predictions next to actuals, with “…_pred” names
actual_plus_pred <- cbind(
  A,
  setNames(P, paste0(common_cols, "_pred"))
)

merged_data <- merge(pred_meta, actual_values, by = "row.names", all.x = TRUE)

#Plot actual vs. predicted outcomes 
cells <- c("Astrocyte","CGE","L2.3.IT","L4.5.6.IT.ET.Car3","MGE",                    
           "Microglia.PVM","NP.CT.L6b","Oligo","Vascular")  # edit to your set

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
  ggsave(file.path(output, paste0("predicted_vs_observed_", cell, ".pdf")), 
         plot = p, width = 7, height = 7, dpi = 300)
}

