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

input <- "C:/Users/audweber/Box/plasma_ptau_proteomics/dataProt_SS-205063.merged.QC.csv" 
output <- "prediction_models"
prots <- "C:/Users/audweber/Box/plasma_ptau_proteomics/prot_names.txt" 
covs <- 'Age,Gender,APOE'   
outcome <- c("Microglia","Astrocyte","Neuron","Oligodendrocyte") 
nameset <- "dataProt_SS-205063.PC.Plasma_pTau217.APOE" 
pcs <- "C:/Users/audweber/Box/plasma_ptau_proteomics/genomic/mypcair.PCs.csv"
prediction_scores <- "C:/Users/audweber/Box/plasma_ptau_proteomics/celltype_scores/dataProt_SS-205063.Sharma.lvl1.normalized_scores.csv"
feature <- c("PlasmaPTau181", "PlasmaAB142P","PlasmaAB140P","PlasmaABRatio", "PlasmapTau217","PlasmapTau217_AB42Ratio","PlasmaGFAP","PlasmaNfL")

covs <- strsplit(covs,",")[[1]]
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
pheno2$Gender <- as.numeric(pheno2$Gender == "M")
pheno2$APOE <- ifelse(grepl("4", pheno2$APOE), 1, 0) 



## outcomes- cell types 

y_cols <- c("Microglia","Astrocyte","Neuron","Oligodendrocyte")  

#predictors

variants <- list(
  ALL             = character(0),              
  NO_PCSandAPOE   = c("APOE","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"),
  Only_biomarkers = c("PlasmaPTau181","PlasmaAB142P","PlasmaAB140P","PlasmaABRatio",
                      "PlasmapTau217","PlasmapTau217_AB42Ratio","PlasmaGFAP","PlasmaNfL")
)

predictor_pool <- setdiff(names(pheno2), y_cols)

tuning_grid <- expand.grid(
  alpha = c(0, 0.1, 0.2, 0.3, 0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
  lambda = seq(0.001, 0.1, length = 20)  # Test 10 lambda values
)


# Safe RSE (SSE/TSS), NA if TSS == 0
rse_safe <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]
  sse <- sum((obs - pred)^2)
  tss <- sum((obs - mean(obs))^2)
  if (tss == 0) return(NA_real_)
  sse / tss
}

# Simple RROC + AOC (for a scalar summary alongside the auditor plot)
# OVER on x-axis, UNDER on y-axis
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
  d <- rroc_df[order(rroc_df$OVER), ]
  if (nrow(d) < 2) return(NA_real_)
  sum((d$UNDER[-1] + d$UNDER[-nrow(d)]) * diff(d$OVER) / 2)
}

## ===== data slices (adjust if needed) =====
y_idx <- 1:4
x_idx <- 5:23

Y <- as.matrix(pheno2[, y_idx, drop = FALSE])
X <- as.matrix(pheno2[, x_idx, drop = FALSE])
if (is.null(colnames(Y))) colnames(Y) <- paste0("y", seq_len(ncol(Y)))
if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))

outcome_names <- colnames(Y)
n <- nrow(Y); m <- ncol(Y)

## ===== 1) pick alphas via cv.joinet (macro MSE on meta) =====
set.seed(42)
foldid.ext <- sample(rep(1:5, length.out = n))   # outer folds (selection)
foldid.int <- sample(rep(1:10, length.out = n))  # inner folds (glmnet tuning)

alphas <- seq(0, 1, 0.1)
grid <- CJ(alpha.base = alphas, alpha.meta = alphas)
res_list <- vector("list", nrow(grid))

for (k in seq_len(nrow(grid))) {
  ab <- grid$alpha.base[k]; am <- grid$alpha.meta[k]
  mcv <- cv.joinet(
    Y = Y, X = X, family = "gaussian", type.measure = "mse",
    alpha.base = ab, alpha.meta = am,
    nfolds.ext = 5, nfolds.int = 10,
    foldid.ext = foldid.ext, foldid.int = foldid.int
  )
  res_list[[k]] <- data.table(
    alpha.base = ab,
    alpha.meta = am,
    method     = rownames(mcv),          # "none","base","meta"
    mse        = rowMeans(mcv)           # macro-average across outcomes
  )
}

sums_wide <- dcast(rbindlist(res_list), alpha.base + alpha.meta ~ method, value.var = "mse")
sums_wide[, `:=`(
  diff      = base - meta,               # >0 means META better than BASE
  base.rmse = sqrt(base),
  meta.rmse = sqrt(meta),
  none.rmse = sqrt(none)
)]
setorder(sums_wide, meta)                 # choose alphas minimizing meta MSE
best_ab <- sums_wide$alpha.base[1]
best_am <- sums_wide$alpha.meta[1]
fwrite(sums_wide, file.path(output, "joinet_alpha_grid.csv"))
message(sprintf("Chosen alphas: alpha.base=%.2f, alpha.meta=%.2f", best_ab, best_am))

## ===== 2) build true OOF predictions for metrics/plots =====
set.seed(42)
foldid <- sample(rep(1:5, length.out = n))  # final 5-fold OOF

oof_pred <- matrix(NA_real_, nrow = n, ncol = m, dimnames = list(NULL, outcome_names))

for (f in sort(unique(foldid))) {
  test_idx  <- which(foldid == f)
  train_idx <- setdiff(seq_len(n), test_idx)
  
  fit <- joinet(
    Y = Y[train_idx, , drop = FALSE],
    X = X[train_idx, , drop = FALSE],
    family = "gaussian",
    nfolds = 10,                     # inner CV for lambdas
    type.measure = "mse",
    alpha.base = best_ab,
    alpha.meta = best_am,
    compare = TRUE
  )
  
  pred <- try(predict(fit, newx = X[test_idx, , drop = FALSE]), silent = TRUE)
  if (inherits(pred, "try-error")) stop("predict(joinet, newx=...) failed; check joinet version.")
  # Typical shape is a list with $meta (n_test × m); fallback to matrix if needed
  if (is.list(pred) && !is.null(pred$meta)) {
    yhat <- as.matrix(pred$meta)
  } else if (is.matrix(pred) && ncol(pred) == m) {
    yhat <- pred
  } else if (is.array(pred) && length(dim(pred)) >= 2) {
    yhat <- as.matrix(pred[, seq_len(m)])
  } else {
    stop("Unrecognized predict(joinet) output shape; run str(pred).")
  }
  
  oof_pred[test_idx, ] <- yhat
}

## ===== 3) per-outcome metrics + auditor-style RROC =====
metrics <- data.table(
  outcome = outcome_names,
  RMSE    = NA_real_,
  RSE     = NA_real_,
  R2      = NA_real_,
  AOC     = NA_real_
)

rroc_pdf <- file.path(output, "joinet_rroc_auditor.pdf")
pdf(rroc_pdf, width = 7, height = 6)
on.exit(try(dev.off(), silent = TRUE), add = TRUE)

for (j in seq_along(outcome_names)) {
  y  <- as.numeric(Y[, j])
  yh <- as.numeric(oof_pred[, j])
  
  rmse_val <- sqrt(mean((y - yh)^2, na.rm = TRUE))
  rse_val  <- rse_safe(y, yh)
  r2_val   <- if (is.na(rse_val)) NA_real_ else (1 - rse_val)
  
  # scalar AOC (our own calc; auditor doesn't export it directly)
  dr   <- rroc_curve(y, yh, n_points = 400)
  aocv <- rroc_aoc(dr)
  
  metrics$RMSE[j] <- rmse_val
  metrics$RSE[j]  <- rse_val
  metrics$R2[j]   <- r2_val
  metrics$AOC[j]  <- aocv
  
  # auditor-style RROC using OOF preds (one page per outcome)
  expl <- DALEX::explain(
    model = list(),  # dummy; we override predict to return OOF preds
    data  = data.frame(.oof_id = seq_along(yh)),
    y     = y,
    predict_function = function(m, newdata) yh[newdata$.oof_id],
    label = sprintf("%s (joinet meta)", outcome_names[j]),
    verbose = FALSE
  )
  resid_obj <- auditor::model_residual(expl)
  p <- auditor::plot_rroc(resid_obj)     # classic auditor look
  print(p)
}

dev.off()

## ===== 4) save metrics + macro =====
fwrite(metrics, file.path(output, "joinet_multivariate_metrics_by_outcome.csv"))

macro <- data.table(
  macro_RMSE = mean(metrics$RMSE, na.rm = TRUE),
  macro_RSE  = mean(metrics$RSE,  na.rm = TRUE),
  macro_R2   = mean(metrics$R2,   na.rm = TRUE),
  macro_AOC  = mean(metrics$AOC,  na.rm = TRUE),
  alpha_base = best_ab,
  alpha_meta = best_am
)
fwrite(macro, file.path(output, "joinet_multivariate_macro.csv"))

cat("\nSaved:\n",
    "- Per-outcome metrics: ", file.path(output, "joinet_multivariate_metrics_by_outcome.csv"), "\n",
    "- Macro metrics:      ", file.path(output, "joinet_multivariate_macro.csv"), "\n",
    "- RROC PDF:           ", rroc_pdf, "\n", sep = "")