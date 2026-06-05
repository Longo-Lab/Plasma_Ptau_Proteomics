



###### Individual glmnets 
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

#Join all predictors from pheno2 and get rid of cell types 
predictor_pool <- setdiff(names(pheno2), y_cols)

#Tuning grid that tests different alphas and lamdas
tuning_grid <- expand.grid(
  alpha = c(0, 0.1, 0.2, 0.3, 0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
  lambda = seq(0.001, 0.1, length = 20)  # Test 10 lambda values
)

#Drops values that are not finite and calculates see and TSS so that we can evalaute RSE
rse_safe <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]
  sse <- sum((obs - pred)^2)
  tss <- sum((obs - mean(obs))^2)
  if (tss == 0) return(NA_real_)
  sse / tss
}

rroc_curve <- function(y, yhat, n_points = 400) {
  stopifnot(length(y) == length(yhat))
  ok <- is.finite(y) & is.finite(yhat)
  y <- y[ok]; yhat <- yhat[ok]
  e <- yhat - y
  s_min <- -max(e); s_max <- -min(e)
  shifts <- seq(s_min - 1e-8, s_max + 1e-8, length.out = max(20, n_points))
  UNDER <- numeric(length(shifts)); OVER <- numeric(length(shifts))
  for (i in seq_along(shifts)) {
    ee <- e + shifts[i]
    OVER[i]  <- sum(ee[ee > 0])     # x-axis
    UNDER[i] <- sum(-ee[ee < 0])    # y-axis
  }
  data.frame(OVER = OVER, UNDER = UNDER, shift = shifts)
}

rroc_aoc <- function(rroc_df) {
  d <- rroc_df[order(rroc_df$OVER), ]
  sum((d$UNDER[-1] + d$UNDER[-nrow(d)]) * diff(d$OVER) / 2)
}

## ---- CV folds (yours) ----
set.seed(42)
K <- 50
folds_train <- createFolds(seq_len(nrow(pheno2)), k = K, returnTrain = TRUE)
index    <- folds_train
indexOut <- lapply(folds_train, function(tr) setdiff(seq_len(nrow(pheno2)), tr))

## ---- loop: outcome × variant with RMSE, RSE, RROC/AOC ----
rows <- list()
coef_txt <- file.path(output, "coefficients_printouts.txt")
if (file.exists(coef_txt)) file.remove(coef_txt)

# RROC PDF (replaces ROC PDF)
rroc_pdf <- file.path(output, "rroc_curves.pdf")
pdf(rroc_pdf, width = 7, height = 6)
on.exit(try(dev.off(), silent = TRUE), add = TRUE)

for (outcome in y_cols) {
  for (tag in names(variants)) {
    
    # Build predictor set for this variant
    if (tag == "Only_biomarkers") {
      keep <- intersect(variants[[tag]], predictor_pool)
      if (length(keep) == 0)
        warning(sprintf("[%s|%s] None of the requested 'Only_biomarkers' found in data.", outcome, tag))
      x_cols <- keep
    } else {
      drops <- variants[[tag]]
      missing_drops <- setdiff(drops, names(pheno2))
      if (length(missing_drops))
        warning(sprintf("[%s|%s] Drop list contains missing columns: %s",
                        outcome, tag, paste(missing_drops, collapse=", ")))
      x_cols <- setdiff(predictor_pool, drops)
    }
    
    # Restrict data to outcome + selected predictors
    df <- pheno2[, c(outcome, x_cols), drop = FALSE]
    
    # Train
    set.seed(42)
    cv_model <- train(
      as.formula(paste(outcome, "~ .")),
      data = df,
      method = "glmnet",
      trControl = trainControl(method = "cv", number = K,
                               index = index, indexOut = indexOut,
                               savePredictions = "final", verboseIter = FALSE),
      tuneGrid = tuning_grid,
      metric = "RMSE"
    )
    
    # Best tune + OOF predictions
    bt  <- cv_model$bestTune
    oof <- subset(cv_model$pred, alpha == bt$alpha & lambda == bt$lambda)
    
    # Metrics: RMSE, RSE (safe), R2
    rmse_val <- RMSE(oof$pred, oof$obs)
    rse_val  <- rse_safe(oof$obs, oof$pred)
    r2_val   <- 1 - rse_val
    
    # ---- RROC (Regression ROC) ----
    plot(
      d_rroc$OVER, d_rroc$UNDER,
      type = "l",
      xlab = "OVER-estimation",
      ylab = "UNDER-estimation",
      main = sprintf(
        "RROC - %s | %s\nRMSE=%.3f  RSE=%.3f  AOC=%.2f",
        outcome, tag, rmse_val, rse_val, aoc_val
      )
    )
    abline(h = 0, v = 0, lty = 3)
    
    
    
    # Collect row
    rows[[length(rows) + 1]] <- data.table(
      outcome = outcome,
      variant = tag,
      best_alpha = bt$alpha,
      best_lambda = bt$lambda,
      rmse_oof = rmse_val,
      r2_oof   = r2_val,
      rse_oof  = rse_val,
      aoc_rroc = aoc_val,
      n_features = length(x_cols)
    )
    
    # ---- Coefficients (raw printout + tidy CSV of non-zeros) ----
    co <- coef(cv_model$finalModel, s = bt$lambda)
    
    cat(sprintf("\n===== %s | %s =====\n", outcome, tag))
    print(co)
    
    cat(sprintf("\n===== %s | %s =====\n", outcome, tag), file = coef_txt, append = TRUE)
    capture.output(print(co), file = coef_txt, append = TRUE)
    
    cz <- as.matrix(co)
    coef_df <- data.table(feature = rownames(cz), coef = as.numeric(cz))
    coef_df_nz <- coef_df[coef != 0]
    if (nrow(coef_df_nz)) {
      coef_df_nz[, `:=`(outcome = outcome, variant = tag)]
      fwrite(coef_df_nz,
             file = file.path(output, "coefficients_nonzero.csv"),
             append = file.exists(file.path(output, "coefficients_nonzero.csv")))
    }
  }
}

## ---- save summaries ----
results <- rbindlist(rows)
fwrite(results, file.path(output, "glmnet_outcomes_by_variant.csv"))

# Macro CSV now includes mean RSE (macro_rse) and mean AOC (macro_aoc)
macro <- results[, .(
  macro_rmse = mean(rmse_oof, na.rm = TRUE),
  macro_r2   = mean(r2_oof,   na.rm = TRUE),
  macro_rse  = mean(rse_oof,  na.rm = TRUE),
  macro_aoc  = mean(aoc_rroc, na.rm = TRUE)
), by = variant][order(macro_rmse)]
fwrite(macro, file.path(output, "macro_by_variant.csv"))

cat("\nSaved:\n",
    "- metrics: ", file.path(output, "glmnet_outcomes_by_variant.csv"), "\n",
    "- macro:   ", file.path(output, "macro_by_variant.csv"), "\n",
    "- coefs:   ", file.path(output, "coefficients_printouts.txt"),
    " (raw coef() printouts)\n",
    "- coefs:   ", file.path(output, "coefficients_nonzero.csv"),
    " (tidy non-zero)\n",
    "- RROC PDF:", rroc_pdf, "\n", sep = "")


dev.list()        
while (!is.null(dev.list())) dev.off()
