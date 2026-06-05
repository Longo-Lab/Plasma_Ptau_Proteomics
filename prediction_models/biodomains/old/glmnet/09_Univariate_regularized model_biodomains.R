#!/usr/bin/env R


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


# Inputs ----------------
inpath <- "C:/Users/audweber/Box/plasma_ptau_proteomics" # don't repeat things you don't have to
input <- file.path(inpath, "dataProt_SS-205063.merged.QC.csv")
prots <- file.path(inpath, "prot_names.txt") 
pcs <- file.path(inpath, "genomic/mypcair.PCs.csv")
prediction_scores <- file.path(inpath, # wrap at 80 characters
  "celltype_scores/dataProt_SS-205063.treatad_ranksum_scaled.normalized_scores.csv"
)
output <- 'prediction_models/biodomains/glmnet'
covs <- 'Age,Gender,APOE,HybControlNormScale'   
outcome <- c("APP Metabolism",	"Apoptosis",	"Autophagy", "Cell Cycle", "Endolysosome",	"Immune Response",	"Lipid Metabolism",	"Metal Binding and Homeostasis",
             "Mitochondrial Metabolism",	"Myelination",	"Oxidative Stress",	"Proteostasis",	"Structural Stabilization",	"Synapse",	"Tau Homeostasis"	,"Vasculature")
 
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

# Functions ------------------
#' Recodes dummy categorical variables and can subset by vector of categories
#'
#' @param dat Data.table with the original data
#' @param dcol String representing the column to be split
#' @param keep Vector of strings that are the categories within dcol to extract. default is all c()
#' @returns modified data.table with new dummy variable columns named "dcol.category" after original column
#' 
recode_dummies = function(dat, dcol, keep=c()){
  dt = dat
  # check input
  if (length(keep) == 0){
    inds <- unique(dt[[dcol]])
  } else if (length(keep) > 0 & all(keep %in% unique(dt[[dcol]]))){
    inds <- unique(dt[[dcol]])[unique(dt[[dcol]]) %in% keep]
  } else {
    stop("Keep must be a vector of names in the original column")
  }
  # split out dummies
  dt[, c(paste(dcol, inds, sep=".")) := lapply(inds, function(x) ifelse(get(dcol) == x, 1, 0))]
  # reorder behind original variable
  v = grep(paste0("^", dcol, "$"), names(dt))[1]
  inds = paste(dcol, inds, sep=".")
  setcolorder(dt, c(names(dt)[1:v], inds, names(dt)[(v + 1):(ncol(dt) - length(inds))]))
  return(dt)
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

# Drops values that are not finite and calculates SSE and TSS 
# so that we can evaluate RSE
rse_safe <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]
  sse <- sum((obs - pred)^2)
  tss <- sum((obs - mean(obs))^2)
  if (tss == 0) return(NA_real_)
  sse / tss
}

# Main ------------------------------------------------------------------------------------------------------------
#Load data as datatable 
a <- fread(input)

# Map in PCs assumes your primary data table exists 
if (pcs != "") {
  print(sprintf("Reading PCs for inclusion from %s", pcs))
  x <- fread(pcs, colClasses = c(PIDN = "character"))
  a <- merge(a, x, by = "PIDN")
  covs <- unique(c(covs, names(a)[names(a) %like% "^PC\\d{1,2}$"])) # flexible for PCs
}

# Add cell scores
if (prediction_scores != "") {
  print(sprintf("Reading Prediction scores for inclusion from %s", prediction_scores))
  x <- fread(prediction_scores, colClasses = c(SampleID = "character"))
  a <- merge(a, x, by = "SampleID")
}

# filter outliars 
drop = a[abs(ConnectivityZscore) > 10, SampleID]
a = a[!SampleID %in% drop]

#Only keep baseline visits where plasma is avaliable 
a = a[order(DS1V)][, head(.SD, 1), by = "PIDN"]

# checking for missed duplicates or bad visit labeling
a[duplicated(PIDN), ..qc_cols]

# All variables need to be numeric & dummy coded
a[, Gender := as.numeric(Gender == "M")]
recode_dummies(a, "APOE", keep = c("e34", "e44"))
covs <- unique(c(covs[covs != "APOE"], "APOE.e34", "APOE.e44"))


# Handling incomplete cases
cols <- c("SampleID", outcome, covs, feature)

# Setting up data for input 
pheno <- a[complete.cases(a[, ..cols]), ..cols] # at this point everything is in pheno and is ordered, just keep calling this
pheno <- `rownames<-`(setDF(copy(a[complete.cases(a[, ..cols]), 
                      setdiff(cols, "SampleID"), with = FALSE])), 
                      make.unique(as.character(a[complete.cases(a[, ..cols]), SampleID])))

#Save a copy of pheno for later 
pheno2 <- pheno
 

names(pheno2) <- names(name_map <- setNames(names(pheno2), make.names(names(pheno2), unique = TRUE)))

# predictors
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
#rename outcomes since they have been changed for readability 
outcome <- c("APP.Metabolism","Apoptosis","Autophagy", "Cell.Cycle","Endolysosome",                 
             "Immune.Response","Lipid.Metabolism","Metal.Binding.and.Homeostasis", "Mitochondrial.Metabolism","Myelination",                  
             "Oxidative.Stress","Proteostasis","Structural.Stabilization","Synapse","Tau.Homeostasis","Vasculature") 

## ---- loop: outcome × variant with RMSE, RSE, RROC/AOC ------------------------------
rows <- list()
coef_txt <- file.path(output, "coefficients_printouts.txt")
if (file.exists(coef_txt)) file.remove(coef_txt)

# PDF device that will collect all auditor RROCs (one page per (outcome×variant))
rroc_pdf <- file.path(output, "rroc_curves.pdf")
pdf(rroc_pdf, width = 7, height = 6)
on.exit(try(dev.off(), silent = TRUE), add = TRUE)

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
    metric = "RMSE"
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
    
saveRDS(bundle, file = file.path(output, sprintf("glmnet_%s_%s.rds", outcome, tag)))
    
oof <- subset(cv_model$pred, alpha == bt$alpha & lambda == bt$lambda)
    
# Guard against empty OOF (can happen if factor levels / NA mismatch)
if (nrow(oof) == 0) {
    warning(sprintf("[%s|%s] No OOF rows after filtering bestTune; skipping.", outcome, tag))
    next
}
    
# Metrics: RMSE, RSE (safe), R2, MAE, MAPE

rmse_val <- RMSE(oof$pred, oof$obs)
rse_val  <- rse_safe(oof$obs, oof$pred)
r2_val   <- 1 - rse_val
mae_val <- mean(abs(oof$pred - oof$obs))
mape_val <- mean(abs((oof$pred - oof$obs) / oof$obs)) * 100
    
# ---- auditor RROC from OOF predictions (saved PNG + printed into the PDF) ----
expl <- DALEX::explain(
    model = list(),  # dummy
    data  = data.frame(.oof_id = seq_len(nrow(oof))),
    y     = oof$obs,
    predict_function = function(m, newdata) oof$pred[newdata$.oof_id],
    label = sprintf("%s | %s (glmnet OOF)", outcome, tag),
    verbose = FALSE
)

resid_obj <- auditor::model_residual(expl)
p_rroc <- auditor::plot_rroc(resid_obj) +
    ggplot2::ggtitle(sprintf("RROC — %s | %s\nRMSE=%.3f  RSE=%.3f",
                        outcome, tag, rmse_val, rse_val))
    
# Save a PNG per (outcome × variant)

png_path <- file.path(output, "rroc", sprintf("rroc_%s_%s.png", outcome, tag))

ggsave(png_path, p_rroc, width = 7, height = 6, dpi = 300)
    
# Also print into the open PDF device
print(p_rroc)

#Scalar AOC from our simple RROC to keep your table column
d_rroc  <- rroc_curve(oof$obs, oof$pred)
aoc_val <- rroc_aoc(d_rroc)
    

# Collect row
rows[[length(rows) + 1]] <- data.table(
    outcome     = outcome,
    variant     = tag,
    best_alpha  = bt$alpha,
    best_lambda = bt$lambda,
    rmse_oof    = rmse_val,
    r2_oof      = r2_val,
    rse_oof     = rse_val,
    aoc_rroc    = aoc_val,
    mae_val     = mae_val,
    mape_val    = mape_val,
    n_features  = length(x_cols)
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
      fwrite(
        coef_df_nz,
        file = file.path(output, "coefficients_nonzero.csv"),
        append = file.exists(file.path(output, "coefficients_nonzero.csv"))
      )
    }
  }
}

## ---- save summaries ----
results <- rbindlist(rows)
fwrite(results, file.path(output, "glmnet_outcomes_by_variant_INT.csv"))

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
    "- RROC PNGs in: ", file.path(output, "rroc"), "\n",
    "- RROC PDF: ", rroc_pdf, "\n", sep = "")

# make sure no open devices remain
dev.list()
while (!is.null(dev.list())) dev.off()


#Setting up outcome and test data 
outcome   <- as.matrix(pheno2[, 11, drop = FALSE])
test_data <- as.matrix(pheno2[, c(17:35)])
#Make sure the samples are in the same order!!
identical(rownames(outcome), rownames(test_data))

oxidative_stress_test <- glmnet(
  x= test_data,
  y=outcome,
  alpha= 0.1,
  lambda=0.0218421052631579
)

saveRDS(oxidative_stress_test, file.path(output, "final_oxidative_stress_glmnet_model.rds"))


## prediction
new_data <- subset(test_data, select = -c(1:9, 13:20))
new_data <- as.matrix(new_data)

oxidative_stress_pred <- data.frame(
  predict(oxidative_stress_test, newx = test_data)
)

              
# ## aligned official outcomes with 
idx <- match(rownames(oxidative_stress_pred), rownames(pheno2))

merged_data <- cbind(
  oxidative_stress_pred,
  pheno2[idx, , drop = FALSE]   # will insert NAs where no match
)

oxidative_stress_metrics <- with(merged_data, {
  data.frame(
    RMSE = RMSE(s0, Oxidative.Stress),
    mae = mean(abs(s0 - Oxidative.Stress), na.rm = TRUE),
    cor = unname(cor.test(x = s0, y = Oxidative.Stress)$estimate),
    p.value = cor.test(x = s0, y = Oxidative.Stress)$p.value,
    r2 = 1 - sum((Oxidative.Stress - s0)^2, na.rm = TRUE) / sum((Oxidative.Stress - mean(Oxidative.Stress, na.rm = TRUE))^2, na.rm = TRUE),
    MAPE = mean(abs((s0 - Oxidative.Stress) / pmax(abs(Oxidative.Stress), .Machine$double.eps)), na.rm = TRUE) * 100
  )
})

#Plot predictions as scatterplot 
oxidative_stress_predictions <- ggplot(data = merged_data, aes(x = Oxidative.Stress, y = s0)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  labs(x = "Observed", y = "Predicted",
       title = paste("Predicted vs Observed —") +
         theme_minimal()
  )

ggsave(file.path(output, paste0("glmnet_oxidative_stress_predicted_vs_observed_",".pdf")),
       plot = oxidative_stress_predictions, width = 7, height = 7, dpi = 300)










