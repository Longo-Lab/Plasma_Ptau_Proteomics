



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

# predictors
variants <- list(
  ALL             = character(0),              
  NO_PCSandAPOE   = c("APOE","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"),
  NO_PCS          = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8")
)

# Join all predictors from pheno2 and get rid of cell types 
predictor_pool <- setdiff(names(pheno2), y_cols)

# Tuning grid that tests different alphas and lambdas
tuning_grid <- expand.grid(
  alpha  = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),  
  lambda = seq(0.001, 0.1, length = 20)
)

# Drops values that are not finite and calculates SSE and TSS so that we can evaluate RSE
rse_safe <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]
  sse <- sum((obs - pred)^2)
  tss <- sum((obs - mean(obs))^2)
  if (tss == 0) return(NA_real_)
  sse / tss
}

# (Keep your scalar RROC/AOC helpers if you still want AOC numbers in the table)
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

## ---- libs for auditor RROC + ensure folder ----
suppressPackageStartupMessages({
  library(DALEX)
  library(auditor)
  library(ggplot2)   # for ggsave()
  library(caret)
  library(data.table)
})
dir.create(file.path(output, "rroc"), showWarnings = FALSE, recursive = TRUE)

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

# PDF device that will collect all auditor RROCs (one page per (outcome×variant))
rroc_pdf <- file.path(output, "rroc_curves.pdf")
pdf(rroc_pdf, width = 7, height = 6)
on.exit(try(dev.off(), silent = TRUE), add = TRUE)

for (outcome in y_cols) {
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
    
    # Metrics: RMSE, RSE (safe), R2
    rmse_val <- RMSE(oof$pred, oof$obs)
    rse_val  <- rse_safe(oof$obs, oof$pred)
    r2_val   <- 1 - rse_val
    
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
    
    # (Optional) Scalar AOC from our simple RROC to keep your table column
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
    "- RROC PNGs in: ", file.path(output, "rroc"), "\n",
    "- RROC PDF: ", rroc_pdf, "\n", sep = "")

# make sure no open devices remain
dev.list()
while (!is.null(dev.list())) dev.off()



####Testing the models##############

b <- readRDS(file.path(output, "glmnet_Neuron_NO_PCSandAPOE.rds"))

new_df <- a_test[, micro_model$x_cols, drop = FALSE] 
a_test <- a
cols <- c("SampleID", new_df)

ids <- a[order(DS1V)][, head(.SD, 1), by = "PIDN"][, SampleID]
f_ups <- a[!SampleID %in% ids]
columns_to_check <- c( "PlasmaPTau181","PlasmaAB142P","PlasmaAB140P","PlasmaABRatio","PlasmapTau217","PlasmapTau217_AB42Ratio",
                       "PlasmaGFAP","PlasmaNfL")
drops <- f_ups[!complete.cases(f_ups[, ..columns_to_check ]), SampleID]
a_test <- f_ups[!SampleID %in% drops]

test_data <- data.frame(a_test[, ..cols], row.names = "SampleID")
test_data <- test_data[complete.cases(test_data),] 


#All variables need to be numeric 
test_data$Gender <- as.numeric(test_data$Gender == "M")



# Setup for baseline visit only analysis ---------------------------
# Correction, keep only first visit where plasma available
a_test <- a_test[Visit == 2, .SD[1L], by = PIDN]

# checking for missed duplicates or bad visit labeling
a_test[duplicated(PIDN), ..qc_cols]


# handling incomplete cases
cols <- c("SampleID", new_df)
drops <- a_test[!complete.cases(a_test[, ..cols]), SampleID]
a_test <- a_test[!SampleID %in% drops]

#Setting up data for input 
pheno2_test <- data.frame(a_test[, ..new_df], row.names = "SampleID")
pheno2_test <- pheno2_test[complete.cases(pheno2_test),] 


#All variables need to be numeric 
pheno2_test$Gender <- as.numeric(pheno2_test$Gender == "M")
pheno2_test$APOE <- ifelse(grepl("4", pheno2_test$APOE), 1, 0) 

new_df <- a_test[, micro_model$x_cols, drop = FALSE] 
ids  <- new_df[["SampleID"]]

# columns to keep (everything except SampleID)
keep <- setdiff(names(new_df), "SampleID")

# data.table way: ensure a 2-D object, then coerce
X <- as.matrix(test_data)
# Predictions of outcomes
pred <- data.table(
  predict(b$model, newx = X, s = b$lambda),
  keep.rownames = "SampleID"
)
names(pred) <- c("SampleID", "outcomes")

# get original scores
OG<- as.data.table(x)

# merge & plot
dt <- merge(pred, OG, by = "SampleID", all.x = TRUE)
nrow(pred) == nrow(dt)

#Plot
test <- ggplot(data = dt, aes(x = Neuron, y = outcomes)) +
  geom_point() +
  geom_abline() +
  labs(x = "Observed", y = "Predicted",
                      title = paste("Predicted vs Observed —") +
  theme_minimal()
)
test
y_true <- dt$Neuron
y_pred <- dt$outcomes
#Assess accuracy 
rmse <- sqrt(mean((y_true - y_pred)^2))
# 0.2288543
# Mean Absolute Error
mae <- mean(abs(y_true - y_pred))
#0.1808942
# R-squared
r2 <- cor(y_true, y_pred)^2
#0.07439681

library(auditor)

micro <- merge(dt, data.table(train_df[, -1], keep.rownames = "SampleID"), by = "SampleID")

micro_exp <- DALEX::explain(b$model, 
  label = "neuron", 
  data = micro[, c(7:ncol(micro)), with = FALSE], 
  y = micro$outcomes
)
score_acc(micro_exp, y = dt$Neuron)
b_mp <- model_performance(b$model)


#testing individual model
outcome <- pheno2[, c(3)]
test_data <- pheno2[, c(5, 6, 16:23)]
neuron_test <- glmnet(
  x= test_data,
  y=outcome,
  alpha= 0,
  lambda=0.0531052631578947 
)

test2_data <- as.matrix(test_data)
pred <- predict(neuron_test, newx= test2_data)
RMSE(exp(pred), exp(outcome))
#testing on baseline: 0.01347016

## new prediction 
a_test[, Gender := as.numeric(Gender == "M")]
testing <- as.matrix(a_test[,c(1,10,11,20:27)], rownames = "SampleID")

pred <- data.table(
  predict(neuron_test, newx = testing),
  keep.rownames = "SampleID"
)
              

# ## aligned offical outcomes with 
merged_data <- merge(pred, x[, .(SampleID, Neuron)])
merged_data[, .(
  RMSE = RMSE(exp(s0), exp(Neuron)),
  mae = mean(abs(s0 - Neuron)),
  cor = cor.test(x = s0, y = Neuron)$estimate,
  p.value = cor.test(x = s0, y = Neuron)$p.value
)]
#0.01517091

mae <- mean(abs(merged_data$s0 - merged_data$Neuron))
#0.2089119
# R-squared
cor.test(x = merged_data$s0, y = merged_data$Neuron)
#0.2458497 

correlation <- cor(pred$s0, merged_data$Neuron)
r2 <- correlation^2

