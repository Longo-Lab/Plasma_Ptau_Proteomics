## ===== deps =====
req <- c("data.table","randomForestSRC","DALEX","auditor","ggplot2")
inst <- req[!req %in% installed.packages()[,"Package"]]
if (length(inst)) install.packages(inst, repos = "https://cloud.r-project.org")
invisible(lapply(req, library, character.only = TRUE))

## ===== config =====
y_names <- c("Microglia","Astrocyte","Neuron","Oligodendrocyte")

variants <- list(
  ALL             = character(0),
  NO_PCSandAPOE   = c("APOE","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"),
  NO_PCS = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8")
)

stopifnot(exists("output"))
dir.create(output, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output, "rroc_mRF"), showWarnings = FALSE)
dir.create(file.path(output, "rroc_mRF_pdf"), showWarnings = FALSE)

## ===== helpers =====
# split covs if user passed a comma string
if (is.character(covs) && length(covs) == 1 && grepl(",", covs)) {
  covs <- trimws(strsplit(covs, ",")[[1]])
}

rse_safe <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred); obs <- obs[ok]; pred <- pred[ok]
  sse <- sum((obs - pred)^2); tss <- sum((obs - mean(obs))^2)
  if (tss == 0) return(NA_real_)
  sse / tss
}

build_predictors <- function(ph, y_names, tag, variants_list) {
  pool  <- setdiff(colnames(ph), y_names)
  drops <- intersect(variants_list[[tag]], colnames(ph))
  setdiff(pool, drops)   # always DROP the listed columns
}

get_oob_preds <- function(fit, y_names) {
  # Robustly extract per-outcome OOB preds across randomForestSRC versions
  lapply(y_names, function(nm) {
    if (!is.null(fit$regrOutput) && !is.null(fit$regrOutput[[nm]]$predicted.oob))
      return(as.numeric(fit$regrOutput[[nm]]$predicted.oob))
    if (!is.null(fit$regrOutput) && !is.null(fit$regrOutput[[nm]]$predicted))
      return(as.numeric(fit$regrOutput[[nm]]$predicted))
    if (!is.null(fit$predicted.oob) && is.matrix(fit$predicted.oob) && nm %in% colnames(fit$predicted.oob))
      return(as.numeric(fit$predicted.oob[, nm]))
    # fallback (not strictly OOB, but keeps the pipeline working)
    as.numeric(predict(fit)$regrOutput[[nm]]$predicted)
  }) -> out
  names(out) <- y_names
  out
}

make_explainer <- function(y, yhat, label) {
  DALEX::explain(
    model = list(),  # dummy; predict_function returns supplied yhat
    data  = data.frame(.id = seq_along(yhat)),
    y     = y,
    predict_function = function(m, newdata) yhat[newdata$.id],
    label = label,
    verbose = FALSE
  )
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
  d <- rroc_df[order(rroc_df$OVER), ]
  if (nrow(d) < 2) return(NA_real_)
  sum((d$UNDER[-1] + d$UNDER[-nrow(d)]) * diff(d$OVER) / 2)
}

## ===== build a MULTIVARIATE pheno2 with ALL outcomes at once =====
needed <- unique(c("SampleID", y_names, covs, feature))
missing <- setdiff(needed, names(a))
if (length(missing)) stop("Missing columns in 'a': ", paste(missing, collapse=", "))

pheno2_mrf <- data.frame(a[, ..needed], row.names = "SampleID", check.names = FALSE)
pheno2_mrf <- pheno2_mrf[complete.cases(pheno2_mrf[, y_names]), ]  # keep rows with all outcomes observed

# Optional encodings (match your previous pipeline)
if ("Gender" %in% names(pheno2_mrf) && !is.numeric(pheno2_mrf$Gender)) {
  pheno2_mrf$Gender <- as.numeric(pheno2_mrf$Gender == "M")
}
if ("APOE" %in% names(pheno2_mrf) && !is.numeric(pheno2_mrf$APOE)) {
  pheno2_mrf$APOE <- ifelse(grepl("4", pheno2_mrf$APOE), 1, 0)
}

# sanity: all outcomes present and numeric
missing_y <- setdiff(y_names, colnames(pheno2_mrf))
if (length(missing_y)) stop("Outcomes not found in pheno2_mrf: ", paste(missing_y, collapse=", "))
if (!all(vapply(pheno2_mrf[, y_names], is.numeric, TRUE))) stop("All outcomes must be numeric.")

## ===== fit one multivariate RF PER VARIANT; compute RMSE/RSE; plot auditor RROC =====
all_variant_metrics <- list()

for (tag in names(variants)) {
  message("== Variant: ", tag)
  
  x_cols <- build_predictors(pheno2_mrf, y_names, tag, variants)
  
  # drop rows with NA in predictors for this variant (keep outcomes intact)
  keep_rows <- complete.cases(pheno2_mrf[, x_cols, drop = FALSE])
  df <- pheno2_mrf[keep_rows, c(y_names, x_cols), drop = FALSE]
  if (nrow(df) < 10) stop("Too few rows after NA filtering for variant ", tag)
  
  # multivariate formula
  rhs <- paste(x_cols, collapse = " + ")
  fml <- as.formula(paste0("Multivar(", paste(y_names, collapse = ", "), ") ~ ", rhs))
  
  set.seed(42)
  fit <- rfsrc(
    formula    = fml,
    data       = df,
    importance = TRUE,
    forest     = TRUE,
    na.action  = "na.omit"
  )
  
  # OOB predictions per outcome
  oob_list <- get_oob_preds(fit, y_names)
  
  # metrics + auditor RROC
  metrics <- data.table(outcome = y_names, RMSE_OOB = NA_real_, RSE_OOB = NA_real_, AOC_RROC = NA_real_)
  plots <- vector("list", length(y_names))
  
  for (j in seq_along(y_names)) {
    y    <- as.numeric(df[[y_names[j]]])
    yhat <- as.numeric(oob_list[[j]])
    
    n <- min(length(y), length(yhat))
    y <- y[seq_len(n)]; yhat <- yhat[seq_len(n)]
    
    rmse <- sqrt(mean((y - yhat)^2))
    rse  <- rse_safe(y, yhat)
    
    dr   <- rroc_curve(y, yhat)
    aoc  <- rroc_aoc(dr)
    
    metrics$RMSE_OOB[j] <- rmse
    metrics$RSE_OOB[j]  <- rse
    metrics$AOC_RROC[j] <- aoc
    
    expl <- make_explainer(y, yhat, sprintf("%s | %s (mRF OOB)", y_names[j], tag))
    resid_obj <- auditor::model_residual(expl)
    p <- auditor::plot_rroc(resid_obj) +
      ggplot2::ggtitle(sprintf("RROC — %s | %s", y_names[j], tag))
    plots[[j]] <- p
    
    # PNG per outcome×variant
    png_path <- file.path(output, "rroc_mRF", sprintf("rroc_mRF_%s_%s.png", y_names[j], tag))
    ggsave(png_path, p, width = 7, height = 6, dpi = 300)
  }
  
  # CSV metrics for this variant
  fwrite(metrics, file.path(output, sprintf("mrf_oob_metrics_%s.csv", tag)))
  
  # PDF with one page per outcome
  pdf_path <- file.path(output, "rroc_mRF_pdf", sprintf("rroc_mRF_%s.pdf", tag))
  pdf(pdf_path, width = 7, height = 6)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)
  for (p in plots) print(p)
  dev.off()
  
  all_variant_metrics[[tag]] <- metrics[, variant := tag]
}

# Macro summary across outcomes per variant
macro <- rbindlist(all_variant_metrics)[, .(
  macro_RMSE_OOB = mean(RMSE_OOB, na.rm = TRUE),
  macro_RSE_OOB  = mean(RSE_OOB,  na.rm = TRUE),
  macro_AOC_RROC = mean(AOC_RROC, na.rm = TRUE) 
), by = variant][order(macro_RMSE_OOB)]
fwrite(macro, file.path(output, "mrf_macro_metrics_by_variant.csv"))


cat("\nSaved:\n",
    "- Per-variant metrics: mrf_oob_metrics_<variant>.csv\n",
    "- Per-variant RROC PDFs: rroc_mRF_pdf/rroc_mRF_<variant>.pdf\n",
    "- Per-outcome RROC PNGs: rroc_mRF/rroc_mRF_<outcome>_<variant>.png\n",
    "- Macro summary: mrf_macro_metrics_by_variant.csv\n", sep = "")


#### Run model on ideal conditions 

fit <- rfsrc(
  formula    = fml,
  data       = df,
  importance = TRUE,
  forest     = TRUE,
  na.action  = "na.omit"
)


##Predict outcomes on validation set 


