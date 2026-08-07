#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 09-WSS_construction.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Generate Weighted Sum of Scores for STRING network modules
#'  
#'
#' Author: Robert B Butler III
#'
#' Date Created: 2026-02-11
#'
#' Copyright (c) 2026
#' Email: rrbutler@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'   v0.0.1 - initial build
#'   
#'   Based on the initial examination of eigengenes, there is too little variance
#'   captured in eigengenes to be very stratifying, so instead, weight by the 
#'   association to core biomarkers and STRING connectivity
#'
#' ---------------------------

#' load up the packages we will need:  (uncomment as required)
library(data.table)
library(optparse)
# library(WGCNA)
library(RNOmni)
# library(ggplot2)
# library(ggridges)
# library(patchwork)
# library(gtools)
# library(igraph)
# library(stringr)


# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "module_scores/STRING_merged_scaled_scores.csv",
    help = "Biomarker cluster results to use for groups"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project [default: Nameset of input]"
  ),
  make_option(
    c("--INT"),
    action = "store_true",
    default = FALSE,
    help = "Inverse rank normalize expression [default: z-scale log2]"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "module_scores",
    help = "File path to output directory [default: %default]"
  ),
  make_option(
    c("-r", "--ref"),
    type = "character",
    default = "dataProt_SS-205063.merged.QC.csv",
    help = "Path to proteomic data [default: %default]"
   ),
  make_option(
    c("-p", "--prots"),
    type = "character",
    default = "dataProt_SS-205063.probe_list.csv",
    help = "Path to probe/GENE table for SOMA [default: %default]"
   ),
  make_option(
    c("-t", "--threads"),
    default = 0,
    help = "num threads [default: parallelly::availableCores()]"
  )
)

args_list <- parse_args(OptionParser(option_list = option_list))

# Nameset default case
if(is.null(args_list$nameset)){
  args_list$nameset <- nameset_default(args_list$input)
}

# read in fields
fields <- c(
  "input",
  "nameset",
  "INT",
  "output",
  "threads",
  "prots",
  "ref"
)

# check for values
for (field in fields) {
  if (!field %in% names(args_list)) {
    cat(str_c('Argument "', field, '" must be provided\n'), fill = TRUE)
    stop(parse_args(OptionParser(option_list = option_list), args = "--help"),
      call. = FALSE
    )
  }

  assign(field, args_list[[field]])
  print(paste0(field, ": ", get(field)))
}

# check the current active plan
if (threads == 0) {
  threads <- as.integer(parallelly::availableCores()) - 2
} else if (is.numeric(threads) & threads != 0) {
  threads <- as.integer(threads)
} else {
  cat("threads must be an integer value")
  stop(parse_args(OptionParser(option_list = option_list), args = "--help"),
    call. = FALSE
  )
}
print(sprintf("Setting up %d threads for samples", threads))

# check outpath status
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}


# Functions -------------------------
#' Calculate Weighted Sum Scores Using Signed Protein Values
#'
#' This function calculates module scores by sign-flipping proteins based on 
#' their biomarker associations (mean_beta) and weighting by alpha centrality.
#' This approach aligns proteins in the direction of core biomarker progression
#' and weights by network importance and biomarker association strength.
#'
#' @param expr_matrix A numeric matrix with samples as rows and genes/proteins 
#'   as columns (INT-normalized). Column names should match the 'symbol' values 
#'   in the result table.
#' @param result A data.table containing module information with required columns:
#'   module, symbol (gene identifiers), mean_beta (for sign-flipping),
#'   mean_alpha_scaled (for weighting).
#'
#' @return A list containing:
#'   \item{scores}{Matrix of weighted module scores (samples x modules)}
#'   \item{summary_stats}{data.table with summary statistics for each module score}
#'   \item{module_sd}{Named vector of standard deviations for each module}
#'   \item{module_proteins}{List of protein names in each module}
#'
#' @details 
#' For each module, proteins are sign-flipped based on sign(mean_beta) to align
#' all proteins in the direction of disease progression. The signed values are
#' then weighted by mean_alpha_scaled (network centrality + association strength)
#' and averaged to produce a module score.
#' 
#' Higher module scores indicate greater module activity in the direction of
#' disease progression.
#'
#' @examples
#' # Calculate weighted module scores
#' module_results <- calculate_weighted_module_scores(expr_matrix, cluster_result)
#' 
#' # Access score matrix
#' scores <- module_results$scores
#' 
#' # View summary statistics
#' print(module_results$summary_stats)
#' 
#' # Check standard deviations
#' print(module_results$module_sd)
#'
#' @import data.table
#' @export
calculate_weighted_module_scores <- function(expr_matrix, result) {
  
  # Get unique modules
  modules <- unique(result$module)
  
  # Initialize output structures
  score_matrix <- matrix(NA, 
                         nrow = nrow(expr_matrix), 
                         ncol = length(modules),
                         dimnames = list(rownames(expr_matrix), modules))
  
  summary_list <- list()
  sd_vector <- numeric(length(modules))
  names(sd_vector) <- modules
  protein_list <- list()
  
  # Calculate scores for each module
  for(mod in modules) {
    
    # Get proteins in this module
    mod_proteins <- result[module == mod]
    protein_ids <- mod_proteins$symbol
    
    # Check which proteins are in expression matrix
    available_proteins <- protein_ids[protein_ids %in% colnames(expr_matrix)]
    
    if(length(available_proteins) == 0) {
      warning(sprintf("Module %s has no proteins in expression matrix", mod))
      next
    }
    
    # Store protein list
    protein_list[[mod]] <- available_proteins
    
    # Subset expression matrix and module info
    mod_expr <- expr_matrix[, available_proteins, drop = FALSE]
    mod_info <- mod_proteins[symbol %in% available_proteins]
    
    # Ensure order matches
    mod_info <- mod_info[match(available_proteins, symbol)]
    
    # Get sign from mean_beta
    sign_vector <- sign(mod_info$mean_beta)
    
    # Get weights from mean_alpha_scaled
    weights <- mod_info$mean_alpha_scaled
    
    # Handle missing or zero weights
    if(any(is.na(weights)) || sum(weights) == 0) {
      warning(sprintf("Module %s has missing or zero weights, using equal weights", mod))
      weights <- rep(1, length(weights))
    }
    
    # Sign-flip proteins to align with disease progression
    signed_expr <- sweep(mod_expr, 2, sign_vector, "*")
    
    # Calculate weighted sum
    weighted_sum <- signed_expr %*% weights
    
    # Normalize by sum of weights
    module_score <- as.vector(weighted_sum / sum(weights))
    
    # Store in matrix
    score_matrix[, mod] <- module_score
    
    # Calculate summary statistics
    score_summary <- summary(module_score)
    summary_list[[mod]] <- as.list(score_summary)
    
    # Calculate SD
    sd_vector[mod] <- sd(module_score)
  }
  
  # Convert summary list to data.table
  summary_dt <- rbindlist(summary_list, idcol = "module")
  
  # Clean up column names if needed
  setnames(summary_dt, 
           old = c("Min.", "1st Qu.", "3rd Qu.", "Max."),
           new = c("Min", "Q1", "Q3", "Max"),
           skip_absent = TRUE)
  
  # Return results
  return(list(
    scores = score_matrix,
    summary_stats = summary_dt,
    module_sd = sd_vector,
    module_proteins = protein_list
  ))
}



# Main -----------------------------------------------------------------------
# Stage data for analysis --------------------------
result <- fread(input)

# for genes in multiple modules, keep by top mean_alpha_specificity
result <- result[order(-mean_alpha_scaled)][, .SD[1], by = symbol]
result[, module := paste(Biomarker, cluster, sep = ".")]

# filter singles/doubles
keep <- result[, .N, by = module][N > 2, module]
result <- result[module %in% keep]

# split string and explode piped terms into separate rows
prot <- fread(prots)

# drop outlier samples
a <- fread(ref)
drop = a[abs(ConnectivityZscore) > 10, SampleID]
a = a[!SampleID %in% drop]

# Post filter simple table
cols <- c(2, 28:ncol(a))
b <- merge(
  x = prot[, .(ix, GENE)], 
  y = transpose(a[, ..cols], make.names = "SampleID", keep.names = "ix"), 
  all.x = TRUE
)
exp <- as.matrix(
  x = transpose(b[, -1], make.names = "GENE", keep.names = "SampleID"), 
  rownames = "SampleID"
)

if (INT) {
  # Baseline is fit in-sample; follow-up is projected onto that frozen
  # baseline reference so no hold-out information leaks into the
  # baseline-derived WSS scores. Cohort membership is taken from the
  # existing baseline/follow-up clinical roster files (SampleID column
  # only -- not their module-score columns, which these scores feed
  # into), not a naive first-visit-per-PIDN rule: the raw merged.QC.csv
  # includes ~100 additional samples (e.g. non-genotyped/non-ADRC) that
  # the rest of the pipeline excludes from the baseline/follow-up cohorts
  # via an inner join against the genomic PCs / cohort roster files, so
  # first-visit-per-PIDN on the full pooled file overcounts both cohorts.
  # See issue #76.
  baseline_roster <- fread("dataProt_SS-205063.baseline_clin.csv", select = "SampleID")$SampleID
  followup_roster <- fread("dataProt_SS-205063.follow-up_clin.csv", select = "SampleID")$SampleID
  baseline_ids <- intersect(baseline_roster, a$SampleID)
  followup_ids <- intersect(followup_roster, a$SampleID)

  # Out-of-sample projection matching RNOmni::RankNorm (k = 0.375,
  # ties = "average") exactly for in-sample points; for a new point,
  # treats it as a single insertion into the frozen baseline reference.
  project_rank_norm <- function(train_vals, new_vals, k = 0.375) {
    n <- length(train_vals)
    sorted_train <- sort(train_vals)
    r_new <- findInterval(new_vals, sorted_train) + 1
    denom <- (n + 1) - 2 * k + 1
    qnorm((r_new - k) / denom)
  }

  exp_baseline <- exp[baseline_ids, , drop = FALSE]
  exp_followup <- exp[followup_ids, , drop = FALSE]
  for (j in seq_len(ncol(exp))) {
    train_vals <- exp_baseline[, j]
    exp_followup[, j] <- project_rank_norm(train_vals, exp_followup[, j])
    exp_baseline[, j] <- RankNorm(train_vals)
  }
  nameset <- paste(nameset, "INT", sep = ".")
} else {
  # z-scaling
  exp <- scale(log2(exp))
  baseline_ids <- rownames(exp)
  followup_ids <- character(0)
  exp_baseline <- exp
  exp_followup <- exp[0, , drop = FALSE]
}

# Run the analysis ----------------------------------
# Calculate weighted sum scores separately per cohort (baseline fit
# in-sample; follow-up projected onto the baseline reference when INT).
WSS_result_baseline <- calculate_weighted_module_scores(exp_baseline, result)
WSS_result_followup <- calculate_weighted_module_scores(exp_followup, result)

# Write baseline and follow-up scores to two separate files (not
# concatenated) so downstream consumers never see a pooled reference.
fwrite(
  x = as.data.table(WSS_result_baseline$scores, keep.rownames = "SampleID"),
  file = file.path(output, paste(nameset, "baseline", "WSS.csv", sep = "."))
)
saveRDS(WSS_result_baseline,
  file = file.path(output, paste(nameset, "baseline", "WSS.rds", sep = "."))
)
fwrite(
  x = as.data.table(WSS_result_followup$scores, keep.rownames = "SampleID"),
  file = file.path(output, paste(nameset, "follow-up", "WSS.csv", sep = "."))
)
saveRDS(WSS_result_followup,
  file = file.path(output, paste(nameset, "follow-up", "WSS.rds", sep = "."))
)

# # Access components
# module_scores <- results$scores  # Sample x Module matrix
# summary_table <- results$summary_stats  # Summary stats per module
# module_sds <- results$module_sd  # SD per module
# protein_lists <- results$module_proteins  # Proteins in each module

# # Quick diagnostics
# print(results$summary_stats)
# print(sort(results$module_sd, decreasing = TRUE))  # Which modules have most variance?
