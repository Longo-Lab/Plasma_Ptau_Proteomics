#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 08-Cell-type_scoring.R
#'
#' Version: 0.0.3
#'
#' Purpose of script: Score celltype for each sample using weighted GSEA method
#'  
#'
#' Author: Robert B Butler III
#'
#' Date Created: 2025-09-30
#'
#' Copyright (c) 2025
#' Email: rrbutler@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'   v0.0.3 - bugfixes for NA spec values (<- 0), messages for thresholds
#'   v0.0.2 - updates to drop for loops and speed processing
#'   
#'   
#'   
#'
#' ---------------------------

#' load up the packages we will need:  (uncomment as required)
library(Rcpp)
library(RcppParallel)
# library(bit)
library(data.table)
library(optparse)
library(ggplot2)
library(patchwork)


# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Biomarker results to use for probe filtering, keeps most significant"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project [default: Nameset of input]"
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
    default = "cell_type_prot/Sharma_etal.EWCE.human.lvl1.spec_matr.txt",
    help = "Path to specificity matrix [default: %default]"
   ),
  make_option(
    c("-t", "--threads"),
    default = 0,
    help = "num threads [default: parallelly::availableCores()]"
   ),
  make_option(
    c("-p", "--plot"), 
    action = "store_true", 
    default = FALSE,
    help = "Plot a random sampling of 10 individuals in output folder"
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
  "output",
  "ref",
  "threads",
  "nameset",
  "plot"
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

# input variables
input <- strsplit(input, ",")[[1]]

# check outpath status
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}


# Functions -------------------------
#' @title Process a single sample for specificity-weighted GSEA with optimizations
#'
#' @description 
#' Internal function that processes a single sample for specificity-weighted GSEA analysis.
#' Implements multiple optimizations including compiled code, matrix operations,
#' bit operations, and algorithmic improvements.
#'
#' @param sample_idx Integer. Index of the sample to process.
#' @param sample_names Character vector. Names of all samples.
#' @param expr_matrix Numeric matrix. Expression data.
#' @param specificity_matrix Numeric matrix. Cell type specificity scores.
#' @param thresholds Numeric vector. Specificity thresholds for each cell type.
#' @param above_threshold_matrix Logical matrix. Which proteins are above threshold.
#' @param cell_types Character vector. Names of cell types to analyze.
#' @param min_size Integer. Minimum number of proteins above threshold.
#' @param p_weight Numeric. Exponent for the ranking metric.
#' @param nperm Integer. Number of permutations for significance testing.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A list containing enrichment results for the specified sample.
#'
#' @examples
#' \dontrun{
#' # This is an internal function not meant to be called directly
#' # See specificityWeightedGSEA_optimized for usage examples
#' 
#' # Example data (used internally)
#' expr_data <- matrix(
#'   c(5.2, 3.8, 0.5, 2.1, 4.5, 0.2, 3.3, 1.2),
#'   ncol = 1,
#'   dimnames = list(
#'     c("ProteinA", "ProteinB", "ProteinC", "ProteinD", 
#'       "ProteinE", "ProteinF", "ProteinG", "ProteinH"),
#'     c("Sample1")
#'   )
#' )
#' 
#' spec_data <- matrix(
#'   c(0.9, 0.1, 0.0, 0.3, 0.1, 0.0, 0.8, 0.2,
#'     0.1, 0.8, 0.0, 0.2, 0.7, 0.1, 0.0, 0.4),
#'   ncol = 2,
#'   dimnames = list(
#'     c("ProteinA", "ProteinB", "ProteinC", "ProteinD",
#'       "ProteinE", "ProteinF", "ProteinG", "ProteinH"),
#'     c("Neuron", "Astrocyte")
#'   )
#' )
#' 
#' # Pre-compute values that would be calculated by the main function
#' cell_types <- colnames(spec_data)
#' thresholds <- c(Neuron = 0.1, Astrocyte = 0.1)
#' above_threshold_matrix <- matrix(spec_data > 0.1, nrow = nrow(spec_data),
#'                              dimnames = list(rownames(spec_data), cell_types))
#' 
#' # Process a single sample (for internal use only)
#' results <- .process_sample_optimized(
#'   sample_idx = 1,
#'   sample_names = "Sample1",
#'   expr_matrix = expr_data,
#'   specificity_matrix = spec_data,
#'   thresholds = thresholds,
#'   above_threshold_matrix = above_threshold_matrix,
#'   cell_types = cell_types,
#'   min_size = 2,  # Small for example
#'   p_weight = 1,
#'   nperm = 10,    # Small for example
#'   verbose = TRUE
#' )
#' }
#'
#' @keywords internal
.process_sample_optimized <- function(sample_idx, sample_names, expr_matrix, specificity_matrix,
                                   thresholds, above_threshold_matrix, cell_types,
                                   min_size, p_weight, nperm, verbose) {
  sample_name <- sample_names[sample_idx]
  if (verbose) message(sprintf("Processing sample %s (%d of %d)", 
                            sample_name, sample_idx, length(sample_names)))
  
  # Get expression values for this sample
  expr_values <- expr_matrix[, sample_idx]
  names(expr_values) <- rownames(expr_matrix)
  
  # Rank proteins based on expression
  protein_ranks <- sort(expr_values, decreasing = TRUE, index.return = TRUE)
  ranked_proteins <- rownames(expr_matrix)[protein_ranks$ix]
  n_proteins <- length(ranked_proteins)
  
  # Create matrix of ranked specificity scores
  spec_scores_ranked <- matrix(0, nrow = n_proteins, ncol = length(cell_types),
                            dimnames = list(ranked_proteins, cell_types))
  
  for (i in seq_along(cell_types)) {
    ct <- cell_types[i]
    spec_scores_ranked[, i] <- specificity_matrix[ranked_proteins, ct]
  }
  
  # REMOVED: Bit operations - just use standard logical matrices instead
  
  # Create matrix of above threshold indicators in ranked order
  above_threshold_ranked <- matrix(FALSE, nrow = n_proteins, ncol = length(cell_types),
                                dimnames = list(ranked_proteins, cell_types))
  
  for (i in seq_along(cell_types)) {
    above_threshold_ranked[, i] <- above_threshold_matrix[ranked_proteins, i]
  }
  
  # Compute normalization factors for each cell type
  N_R_values <- numeric(length(cell_types))
  names(N_R_values) <- cell_types
  
  N_non_specific_values <- numeric(length(cell_types))
  names(N_non_specific_values) <- cell_types
  
  valid_cell_types <- logical(length(cell_types))
  names(valid_cell_types) <- cell_types
  
  # Pre-compute powers if p_weight is not 1
  if (p_weight != 1) {
    spec_powers <- spec_scores_ranked^p_weight
  } else {
    spec_powers <- spec_scores_ranked
  }
  
  for (i in seq_along(cell_types)) {
    proteins_above <- above_threshold_ranked[, i]
    count_above <- sum(proteins_above)
    
    if (count_above >= min_size) {
      if (p_weight != 1) {
        N_R_values[i] <- sum(spec_powers[proteins_above, i])
      } else {
        N_R_values[i] <- sum(spec_scores_ranked[proteins_above, i])
      }
      N_non_specific_values[i] <- n_proteins - count_above
      valid_cell_types[i] <- TRUE
    } else {
      N_R_values[i] <- NA
      N_non_specific_values[i] <- NA
      valid_cell_types[i] <- FALSE
    }
  }
  
  # Initialize results containers
  enrichment_scores <- numeric(length(cell_types))
  names(enrichment_scores) <- cell_types
  
  running_scores <- matrix(0, nrow = n_proteins, ncol = length(cell_types))
  rownames(running_scores) <- ranked_proteins
  colnames(running_scores) <- cell_types
  
  leading_edge <- vector("list", length(cell_types))
  names(leading_edge) <- cell_types
  
  max_positions <- numeric(length(cell_types))
  names(max_positions) <- cell_types
  
  # Calculate enrichment scores for valid cell types using Rcpp
  for (i in which(valid_cell_types)) {
    ct <- cell_types[i]
    
    # KEEP: Use Rcpp for running sum calculation
    rs <- calcRunningSum(
      spec_scores_ranked[, i], 
      above_threshold_ranked[, i], 
      N_R_values[i], 
      N_non_specific_values[i], 
      p_weight
    )
    
    # Store running score
    running_scores[, i] <- rs
    
    # Find maximum deviation
    max_idx <- which.max(abs(rs))
    max_val <- rs[max_idx]
    
    # Store results
    enrichment_scores[i] <- max_val
    max_positions[i] <- max_idx
    
    # Store leading edge proteins
    if (max_val > 0) {
      leading_edge[[i]] <- ranked_proteins[1:max_idx]
    } else {
      leading_edge[[i]] <- ranked_proteins[max_idx:n_proteins]
    }
  }
  
  # Permutation testing for each cell type
  p_values <- rep(NA, length(cell_types))
  names(p_values) <- cell_types
  normalized_scores <- enrichment_scores
  
  if (nperm > 0) {
    # Matrix operations for permutation testing
    all_perms <- matrix(0, nrow = nperm, ncol = n_proteins)
    for (p in 1:nperm) all_perms[p,] <- sample(n_proteins)
    
    # Pre-allocate matrix for null distribution
    valid_indices <- which(valid_cell_types)
    null_dist_matrix <- matrix(0, nrow = nperm, ncol = length(valid_indices))
    colnames(null_dist_matrix) <- cell_types[valid_indices]
    
    # Process each valid cell type
    for (idx in seq_along(valid_indices)) {
      i <- valid_indices[idx]
      ct <- cell_types[i]
      
      # Get pre-computed values for this cell type
      ct_spec_scores <- spec_scores_ranked[, i]
      ct_above_threshold <- above_threshold_ranked[, i]
      ct_N_R <- N_R_values[i]
      ct_N_non_specific <- N_non_specific_values[i]
      
      # Process all permutations for this cell type
      for (p in 1:nperm) {
        perm_idx <- all_perms[p,]
        
        # KEEP: Use Rcpp for max deviation calculation
        null_dist_matrix[p, idx] <- calcMaxDeviation(
          ct_spec_scores[perm_idx],
          ct_above_threshold[perm_idx],
          ct_N_R,
          ct_N_non_specific,
          p_weight
        )
      }
    }
    
    # Calculate p-values and normalize scores
    for (idx in seq_along(valid_indices)) {
      i <- valid_indices[idx]
      es <- enrichment_scores[i]
      null_dist <- null_dist_matrix[, idx]
      
      # Calculate p-value with pseudocount
      if (es >= 0) {
        p_values[i] <- (sum(null_dist >= es) + 1) / (nperm + 1)
      } else {
        p_values[i] <- (sum(null_dist <= es) + 1) / (nperm + 1)
      }
      
      # Normalize enrichment score
      null_mean <- mean(null_dist)
      null_sd <- sd(null_dist)
      if (null_sd > 0) {
        normalized_scores[i] <- (es - null_mean) / null_sd
      }
    }
  }
  
  # Adjust p-values for multiple testing
  fdr <- p.adjust(p_values, method = "BH")
  
  # Return results for this sample
  return(list(
    sample = sample_name,
    enrichment_scores = enrichment_scores,
    normalized_scores = normalized_scores,
    p_values = p_values,
    fdr = fdr,
    leading_edge = leading_edge,
    running_scores = running_scores,
    ranked_proteins = ranked_proteins,
    thresholds_used = thresholds[names(enrichment_scores)]
  ))
}


#' @title Calculate Cell Type Enrichment Using Optimized Specificity-Weighted GSEA
#' 
#' @description
#' This function analyzes all samples in an expression matrix using highly optimized
#' specificity-weighted GSEA to calculate cell type enrichment scores. Includes
#' optimizations with compiled code, matrix operations, bit operations, and 
#' algorithmic improvements.
#'
#' @param expr_matrix Numeric matrix of protein expression values. Rows are proteins,
#'        columns are samples.
#' @param specificity_matrix Numeric matrix of cell type specificity scores (0-1 scale).
#'        Rows are proteins, columns are cell types.
#' @param rank_metric Character. Method for ranking proteins: "expression" (default) or "fold_change"
#' @param nperm Integer. Number of permutations for significance testing. Default is 1000.
#' @param p_weight Numeric. Exponent for the ranking metric. Default is 1.
#' @param specificity_threshold Numeric or character. Threshold for considering proteins as "specific".
#'        Can be a numeric value (0-1) or "median". Default is "median".
#' @param min_size Integer. Minimum number of proteins above threshold. Default is 5.
#' @param verbose Logical. Whether to print progress messages. Default is TRUE.
#' @param parallel Logical. Whether to use parallel processing for samples. Default is TRUE.
#' @param n_cores Integer. Number of CPU cores to use. Default is NULL (auto-detect).
#'
#' @return A named list where each element contains the enrichment analysis results for one sample.
#'         Each sample's results include:
#'         \itemize{
#'           \item{\code{sample}}: {Name of the sample}
#'           \item{\code{enrichment_scores}}: {Vector of enrichment scores for each cell type}
#'           \item{\code{normalized_scores}}: {Vector of normalized enrichment scores}
#'           \item{\code{p_values}}: {Vector of p-values from permutation testing}
#'           \item{\code{fdr}}: {Vector of FDR-adjusted p-values}
#'           \item{\code{leading_edge}}: {List of leading edge proteins for each cell type}
#'           \item{\code{running_scores}}: {Matrix of running scores for visualization}
#'           \item{\code{ranked_proteins}}: {Vector of proteins in ranked order}
#'           \item{\code{thresholds_used}}: {Thresholds used for each cell type}
#'         }
#'
#' @examples
#' # Create example expression dataset (proteins x samples)
#' expr_data <- matrix(
#'   c(5.2, 3.8, 0.5, 2.1, 4.5, 0.2, 3.3, 1.2,  # Sample 1
#'     4.8, 4.2, 0.4, 1.9, 5.1, 0.3, 3.5, 1.0,  # Sample 2
#'     6.1, 3.5, 0.6, 2.5, 4.0, 0.1, 3.8, 0.8), # Sample 3
#'   ncol = 3,
#'   byrow = FALSE,
#'   dimnames = list(
#'     c("ProteinA", "ProteinB", "ProteinC", "ProteinD", 
#'       "ProteinE", "ProteinF", "ProteinG", "ProteinH"),
#'     c("Sample1", "Sample2", "Sample3")
#'   )
#' )
#' 
#' # Create example cell type specificity dataset (proteins x cell types)
#' spec_data <- matrix(
#'   c(0.9, 0.1, 0.0, 0.3, 0.1, 0.0, 0.8, 0.2,  # Neuron
#'     0.1, 0.8, 0.0, 0.2, 0.7, 0.1, 0.0, 0.4,  # Astrocyte
#'     0.0, 0.2, 0.9, 0.5, 0.1, 0.8, 0.0, 0.3), # Microglia
#'   ncol = 3,
#'   byrow = FALSE,
#'   dimnames = list(
#'     c("ProteinA", "ProteinB", "ProteinC", "ProteinD",
#'       "ProteinE", "ProteinF", "ProteinG", "ProteinH"),
#'     c("Neuron", "Astrocyte", "Microglia")
#'   )
#' )
#' 
#' # Run analysis with minimal permutations for the example
#' \dontrun{
#' results <- specificityWeightedGSEA_optimized(
#'   expr_matrix = expr_data,
#'   specificity_matrix = spec_data,
#'   nperm = 100,  # Reduced for quick example
#'   specificity_threshold = "median",
#'   verbose = TRUE
#' )
#' 
#' # Access results for a specific sample
#' sample1_results <- results[["Sample1"]]
#' 
#' # Extract enrichment scores for all samples
#' all_scores <- sapply(results, function(x) x$enrichment_scores)
#' 
#' # Find significant cell types (FDR < 0.05) for each sample
#' significant_cell_types <- lapply(results, function(x) {
#'   cell_types <- names(x$fdr)[x$fdr < 0.05]
#'   if(length(cell_types) > 0) cell_types else "None"
#' })
#' 
#' # Extract the leading edge proteins for a specific cell type in Sample1
#' if("Neuron" %in% names(results[["Sample1"]]$leading_edge)) {
#'   neuron_proteins <- results[["Sample1"]]$leading_edge[["Neuron"]]
#'   print(paste("Leading edge proteins for Neuron in Sample1:", 
#'              paste(neuron_proteins, collapse=", ")))
#' }
#' }
#'
#' @importFrom parallel makeCluster stopCluster parLapply detectCores clusterExport
#' @importFrom stats median p.adjust
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @export
specificityWeightedGSEA <- function(expr_matrix, specificity_matrix, 
                                   rank_metric = "expression", nperm = 1000, 
                                   p_weight = 1, specificity_threshold = "median",
                                   min_size = 5, verbose = TRUE, 
                                   parallel = TRUE, n_cores = NULL) {
  
  start_time <- Sys.time()
  
  # Check for required packages
  if (!requireNamespace("Rcpp", quietly = TRUE) || 
      !requireNamespace("RcppParallel", quietly = TRUE)) {
    stop("Packages 'Rcpp' and 'RcppParallel' are required")
  }
  
  # Input validation
  if (!is.matrix(expr_matrix) || !is.matrix(specificity_matrix)) {
    stop("Expression and specificity inputs must be matrices")
  }
  
  # Check for matching protein names
  if (is.null(rownames(expr_matrix)) || is.null(rownames(specificity_matrix))) {
    stop("Both matrices must have rownames identifying proteins")
  }
  
  # Find common proteins and subset matrices
  common_proteins <- intersect(rownames(expr_matrix), rownames(specificity_matrix))
  if (length(common_proteins) == 0) {
    stop("No common proteins found between expression and specificity matrices")
  }
  expr_matrix <- expr_matrix[common_proteins, , drop = FALSE]
  specificity_matrix <- specificity_matrix[common_proteins, , drop = FALSE]
  
  # Get all sample names
  sample_names <- colnames(expr_matrix)
  n_samples <- length(sample_names)
  
  if (verbose) {
    message(sprintf("Processing %d samples across %d proteins and %d cell types", 
                  n_samples, nrow(expr_matrix), ncol(specificity_matrix)))
  }
  
  # Set up parallel processing if requested
  if (parallel) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning("Package 'parallel' not available. Using sequential processing.")
      parallel <- FALSE
    } else {
      if (is.null(n_cores)) {
        n_cores <- max(1, parallel::detectCores() - 2)
      }
      n_cores <- min(n_cores, n_samples) # Don't use more cores than samples
      if (verbose) message(sprintf("Using parallel processing with %d cores", n_cores))
    }
  }
  
  # Pre-compute thresholds and other invariant values once for all samples
  if (verbose) message("Pre-computing cell type thresholds...")
  
  cell_types <- colnames(specificity_matrix)
  n_cell_types <- length(cell_types)
  
  # Calculate thresholds - this only needs to be done once
  thresholds <- numeric(n_cell_types)
  names(thresholds) <- cell_types
  
  for (i in seq_along(cell_types)) {
    ct <- cell_types[i]
    if (is.character(specificity_threshold) && specificity_threshold == "median") {
      non_zero_specs <- specificity_matrix[specificity_matrix[, ct] > 0, ct]
      thresholds[i] <- if (length(non_zero_specs) > 0) median(non_zero_specs) else 0
    } else {
      thresholds[i] <- specificity_threshold
    }
  }
  if (verbose) message(paste(capture.output(print(thresholds)), collapse = "\n"))
  
  # Create above_threshold matrix - this also only needs to be done once
  above_threshold_matrix <- matrix(FALSE, nrow = nrow(specificity_matrix), 
                                ncol = n_cell_types,
                                dimnames = list(rownames(specificity_matrix), cell_types))
  
  for (i in seq_along(cell_types)) {
    above_threshold_matrix[, i] <- specificity_matrix[, cell_types[i]] > thresholds[i]
  }
  
  # Process all samples, either in parallel or sequentially
  if (parallel) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    
    # Load RcppParallel on all worker nodes
    parallel::clusterEvalQ(cl, {
      library(Rcpp)
      library(RcppParallel)
      RcppParallel::setThreadOptions()
    })
    
  # Either export the functions, or source the C++ file on each worker
  if (file.exists("gsea_functions.cpp")) {
    parallel::clusterCall(cl, function() {
      Rcpp::sourceCpp("gsea_functions.cpp")
    })
  } else {
    parallel::clusterExport(cl, c("calcRunningSum", "calcMaxDeviation"))
  }
  
    # Export necessary variables and functions to the cluster
    parallel::clusterExport(cl, c("expr_matrix", "specificity_matrix", "sample_names",
                               "thresholds", "above_threshold_matrix", "cell_types",
                               "min_size", "p_weight", "nperm", "verbose", 
                               ".process_sample_optimized"),
                          envir = environment())
    
    # Process samples in parallel
    results <- parallel::parLapply(cl, 1:n_samples, function(sample_idx) {
      .process_sample_optimized(
        sample_idx = sample_idx,
        sample_names = sample_names,
        expr_matrix = expr_matrix,
        specificity_matrix = specificity_matrix,
        thresholds = thresholds,
        above_threshold_matrix = above_threshold_matrix,
        cell_types = cell_types,
        min_size = min_size,
        p_weight = p_weight,
        nperm = nperm,
        verbose = verbose
      )
    })
    names(results) <- sample_names
  } else {
    # Process samples sequentially
    results <- lapply(1:n_samples, function(sample_idx) {
      .process_sample_optimized(
        sample_idx = sample_idx,
        sample_names = sample_names,
        expr_matrix = expr_matrix,
        specificity_matrix = specificity_matrix,
        thresholds = thresholds,
        above_threshold_matrix = above_threshold_matrix,
        cell_types = cell_types,
        min_size = min_size,
        p_weight = p_weight,
        nperm = nperm,
        verbose = verbose
      )
    })
    names(results) <- sample_names
  }
  
  end_time <- Sys.time()
  if (verbose) {
    time_taken <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
    message(sprintf("Analysis completed in %s seconds", time_taken))
  }
  
  return(results)
}


#' Plot Cell Type Enrichment Results
#'
#' Creates a GSEA-style enrichment plot for a specified cell type from
#' specificityWeightedGSEA results.
#'
#' @param gsea_results List. Results from specificityWeightedGSEA function.
#' @param cell_type Character. Name of cell type to plot.
#' @param show_pval Logical. Whether to display p-value on plot. Default is TRUE.
#' @param title Character. Optional title for the plot.
#'
#' @return A ggplot2 object with the enrichment plot.
#'
#' @import ggplot2
#' @export
plotEnrichment <- function(gsea_results, cell_type, show_pval = TRUE, 
                          title = NULL) {
  # Check if required packages are installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is needed for this function")
  }
  
  # Input validation
  if (!cell_type %in% names(gsea_results[['fdr']])) {
    stop("Cell type not found in results")
  }
  
  # Get cell type index
  ct_idx <- which(names(gsea_results[['fdr']]) == cell_type)
  
  # Extract running scores
  run_scores <- gsea_results$running_scores[, ct_idx]
  
  # Create data frame for plotting
  plot_data <- data.table(
    position = 1:length(run_scores),
    running_score = run_scores
  )
  
  # Create annotation text
  if (show_pval) {
    pval <- gsea_results$p_values[ct_idx]
    padj <- gsea_results$fdr[ct_idx]
    pval_text <- paste("p-value =", format(pval, digits = 3),
                       "FDR =", format(padj, digits = 3))
    es <- gsea_results$enrichment_scores[ct_idx]
    nes <- gsea_results$normalized_scores[ct_idx]
    es_text <- paste("ES =", format(es, digits = 3), 
                    "NES =", format(nes, digits = 3))
  }
  
  # Set plot title
  if (is.null(title)) {
    plot_title <- paste("Enrichment plot:", cell_type)
  } else {
    plot_title <- title
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = position, y = running_score)) +
    ggplot2::geom_line(color = "#40004b", size = 1) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(
      xintercept = plot_data[running_score == es, position], 
      linetype = "dashed",
      color = "#1b7837"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = plot_title,
      x = "Protein Rank",
      y = "Enrichment Score"
    )
  
  # Add p-value annotation
  if (show_pval) {
    p <- p + ggplot2::annotate(
      "text", 
      x = max(plot_data$position) * 0.8, 
      y = max(abs(plot_data$running_score)) * 0.8,
      label = paste(pval_text, "\n", es_text)
    )
  }
  
  return(p)
}


#' Print multiple plots across pages in a grid layout
#'
#' This function takes a list of ggplot2 objects and arranges them in a 
#' specified grid layout per page, then prints them to a PDF file. By default, 
#' plots are arranged in a 2x2 grid (4 plots per page).
#'
#' @param plot_list A list of ggplot2 objects to be arranged and printed
#' @param plots_per_page Integer specifying number of plots per page (default: 4)
#' @param file Character string for the output PDF filename (default: "plots_by_page.pdf")
#' @param width PDF width in inches (default: 10)
#' @param height PDF height in inches (default: 8)
#' @param title Optional overall title for each page (default: NULL)
#' @param tag_levels Whether to add tags (labels) to plots. Can be 'A' for uppercase letters,
#'        'a' for lowercase letters, '1' for numbers, or NULL for no tags (default: NULL)
#'
#' @return Invisibly returns TRUE if successful
#'
#' @details
#' Currently only supports layouts with 4 plots per page in a 2x2 grid.
#' Empty plots are added automatically if the number of plots doesn't fill 
#' the last page completely.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' 
#' # Create a list of 10 plots
#' plot_list <- list()
#' for (i in 1:10) {
#'   plot_list[[i]] <- ggplot(mtcars, aes(hp, mpg)) + 
#'     geom_point() + 
#'     ggtitle(paste("Plot", i))
#' }
#' 
#' # Print plots 4 per page
#' print_plots_per_page(plot_list)
#' 
#' # With custom file name and page labels
#' print_plots_per_page(plot_list, 
#'                      file = "my_analysis_results.pdf", 
#'                      tag_levels = 'A')
#' }
#'
#' @import ggplot2
#' @import patchwork
#' @export
print_plots_per_page <- function(plot_list, 
                                 plots_per_page = 4, 
                                 file = "plots_by_page.pdf",
                                 width = 11, 
                                 height = 8.5,
                                 title = NULL,
                                 tag_levels = NULL) {
  
  # Check inputs
  if (!inherits(plot_list, "list")) {
    stop("plot_list must be a list of ggplot objects")
  }
  
  if (plots_per_page != 4) {
    stop("Currently only supports 4 plots per page")
  }
  
  if (length(plot_list) == 0) {
    warning("Empty plot list provided, no output generated")
    return(invisible(FALSE))
  }
  
  # Calculate pages needed
  n_plots <- length(plot_list)
  n_pages <- ceiling(n_plots / plots_per_page)
  
  # Open PDF device
  pdf(file, width = width, height = height)
  
  # For each page
  for (i in 1:n_pages) {
    # Get plot indices for this page
    start_idx <- (i-1) * plots_per_page + 1
    end_idx <- min(i * plots_per_page, n_plots)
    
    # Select plots for current page
    current_plots <- plot_list[start_idx:end_idx]
    
    # Fill with empty plots if needed to complete the grid
    if (length(current_plots) < plots_per_page) {
      empty_plots <- plots_per_page - length(current_plots)
      for (j in 1:empty_plots) {
        current_plots <- c(current_plots, list(ggplot() + theme_void()))
      }
    }
    
    # Create a 2x2 grid for this page
    layout <- (current_plots[[1]] + current_plots[[2]]) / 
              (current_plots[[3]] + current_plots[[4]])
    
    # Add annotations if requested
    if (!is.null(title) || !is.null(tag_levels)) {
      layout <- layout + patchwork::plot_annotation(
        title = title,
        tag_levels = tag_levels
      )
    }
    
    # Print the current page
    print(layout)
  }
  
  # Close the device
  dev.off()
  
  # Return success invisibly
  return(invisible(TRUE))
}


# Main -----------------------------------------------------------------------

# Stage data for analysis --------------------------
# DE file locations
fh <- file.path("baseline_PCs_relatedness", input,
  paste("dataProt_SS-205063.PC.related", input, "csv", sep = ".")
)
result <- lapply(fh, fread)
setattr(result, 'names', input)

# merge
result <- rbindlist(result)
input <- gsub(".APOE", "", input)

# filtered set of genes with at least one biomarker testing, top sig duplicate probe
all_genes <- unique(result[order(padj), .(Symbol, ix, padj)], by = 'Symbol')
all_genes <- all_genes[, .(GENE = unlist(tstrsplit(Symbol, "\\|"))), by = .(ix, padj)]
all_genes <- all_genes[!duplicated(GENE)]
# bio_genes <- unique(all_genes[padj < 0.05, .(Symbol, ix)], by = 'Symbol')
# fwrite(all_genes, file = "dataProt_SS-205063.probe_list.csv")

# get specificity index
spec <- fread(ref)
cols <- names(spec)[-1]
# spec <- as.matrix(spec[!GENE %in% bio_genes$Symbol, (cols) := 0], rownames = 'GENE')
spec <- as.matrix(spec, rownames = 'GENE')
# Replace NA values with 0
spec[is.na(spec)] <- 0

# get expression data and filter, rename, flip
prots <- fread("dataProt_SS-205063.merged.QC.csv")
cols <- c("SampleID", all_genes$ix)
prots <- prots[, ..cols]
setnames(prots, old = all_genes$ix, new = all_genes$Symbol)
prots <- transpose(prots, keep.names = 'Symbol', make.names = 'SampleID')
# prots <- as.matrix(prots[Symbol %in% bio_genes$Symbol], rownames = 'Symbol')
prots <- as.matrix(prots, rownames = 'Symbol')


# Build result table(s) & run function -------------
results <- specificityWeightedGSEA(
  expr_matrix = prots,
  specificity_matrix = spec,
  specificity_threshold = "median",
  nperm = 1000, 
  verbose = TRUE,
  n_cores = threads
)

# join and save files ------------------------------
saveRDS(results, file = file.path(output, paste0(nameset, ".scores.rds")))

if (plot) {
  # Set seed for reproducible sampling
  set.seed(123)
  sampled <- sample(results, size = 10)
  lapply(sampled, function(i) {
    plist <- lapply(names(i[['fdr']]), function(j) {
      plotEnrichment(i, j, title = paste(j, "in", i[['sample']]))
    })
    print_plots_per_page(plist, 
      file = file.path(output, paste(nameset, i[['sample']], "plots.pdf", sep = "."))
    )
  })
}

