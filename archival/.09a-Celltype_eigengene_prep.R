#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 09a-Cell-type_eigengene_prep.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Modify the cell type specificities for eigengene creation
#'  
#'
#' Author: Robert B Butler III
#'
#' Date Created: 2025-12-13
#'
#' Copyright (c) 2025
#' Email: rrbutler@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'   v0.0.1 - initial build
#'   
#'   
#'   
#'
#' ---------------------------

#' load up the packages we will need:  (uncomment as required)
library(data.table)
library(optparse)
library(WCluster)
library(RNOmni)


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
    default = "cell_type_prot",
    help = "File path to output directory [default: %default]"
  ),
  make_option(
    c("--INT"),
    action = "store_true",
    default = FALSE,
    help = "Inverse rank normalize expression [default: log2]"
  ),
  make_option(
    c("-s", "--spec"),
    type = "character",
    default = "cell_type_prot/Sharma_etal.EWCE.human.lvl1.spec_matr.txt",
    help = "Path to specificity matrix [default: %default]"
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
  "output",
  "INT",
  "threads",
  "prots",
  "spec",
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

# input variables
input <- strsplit(input, ",")[[1]]

# check outpath status
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}


# Functions -------------------------
#' Write Module Eigengenes Matrix to File
#'
#' Extracts the eigengenes matrix from moduleEigengenes results and writes it
#' to a tab-delimited file as-is.
#'
#' @param ME_result The moduleEigengenes result object containing the eigengenes matrix.
#' @param filename A character string specifying the output file path.
#'
#' @return Invisibly returns the eigengenes matrix that was written to file.
#'
#' @details 
#' The output file contains the eigengenes matrix exactly as returned by 
#' moduleEigengenes, with samples as rows and module eigengenes as columns.
#' The file is written in tab-delimited format with headers.
#'
#' @examples
#' # Calculate eigengenes
#' ME_result <- calculate_module_eigengenes_hybrid(expr_matrix, cluster_result)
#' 
#' # Write to file
#' write_eigengenes_matrix(ME_result, "module_eigengenes.txt")
#' 
#' # Get the matrix for further processing
#' eigengenes_matrix <- write_eigengenes_matrix(ME_result, "output.txt")
#'
#' @seealso 
#' \code{\link{calculate_module_eigengenes_hybrid}} for eigengene calculation
#'
#' @import data.table
#' @export
write_eigengenes_matrix <- function(ME_result, filename) {
  
  # Extract eigengenes matrix
  eigengenes_matrix <- ME_result$eigengenes
  
  # Write to file
  fwrite(
    x = as.data.table(eigengenes_matrix, keep.rownames = "SampleID"), 
    file = filename
  )
  
  # Return the matrix invisibly
  invisible(eigengenes_matrix)
}

#' Calculate Weighted Module Eigengenes
#'
#' Computes module eigengenes (first principal components) using weighted PCA,
#' where genes are weighted by their specificity scores for each module.
#' This allows genes with higher module-specific importance to contribute
#' more to the eigengene calculation.
#'
#' @param expr_data A numeric matrix of gene expression data with samples as rows
#'   and genes as columns. Row names should be sample IDs and column names
#'   should be gene IDs.
#' @param specificity_matrix A numeric matrix of specificity scores with genes as
#'   rows and modules as columns. Scores should typically range from 0 to 1,
#'   where higher values indicate greater specificity. Row names must be gene IDs
#'   matching those in \code{expr_data}.
#' @param min_weight Threshold for including genes in eigengene calculation.
#'   Either a numeric value (genes with specificity > min_weight are included)
#'   or the string "median" (genes above the median non-zero specificity for
#'   each module are included). Default is 0 (all genes with non-zero weights).
#' @param verbose Logical. If TRUE, prints progress information including the
#'   threshold and number of genes used for each module. Default is FALSE.
#' @param scale_power Numeric. Raises weights to this power to amplify differences.
#'   Values > 1 emphasize high-specificity genes more. Default is 1 (no scaling).
#'
#' @return A list with four elements:
#'   \describe{
#'     \item{eigengenes}{A data frame with samples as rows and module eigengenes
#'       as columns. Column names are formatted as "ME" + module name.}
#'     \item{varExplained}{A named numeric vector indicating the proportion of
#'       variance explained by the first principal component for each module.}
#'     \item{thresholds}{A named numeric vector showing the weight threshold
#'       used for each module.}
#'     \item{nGenesUsed}{A named numeric vector indicating the number of genes
#'       used in the eigengene calculation for each module.}
#'     \item{geneOverlap}{A matrix showing the number of shared genes between modules.}
#'   }
#'
#' @details
#' The function uses weighted PCA to calculate module eigengenes. For each module:
#' \enumerate{
#'   \item Genes are filtered based on \code{min_weight}
#'   \item Data is standardized (scaled to mean 0, sd 1)
#'   \item Genes are weighted by multiplying by their specificity scores
#'   \item SVD is performed to extract the first principal component
#'   \item The sign is adjusted so the eigengene positively correlates with
#'     module genes on average (following WGCNA convention)
#' }
#'
#' Modules with fewer than 2 genes passing the weight threshold are skipped
#' with a warning.
#'
#' @examples
#' \dontrun{
#' # Create example data
#' expr_data <- matrix(rnorm(1000), 100, 10)
#' colnames(expr_data) <- paste0("Gene", 1:10)
#' rownames(expr_data) <- paste0("Sample", 1:100)
#'
#' specificity_matrix <- matrix(runif(30), 10, 3)
#' rownames(specificity_matrix) <- paste0("Gene", 1:10)
#' colnames(specificity_matrix) <- c("Module1", "Module2", "Module3")
#'
#' # Calculate weighted eigengenes with default settings
#' result <- weighted_module_eigengenes(expr_data, specificity_matrix)
#'
#' # Use median threshold per module
#' result <- weighted_module_eigengenes(expr_data, specificity_matrix,
#'                                     min_weight = "median", verbose = TRUE)
#'
#' # Emphasize high-specificity genes more
#' result <- weighted_module_eigengenes(expr_data, specificity_matrix,
#'                                     min_weight = "median", scale_power = 2)
#'
#' # Access results
#' MEs <- result$eigengenes
#' var_exp <- result$varExplained
#' 
#' # Check module gene overlap
#' print(result$geneOverlap)
#' 
#' # Check correlation between eigengenes
#' cor(result$eigengenes)
#' }
#'
#' @seealso \code{\link[WGCNA]{moduleEigengenes}} for standard (unweighted)
#'   module eigengene calculation
#'
#' @export
weighted_module_eigengenes <- function(expr_data, specificity_matrix, 
                                      min_weight = 0, 
                                      verbose = FALSE,
                                      scale_power = 1) {
  
  # Input validation
  common_genes <- intersect(colnames(expr_data), rownames(specificity_matrix))
  
  if(length(common_genes) == 0) {
    stop("No common genes found between expr_data and specificity_matrix")
  }
  
  # Subset to common genes
  expr_data <- expr_data[, common_genes, drop = FALSE]
  specificity_matrix <- specificity_matrix[common_genes, , drop = FALSE]
  
  modules <- colnames(specificity_matrix)
  
  all_MEs <- list()
  var_exp <- numeric(length(modules))
  thresholds_used <- numeric(length(modules))
  genes_used <- numeric(length(modules))
  genes_per_module <- list()
  names(var_exp) <- modules
  names(thresholds_used) <- modules
  names(genes_used) <- modules
  
  for(i in seq_along(modules)) {
    mod <- modules[i]
    
    # Get specificity scores for this module
    weights <- specificity_matrix[, mod]
    
    # Determine threshold
    if(is.character(min_weight) && min_weight == "median") {
      # Use median of non-zero weights for this module
      non_zero_weights <- weights[weights > 0]
      if(length(non_zero_weights) == 0) {
        warning(sprintf("Module %s has no non-zero weights, skipping", mod))
        next
      }
      threshold <- median(non_zero_weights)
    } else if(is.numeric(min_weight)) {
      threshold <- min_weight
    } else {
      stop("min_weight must be numeric or 'median'")
    }
    
    thresholds_used[mod] <- threshold
    
    # Filter to genes with weight > threshold
    genes_to_use <- names(weights)[weights > threshold]
    genes_used[mod] <- length(genes_to_use)
    genes_per_module[[mod]] <- genes_to_use
    
    if(verbose) {
      cat(sprintf("Module %s: threshold=%.3f, genes=%d, weight range=[%.3f, %.3f]\n", 
                  mod, threshold, length(genes_to_use),
                  min(weights[genes_to_use]), max(weights[genes_to_use])))
    }
    
    if(length(genes_to_use) < 2) {
      warning(sprintf("Module %s has < 2 genes with weight > %.3f, skipping", 
                      mod, threshold))
      next
    }
    
    # Subset expression data and weights
    module_expr <- expr_data[, genes_to_use, drop = FALSE]
    module_weights <- weights[genes_to_use]
    
    # Apply power scaling to amplify weight differences
    if(scale_power != 1) {
      module_weights <- module_weights^scale_power
    }
    
    # Manual weighted PCA using SVD
    # 1. Scale the data (standardize each gene)
    module_expr_scaled <- scale(module_expr)
    
    # 2. Weight genes directly
    module_expr_weighted <- t(t(module_expr_scaled) * module_weights)
    
    # 3. Perform SVD
    svd_result <- svd(module_expr_weighted)
    
    # 4. Get PC1 scores (left singular vectors * first singular value)
    eigengene <- svd_result$u[, 1] * svd_result$d[1]
    
    # Calculate variance explained
    total_var <- sum(svd_result$d^2)
    var_explained <- svd_result$d[1]^2 / total_var
    
    # Adjust sign so average correlation is positive
    if(mean(cor(module_expr, eigengene)) < 0) {
      eigengene <- -eigengene
    }
    
    # Format output
    ME_df <- data.frame(eigengene)
    colnames(ME_df) <- paste0("ME", mod)
    rownames(ME_df) <- rownames(module_expr)
    
    all_MEs[[i]] <- ME_df
    var_exp[mod] <- var_explained
  }
  
  # Combine all eigengenes
  MEs <- do.call(cbind, all_MEs)
  
  # Calculate gene overlap between modules
  overlap_matrix <- matrix(0, length(modules), length(modules),
                          dimnames = list(modules, modules))
  for(i in seq_along(modules)) {
    for(j in seq_along(modules)) {
      if(modules[i] %in% names(genes_per_module) && 
         modules[j] %in% names(genes_per_module)) {
        overlap_matrix[i, j] <- length(intersect(genes_per_module[[modules[i]]], 
                                                  genes_per_module[[modules[j]]]))
      }
    }
  }
  
  return(list(
    eigengenes = MEs,
    varExplained = var_exp,
    thresholds = thresholds_used,
    nGenesUsed = genes_used,
    geneOverlap = overlap_matrix
  ))
}


# Main -----------------------------------------------------------------------
# need 
# (1) wide to long flip 
# (2) mean beta by gene x cluster
# (3) zero out and drop unincluded (non-zero median specificity)
# Stage data for analysis --------------------------
# DE file locations
fh <- file.path("baseline_PCs_relatedness", input,
  paste("dataProt_SS-205063.PC.related", input, "csv", sep = ".")
)
result <- lapply(fh, fread)
setattr(result, 'names', input)

# merge
result <- rbindlist(result)
result

# filtered set of genes with at least one biomarker testing
prot <- fread(prots)

# calculate mean beta
prot <- merge(prot, result[, .(mean_beta = mean(beta)), by = ix], all.x = TRUE)
prot[, Biomarker := "celltype"]

# Load specificity scores and trim to non-zero median -------------------------
cells <- fread(spec)
scores <- melt(cells, 
  id.vars = "GENE", 
  variable.name = "cluster", 
  value.name = "specificity"
)
scores <- merge(scores, prot)
setnames(scores, "GENE", "symbol")
scores[, term := cluster] # at this point this is a placeholder for Celltypes to be labeled correctly in downstream plots

# # STRING db Symbol error, labeled CYRIB as FAM49B
# scores[symbol == "FAM49B", symbol := "CYRIB"]
# # This one is on the SomaLogic folks (should be COXFA4L3)
# scores[symbol == "C15orf48", symbol := "NMES1"]

# Writing scores for standard eigengene quantiles -----------------------------
cols <- c("Biomarker", "symbol", "cluster", "ix", "mean_beta", "specificity", "term")
quants <- c(0.5, 0.6, 0.7, 0.8, 0.9)

lapply(quants, function(i) {
  fwrite(scores[specificity != 0][, .SD[specificity >= quantile(specificity, i)], by = cluster], 
    file = file.path(output, paste(nameset, i, "quant_scores.csv", sep = "_"))
  )
})


# # Calculating specificity-weighted PCA eigengenes -----------------------------
# # load in expression matrix
# # drop outlier samples
# a <- fread(ref)
# drop = a[abs(ConnectivityZscore) > 10, SampleID]
# a = a[!SampleID %in% drop]

# # Post filter simple table
# cols <- c(2, 28:ncol(a))
# b <- merge(
  # x = prot[, .(ix, GENE)], 
  # y = data.table::transpose(a[, ..cols], make.names = "SampleID", keep.names = "ix"), 
  # all.x = TRUE
# )
# exp <- as.matrix(
  # x = data.table::transpose(b[, -1], make.names = "GENE", keep.names = "SampleID"), 
  # rownames = "SampleID"
# )

# if (INT) {
  # # Inverse rank normalization (INT) to z-scale
  # exp <- apply(exp, 2, RankNorm)
  # nameset <- paste(nameset, "INT", sep = ".")
# } else {
  # # Log transformation
  # exp <- log2(exp)
# }

# # Calculate eigengenes with hybrid signed approach
# ME_result <- weighted_module_eigengenes(exp, 
  # as.matrix(cells, rownames = "GENE"), 
  # min_weight = "median",
  # scale_power = 2
# )

# # Write eigengenes matrix to file
# write_eigengenes_matrix(ME_result, 
  # filename = file.path(output, paste(nameset, "weighted_eigengenes.csv", sep = "."))
# )
# saveRDS(ME_result,   
  # file = file.path(output, paste(nameset, "weighted_eigengenes.rds", sep = "."))
# )

# # # Access components if needed
# # eigengenes <- ME_result$eigengenes
# # valid_colors <- ME_result$validColors
# # all_ok <- ME_result$allOK




