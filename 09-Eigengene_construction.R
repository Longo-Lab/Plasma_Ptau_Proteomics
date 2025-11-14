#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 09-Eigengene_construction.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Generate eigengenes for STRING network modules
#'  
#'
#' Author: Robert B Butler III
#'
#' Date Created: 2025-11-12
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
library(WGCNA)
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
    default = "cell_type_prot/STRING_merged_scaled_scores.csv",
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
    help = "Inverse rank normalize expression [default: log2]"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "celltype_scores",
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
#' Calculate Module Eigengenes Using Hybrid Signed Approach
#'
#' This function calculates module eigengenes for predefined clusters using a 
#' hybrid signed network approach, where only genes with positive mean_beta 
#' (positive trait correlation) are included in module eigengene calculation.
#' Genes with negative correlations are excluded from modules and assigned to grey.
#'
#' @param expr_matrix A numeric matrix with samples as rows and genes/proteins 
#'   as columns. Column names should match the 'symbol' values in the result table.
#' @param result A data.table containing cluster information with required columns:
#'   Biomarker, cluster, symbol (gene identifiers), mean_beta (trait correlations).
#'
#' @return The complete moduleEigengenes result object from WGCNA containing
#'   eigengenes matrix, validColors, allOK, and other WGCNA output components.
#'
#' @details 
#' The hybrid signed approach excludes genes with negative trait correlations
#' (mean_beta < 0) from module eigengene calculation. These genes are assigned
#' to the "grey" module and effectively excluded from downstream analyses.
#' 
#' Module colors are created by concatenating Biomarker and cluster with "-".
#' Only genes with mean_beta >= 0 are assigned to their respective modules.
#'
#' @examples
#' # Calculate eigengenes with hybrid signed approach
#' ME_result <- calculate_module_eigengenes_hybrid(expr_matrix, cluster_result)
#' 
#' # Access eigengenes matrix
#' eigengenes <- ME_result$eigengenes
#' 
#' # Check valid colors
#' valid_colors <- ME_result$validColors
#'
#' @seealso 
#' \code{\link[WGCNA]{moduleEigengenes}} for the underlying eigengene calculation
#' \code{\link{write_eigengenes_matrix}} for writing results to file
#'
#' @import data.table
#' @import WGCNA
#' @export
calculate_module_eigengenes_hybrid <- function(expr_matrix, result) {
  
  # Create module colors from Biomarker-cluster combinations
  result[, module_color := paste(Biomarker, cluster, sep = "-")]
  
  # For hybrid signed: keep only genes with positive mean_beta
  result_positive <- result[mean_beta >= 0]
  
  # Get gene names from expression matrix
  gene_names <- colnames(expr_matrix)
  
  # Initialize all genes as "grey" (unassigned/excluded)
  color_vector <- rep("grey", length(gene_names))
  names(color_vector) <- gene_names
  
  # Assign module colors only to genes with positive mean_beta
  for(i in 1:nrow(result_positive)) {
    gene_id <- result_positive$symbol[i]
    if(gene_id %in% gene_names) {
      color_vector[gene_id] <- result_positive$module_color[i]
    }
  }
  
  # Calculate and return module eigengenes
  moduleEigengenes(
    expr = expr_matrix,
    colors = color_vector,
    nPC = 1,
    align = "along average",
    excludeGrey = FALSE,
    grey = "grey"
  )
}

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


# Main -----------------------------------------------------------------------
# Stage data for analysis --------------------------
result <- fread(input)

# split string and explode piped terms into separate rows
prot <- fread(prots)
prot <- prot[, .(GENE = unlist(tstrsplit(Symbol, "\\|"))), by = .(ix)]

# drop outlier samples
a <- fread(ref)
drop = a[abs(ConnectivityZscore) > 10, SampleID]
a = a[!SampleID %in% drop]

# Post filter simple table
cols <- c(2, 28:ncol(a))
b <- merge(
  x = prot, 
  y = transpose(a[, ..cols], make.names = "SampleID", keep.names = "ix"), 
  all.x = TRUE
)
exp <- as.matrix(
  x = transpose(b[, -1], make.names = "GENE", keep.names = "SampleID"), 
  rownames = "SampleID"
)

if (INT) {
  # Inverse rank normalization (INT) to z-scale
  exp <- apply(exp, 2, RankNorm)
  nameset <- paste(nameset, "INT", sep = ".")
} else {
  # Log transformation
  exp <- log2(exp)
}

# Run the analysis ----------------------------------
# Calculate eigengenes with hybrid signed approach
ME_result <- calculate_module_eigengenes_hybrid(exp, result)

# Write eigengenes matrix to file
write_eigengenes_matrix(ME_result, 
  filename = file.path(output, paste(nameset, "eigengenes.csv", sep = "."))
)
saveRDS(ME_result,   
  file = file.path(output, paste(nameset, "eigengenes.rds", sep = "."))
)

# # Access components if needed
# eigengenes <- ME_result$eigengenes
# valid_colors <- ME_result$validColors
# all_ok <- ME_result$allOK
