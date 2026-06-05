#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 09-Eigengene_construction.R
#'
#' Version: 0.0.2
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
#'   v0.0.2 - pre-implementation of gene-scaled eigengenes
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
  result_positive <- result[mean_beta > 0]
  
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

#' #########NOT IMPLEMENTED Scale Module Eigengenes Back to Original Expression Scale
#'
#' Converts module eigengenes from standardized PCA scale back to the original 
#' expression scale by matching the mean and standard deviation of genes within 
#' each module.
#'
#' @param ME_result The moduleEigengenes result object containing eigengenes matrix
#' @param expr_matrix The original expression matrix (samples × genes)
#' @param color_vector Named vector assigning genes to module colors
#' @param method Character string specifying scaling method:
#'   "module_stats" (default) - match mean/SD of genes in the module
#'   "representative_gene" - match the gene with highest module membership
#'   "module_average" - use simple average instead of eigengene
#'
#' @return A matrix of scaled eigengenes in original expression units
#'
#' @details 
#' The "module_stats" method calculates the mean and standard deviation of all 
#' genes within each module across all samples, then scales the eigengene to 
#' have the same mean and standard deviation. This preserves the eigengene's 
#' patterns while returning it to interpretable expression units.
#'
#' @examples
#' # Calculate eigengenes
#' ME_result <- calculate_module_eigengenes_hybrid(expr_matrix, result)
#' 
#' # Scale back to original units
#' scaled_eigengenes <- scale_eigengenes_to_original(ME_result, expr_matrix, color_vector)
#' 
#' # Compare scales
#' summary(ME_result$eigengenes[,1])  # Standardized scale
#' summary(scaled_eigengenes[,1])     # Original expression scale
#'
#' @import data.table
#' @export
scale_eigengenes_to_original <- function(ME_result, expr_matrix, color_vector, 
                                        method = "module_stats") {
  
  # Extract eigengenes matrix
  eigengenes <- ME_result$eigengenes
  
  # Initialize scaled eigengenes matrix
  scaled_eigengenes <- eigengenes
  
  # Get unique module colors (excluding grey if desired)
  module_colors <- colnames(eigengenes)
  module_colors <- gsub("^ME", "", module_colors)  # Remove "ME" prefix
  
  for (module_color in module_colors) {
    
    # Skip grey module if it exists but has no specific genes assigned
    if (module_color == "grey") next
    
    # Get genes in this module
    genes_in_module <- names(color_vector)[color_vector == module_color]
    genes_in_expr <- intersect(genes_in_module, colnames(expr_matrix))
    
    if (length(genes_in_expr) == 0) {
      warning("No genes found for module ", module_color, " in expression matrix")
      next
    }
    
    # Current eigengene column
    eigengene_col <- paste0("ME", module_color)
    if (!eigengene_col %in% colnames(eigengenes)) {
      warning("Eigengene column ", eigengene_col, " not found")
      next
    }
    
    current_eigengene <- eigengenes[, eigengene_col]
    
    if (method == "module_stats") {
      # Calculate mean and SD of genes in this module
      module_expr <- expr_matrix[, genes_in_expr, drop = FALSE]
      
      # Calculate overall statistics across all genes in module
      all_values <- as.vector(module_expr)
      target_mean <- mean(all_values, na.rm = TRUE)
      target_sd <- sd(all_values, na.rm = TRUE)
      
      # Scale eigengene to match these statistics
      current_mean <- mean(current_eigengene, na.rm = TRUE)
      current_sd <- sd(current_eigengene, na.rm = TRUE)
      
      if (current_sd > 0) {
        scaled_eigengenes[, eigengene_col] <- 
          (current_eigengene - current_mean) * (target_sd / current_sd) + target_mean
      }
      
    } else if (method == "representative_gene") {
      # Find gene with highest correlation to eigengene
      correlations <- cor(current_eigengene, expr_matrix[, genes_in_expr, drop = FALSE], 
                         use = "pairwise.complete.obs")
      best_gene <- genes_in_expr[which.max(abs(correlations))]
      
      # Scale to match this gene's distribution
      target_values <- expr_matrix[, best_gene]
      target_mean <- mean(target_values, na.rm = TRUE)
      target_sd <- sd(target_values, na.rm = TRUE)
      
      current_mean <- mean(current_eigengene, na.rm = TRUE)
      current_sd <- sd(current_eigengene, na.rm = TRUE)
      
      if (current_sd > 0) {
        scaled_eigengenes[, eigengene_col] <- 
          (current_eigengene - current_mean) * (target_sd / current_sd) + target_mean
      }
      
    } else if (method == "module_average") {
      # Replace eigengene with simple average of module genes
      module_expr <- expr_matrix[, genes_in_expr, drop = FALSE]
      scaled_eigengenes[, eigengene_col] <- rowMeans(module_expr, na.rm = TRUE)
    }
  }
  
  return(scaled_eigengenes)
}

#' #########NOT IMPLEMENTED Calculate and Scale Module Eigengenes in One Step
#'
#' Combines eigengene calculation with scaling back to original expression units.
#'
#' @param expr_matrix Expression matrix (samples × genes)
#' @param result Cluster result table with Biomarker, cluster, ix, mean_beta columns
#' @param scale_method Method for scaling back to original units
#'
#' @return List containing ME_result and scaled_eigengenes
#'
#' @examples
#' # Calculate and scale eigengenes in one step
#' eigen_results <- calculate_and_scale_eigengenes(expr_matrix, result)
#' 
#' # Access results
#' original_eigengenes <- eigen_results$ME_result$eigengenes  # Standardized
#' scaled_eigengenes <- eigen_results$scaled_eigengenes       # Original scale
#'
#' @export
calculate_and_scale_eigengenes <- function(expr_matrix, result, scale_method = "module_stats") {
  
  # Calculate eigengenes with hybrid signed approach
  ME_result <- calculate_module_eigengenes_hybrid(expr_matrix, result)
  
  # Create color vector for scaling
  result[, module_color := paste(Biomarker, cluster, sep = "-")]
  result_positive <- result[mean_beta >= 0]
  
  gene_names <- colnames(expr_matrix)
  color_vector <- rep("grey", length(gene_names))
  names(color_vector) <- gene_names
  
  for(i in 1:nrow(result_positive)) {
    gene_id <- result_positive$ix[i]
    if(gene_id %in% gene_names) {
      color_vector[gene_id] <- result_positive$module_color[i]
    }
  }
  
  # Scale eigengenes back to original units
  scaled_eigengenes <- scale_eigengenes_to_original(ME_result, expr_matrix, 
                                                   color_vector, method = scale_method)
  
  return(list(
    ME_result = ME_result,
    scaled_eigengenes = scaled_eigengenes,
    color_vector = color_vector
  ))
}




# Main -----------------------------------------------------------------------
# Stage data for analysis --------------------------
result <- fread(input)

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
