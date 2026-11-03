#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 03a-GNPC_imputation.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Impute missing values in GNPC using kNN
#'  
#'
#' Author: Robert B Butler III
#'
#' Date Created: 2025-11-16
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
library(foreach)
library(doParallel)
library(VIM)


# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "GNPC_table.csv",
    help = "Data file to be imputed"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = ".",
    help = "File path to output directory [default: %default]"
  ),
  make_option(
    c("-n", "--nameset"),
    default = "GNPC_table.imputed",
    type = "character",
    help = "Name of project [default: Nameset of input]"
  ),
  make_option(
    c("-t", "--threads"),
    default = 2,
    help = "num threads [default: parallelly::availableCores()]"
  )
)

args_list <- parse_args(OptionParser(option_list = option_list))

# read in fields
fields <- c(
  "input",
  "nameset",
  "output",
  "threads"
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
knn_imputation <- function(dt, k = 5, gene_chunk_size = 2000, n_cores = NULL) {

  # Set up parallel backend
  if(is.null(n_cores)) {
    n_cores <- detectCores() - 1
  }
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  cat("Transposing dataset...\n")

  # Single transpose: genes as rows → samples as rows
  mat <- as.matrix(
    x = transpose(dt, keep.names = 'Samp', make.names = 'GENE'), 
    rownames = 'Samp'
  )
  print(cat("Rownames:", head(rownames(mat)), sep = " "))

  # Calculate chunks
  total_genes <- ncol(mat)
  n_chunks <- ceiling(total_genes / gene_chunk_size)
  cat("Processing", n_chunks, "gene chunks of size", gene_chunk_size, "using", n_cores, "cores\n")

  # Create chunk indices
  chunk_indices <- lapply(1:n_chunks, function(i) {
    start_col <- (i - 1) * gene_chunk_size + 1
    end_col <- min(i * gene_chunk_size, total_genes)
    start_col:end_col
  })

  # Process chunks in parallel
  results <- foreach(indices = chunk_indices,
                     .packages = c("VIM")) %dopar% {

    # Extract gene chunk
    chunk_mat <- mat[, indices, drop = FALSE]
    original_genes_chunk <- colnames(chunk_mat)

    # kNN imputation on chunk
    sink(nullfile())
    imputed_chunk <- VIM::kNN(chunk_mat, k = k)
    sink()
    
    # Return only original columns (remove VIM indicator columns)
    imputed_chunk[, original_genes_chunk, drop = FALSE]
  }

  # Reduce overhead
  sample_names <- rownames(mat)
  remove(mat)
  gc()

  # Stop cluster
  stopCluster(cl)

  # Combine chunks back into full matrix
  cat("Combining chunks...\n")
  imputed <- do.call(cbind, results)
  rownames(imputed) <- sample_names
  remove(results)
  gc()
  print(cat("Rownames:", head(rownames(imputed)), sep = " "))

  # Transpose back to original structure
  cat("Transposing back to original structure...\n")
  return(
    transpose(as.data.table(imputed, keep.rownames = TRUE), 
      keep.names = "GENE", 
      make.names = 'rn'
    )
  )
}

analyze_missing_values <- function(dt, top_n = 20) {
  # Calculate missing proportion for each protein (row)
  tot <- ncol(dt) - 1
  protein_missing <- dt[, .(missingness = sum(is.na(.SD)) / tot), by = GENE]
  setorder(protein_missing, -missingness)

  # Overall statistics
  overall_missing <- sum(is.na(dt)) / (nrow(dt) * ncol(dt))

  cat("=== PROTEIN MISSING VALUE ANALYSIS ===\n")
  cat(sprintf("Overall missing proportion: %.2f%%\n", overall_missing * 100))
  cat(sprintf("Total proteins: %d\n", nrow(dt)))
  cat(sprintf("Total samples: %d\n", ncol(dt)))

  # Top proteins with highest missingness
  cat(sprintf("\nTop %d proteins with highest missingness:\n", top_n))
  
  for(i in 1:top_n) {
    cat(sprintf("%2d. %s: %.1f%%\n", i, protein_missing[i, GENE], protein_missing[i, missingness] * 100))
  }

  # Summary statistics for proteins
  cat(sprintf("\nProtein missingness summary:\n"))
  cat(sprintf("Median protein missing: %.1f%%\n", median(protein_missing$missingness) * 100))
  cat(sprintf("Proteins with >50%% missing samples: %d\n", sum(protein_missing$missingness > 0.5)))
  cat(sprintf("Proteins with >80%% missing samples: %d\n", sum(protein_missing$missingness > 0.8)))

  return(protein_missing)
}

# Main -----------------------------------------------------------------------
# Stage data for analysis 
dt <- fread(input)
dt[, column_name := NULL]

# Convert -1 to NA
dt[dt == -1] <- NA

# Analyze missing values
# missing_analysis <- analyze_missing_values(dt, top_n = 20)

# # Get proteins with extreme missingness if you want to remove them
# extreme_missing <- missing_analysis[missing_analysis > 0.9]
# cat(sprintf("\nProteins with >90%% missing (%d proteins):\n", length(extreme_missing)))
# print(names(extreme_missing))

# impute
dt_imputed <- knn_imputation(dt, gene_chunk_size = 3194, n_cores = threads)
fwrite(dt_imputed, file = file.path(output, paste(nameset, "csv", sep = ".")))



