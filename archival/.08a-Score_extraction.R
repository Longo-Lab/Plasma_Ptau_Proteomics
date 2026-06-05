#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 08a-Score_extraction.R
#'
#' Version: 0.0.3
#'
#' Purpose of script: Get celltype matricies from results for Score celltype 
#'  
#'
#' Author: Robert B Butler III
#'
#' Date Created: 2025-10-03
#'
#' Copyright (c) 2025
#' Email: rrbutler@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'   Updates:
#'   v0.0.3 - updating for c version of 8 to handle new data structure
#'   v0.0.2 - bugfixes for outpath typos and help menu
#'   
#'   
#'
#' ---------------------------

#' load up the packages we will need:  (uncomment as required)
library(data.table)
library(optparse)

# Nameset helper function
nameset_default <- function(input){
  file_path_sans_ext(basename(input))
}

# read in args
option_list <- list(
  make_option(
    c("-i", "--inpath"),
    type = "character",
    default = "module_scores",
    help = "directory with results rds [default: %default]"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name prefix of the *scores.rds project file"
  ),
  make_option(
    c("-s", "--statistic"),
    dest = "stat",
    type = "character",
    help = paste(
      "statistic(s) to retrieve in flat file [stat1,stat2,...] options are:\n",
      "enrichment_scores, normalized_scores, p_values, fdr_values, thresholds used"
    )
  ),
  make_option(
    c("-o", "--outpath"),
    type = "character",
    default = "module_scores",
    help = "File path to output directory [default: %default]"
  )
)

args_list <- parse_args(OptionParser(option_list = option_list))

# Nameset default case
if(is.null(args_list$nameset)){
  args_list$nameset <- nameset_default(args_list$input)
}

# read in fields
fields <- c(
  "inpath",
  "outpath",
  "nameset",
  "stat"
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

# check outpath status
if (!dir.exists(outpath)) {
  dir.create(outpath, recursive = TRUE)
}

# process stat to list
stat <- strsplit(stat, ",")[[1]]


# Functions -------------------------
#' Extract nested list values from a data.table and reshape to wide format
#'
#' This function extracts values from a nested list column in a data.table and
#' reshapes the data to a wide format where the names become column headers.
#'
#' @param results A results object containing listed samples
#' @param id_col Character string specifying the ID column name (e.g., "SampleID")
#' @param name_col Character string specifying the column containing names that will
#'   become column headers in the result (e.g., "cell_types")
#' @param value_col Character string specifying the column containing values to extract
#'
#' @return A data.table in wide format where:
#'   - Each row is identified by the id_col
#'   - Each column corresponds to a unique name from name_col
#'   - The values come from the value_col
#'   - Any missing combinations are filled with NA
#'
#' @importFrom data.table data.table dcast :=
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
#' @export
extract_simple_results <- function(x, id_col, name_col, value_col) {
}



# Main -----------------------------------------------------------------------
# results file location
results <- readRDS(file.path(inpath, paste0(nameset, ".scores.rds")))

# loop through requested stats
lapply(stat, function(i) {
  fwrite(
    x = data.table(
      t(sapply(results, function(x) x[[i]])), 
      keep.rownames = "SampleID"
    ),
    file = file.path(outpath, paste(nameset, i, "csv", sep = "."))
  )
})


