#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 08a-Score_extraction.R
#'
#' Version: 0.0.2
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
    default = "celltype_scores",
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
    default = "celltype_scores",
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
#' @param dt A data.table containing nested list columns
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
#' \dontrun{
#' # Extract enrichment scores
#' scores_dt <- extract_nested_results(
#'   results_dt,
#'   id_col = "SampleID", 
#'   name_col = "cell_types", 
#'   value_col = "enrichment_scores"
#' )
#' 
#' # Extract p-values with different column names
#' pval_dt <- extract_nested_results(
#'   other_results,
#'   id_col = "PatientID", 
#'   name_col = "gene_sets", 
#'   value_col = "p_values"
#' )
#' }
#'
#' @export
extract_nested_results <- function(dt, id_col, name_col, value_col) {
  # Input validation
  if (!data.table::is.data.table(dt)) {
    stop("Input must be a data.table")
  }
  
  req_cols <- c(id_col, name_col, value_col)
  missing_cols <- setdiff(req_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Create a formula for dcast as string and then convert to formula
  formula_str <- paste(id_col, "~ name")
  cast_formula <- as.formula(formula_str)
  
  # Extract data in long format using data.table syntax
  long_result <- dt[, {
    # Get the names and values from specified columns
    names_list <- get(name_col)
    values_list <- get(value_col)
    
    # Handle both when these are plain vectors or nested lists
    if (is.list(names_list) && !is.data.frame(names_list)) {
      names_vec <- unlist(names_list)
    } else {
      names_vec <- names_list
    }
    
    if (is.list(values_list) && !is.data.frame(values_list)) {
      values_vec <- unlist(values_list)
    } else {
      values_vec <- values_list
    }
    
    # Check for valid data
    if (length(names_vec) != length(values_vec)) {
      warning(paste("Length mismatch between names and values for ID:", get(id_col)))
      return(data.table())
    }
    
    data.table(name = names_vec, value = values_vec)
  }, by = id_col]
  
  # If we have no valid data, return empty result with structure
  if (nrow(long_result) == 0) {
    return(data.table::data.table(setNames(list(character()), id_col)))
  }
  
  # Cast from long to wide format
  result_dt <- data.table::dcast(long_result, cast_formula, value.var = "value", fill = NA)
  
  return(result_dt)
}



# Main -----------------------------------------------------------------------
# results file location
result <- readRDS(file.path(inpath, paste0(nameset, ".scores.rds")))


# loop through requested stats
lapply(stat, function(i) {
  x <- extract_nested_results(result, 
    id_col = "SampleID",
    name_col = "cell_types", 
    value_col = i
  )
  
  fwrite(x, file = file.path(outpath, paste(nameset, i, "csv", sep = ".")))
})


