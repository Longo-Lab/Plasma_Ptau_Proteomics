#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 02-Combined-TREAT-AD_GSEA.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Pairwise GSEA correlations of TREAT-AD in proteomic associations 
#'
#' Author: Robert B Butler III, Michael Peter Brown
#'
#' Date Created: 2024-11-18
#'
#' Copyright (c) 2024
#' Email: rrbutler@stanford.edu, mpb@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'
#'
#' ---------------------------
#' load up the packages we will need:  (uncomment as required)
library(optparse)
library(ggplot2)
library(ggpubr)
library(data.table)
library(stringr)
library(plyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(GeneOverlap)
library(tools)

# Nameset helper function
nameset_default <- function(input){
  file_path_sans_ext(basename(input))
}

# read in args
option_list <- list(
  make_option(
    c("-m", "--modules"),
    type = "character",
    help = "Name of module script set desired to run"
  ),
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "List of input datasets comma separated"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = ".",
    help = "File path to output directory [default: %default]"
  ),
  make_option(
    c("-r", "--ref"),
    type = "character",
    default = "../biodomains",
    help = "File path to directory of reference data [default: %default]"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project [default: Nameset of input]"
  )
)

args_list <- parse_args(OptionParser(option_list = option_list))

# Nameset default case
if(is.null(args_list$nameset)){
  args_list$nameset <- nameset_default(args_list$input)
}

# read in fields
fields <- c(
  "modules",
  "input",
  "output",
  "ref",
  "nameset"
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

# Source() file path
source_file <- file.path(ref, paste0("biodomain_", modules, ".R"))

# Source() error/warning handling
tryCatch(
  expr = {
    source(source_file)
  },
  error = function(e){
    message('Error with source file:')
    print(e)
  },
  warning = function(w){
    message('Warning messages from source file:')
    print(w)
  }
)

# input variables
input <- strsplit(input, ",")[[1]]


# Functions ---------------------------

#' for a string of ens ids: explode, look up symbols, return symbols as string
#'
#' @param str string of comma separated ens_ids
#' @param result data.table with ens_ids in GENE col
#' @param col colname to grab (i.e. Symbol)
#'
#' @return concatenated string of lookup terms 
get_mcols <- function(str, result, col) {
  require("stringr")
  y <- result[[col]][match(unlist(str_split(str, ",")), result$GENE)]
  return(paste0(y, collapse = ","))
}



# Main ------------------------------------------
# DE file locations
fh <- file.path("../..", "baseline_PCs_relatedness", input,
  paste("dataProt_SS-205063.PC.related", input, "csv", sep = ".")
)
result <- lapply(fh, fread)
setattr(result, 'names', input)

# GSEA testing for enrichment terms from tstat ranking ----------------


# Pairwise Kendall's Taus (& Plot) ------------------------------------


# Possibly focus on sub-domains ---------------------------------------


# Driver domains/proteins ---------------------------------------------

