#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 03-Long_plots.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Plots of biomarker characteristics in Fig 1
#'
#' Author: Robert B Butler III
#'
#' Date Created: 2025-02-24
#'
#' Copyright (c) 2025
#' Email: rrbutler@stanford.edu
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
library(RColorBrewer)
library(data.table)
library(patchwork)

# read in args
option_list <- list(
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
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project [default: Nameset of input]"
  )
)

args_list <- parse_args(OptionParser(option_list = option_list))

# Nameset default case
if(is.null(args_list$nameset)){
  args_list$nameset <- file_path_sans_ext(basename(args_list$input))
}

# read in fields
fields <- c(
  "input",
  "output",
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

# input variables
input <- strsplit(input, ",")[[1]]


# Functions ---------------------------


# Main ------------------------------------------
# DE file locations ########################just need to load overall and filter samples the same.
fh <- file.path("../..", "baseline_PCs_relatedness", input,
  paste("dataProt_SS-205063.PC.related", input, "csv", sep = ".")
)
result <- lapply(fh, fread)
setattr(result, 'names', input)

# merge
result <- rbindlist(result)

########################check the n for each biomarker match up to verify


######plot vs Age (consider MOCA) scatter w/trends
######plot vs diagnosis boxplot

p <- ggpairs(
  data = dcast(result, NCBI_ID ~ dx, value.var = "t stat"),
  upper = list(continuous = upper_fn),
  diag = list(continuous = function(...) ggally_densityDiag(...) + theme_linedraw()),
  lower = list(continuous = lower_fn),
  columns = 2:9
) 

ggsave(p,
  filename = file.path(output, paste(nameset, "correlations.pdf", sep=".")), 
  device = cairo_pdf,
  width = 12, 
  height = 12
)
ggsave(p,
  filename = file.path(output, paste(nameset, "correlations.png", sep=".")), 
  width = 12, 
  height = 12,
  units = "in", 
  dpi = 144, 
  type = "cairo"
)
