#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 03-EWCE_enrichment.R
##
## Version: 0.0.2
##
## Purpose of script: Show EWCE enrichments for profiles with biomarker genesets
##
## Author: Robert B Butler III
##
## Date Created: 2025-09-22
##
## Copyright (c) 2025
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##   v0.0.2 - updated to handle flexible filepaths/nameset/suffix
##   Currently hardcodes the input files (good target for options funtion)
##   paper DOI: none, need to publish
##
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(data.table)
library(ggplot2)
library(stringr)
library(optparse)
library(Cairo)
library(RColorBrewer)
library(patchwork)
library(readxl)
library(EWCE)
future::plan('multicore')

# input variables
fdr_cutoff = 0.05
# threads = 0


# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Biomarker or group names list to analyze"
  ),
  make_option(
    c("--inpath"),
    type = "character",
    default = "../baseline_PCs_relatedness",
    help = "File path to DE input file [default: %default]"
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
    default = ".",
    help = "File path to directory of reference data [default: %default]"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    default = "dataProt_SS-205063.PC.related",
    help = "filename prefix of input project [default: %default]"
  ),
  make_option(
    c("-s", "--suffix"),
    type = "character",
    default = "all_bios",
    help = "output suffix to add to nameset [default: %default]"
  ),
  make_option(
    c("-t", "--threads"),
    default = 0,
    help = "num threads [default: parallelly::availableCores()]"
   )
)

args_list <- parse_args(OptionParser(option_list = option_list))


# read in fields
fields <- c(
  "input",
  "inpath", 
  "output",
  "ref",
  "nameset",
  "suffix",
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
  threads <- as.integer(parallelly::availableCores())
} else if (is.numeric(threads)) {
  threads <- as.integer(threads)
} else {
  cat("threads must be an integer value")
  stop(parse_args(OptionParser(option_list = option_list), args = "--help"),
    call. = FALSE
  )
}
print(sprintf("Setting threads to %d", threads))


# input variables
input <- strsplit(input, ",")[[1]]


# Functions ------------------------------------------------------------
#' Run the function for generating EWCE results and plotting
#'
#' @param x ctd dataset
#' @param gene_sets named list, gene sets to run with annotation as name
#' @param set_name string, name of gene_sets list
#' @param ref_name string, name of single-cell reference
#' @param reps integer, number of permutations to run
#' @param annot integer, annotation level to use: 1 or 2
#' @param spec string, reference species
#'
#' @return none, saves full results file and images
plot_ewce <- function(x, 
                      gene_sets, 
                      set_name, 
                      ref_name, 
                      reps = 10000, 
                      annot = 1, 
                      w = 7, 
                      h = 12,
                      spec = "mouse",
                      thread = threads
  ) {

  outname <- paste(
    nameset, suffix, set_name, ref_name, paste0("lvl", annot), 
    sep = "."
  )
  full_results <- lapply(names(gene_sets), function(i) {
    bootstrap_enrichment_test(
      sct_data = x,
      sctSpecies = spec,
      genelistSpecies = "human",
      hits = gene_sets[[i]], 
      reps = reps,
      annotLevel = annot,
      no_cores = thread
    )
  })
  setattr(full_results, 'names', names(gene_sets))
  saveRDS(full_results, file = paste(outname, "rds", sep = "."))

  plot_results <- lapply(seq(gene_sets), function(i) {
    data.frame(full_results[[i]][['results']], list = names(gene_sets)[i])
  })
  plot_results <- rbindlist(plot_results)

  try({
    plot_list <- ewce_plot(total_res = plot_results, mtc_method = "BH")
    ggsave(
      plot = plot_list$plain, 
      filename = paste(outname, "pdf", sep = "."),
      width = w * annot,
      height = h,
      device = cairo_pdf
    )
  })
}

# Main ----------------------------------------------------------------
# DE file locations
fh <- file.path(inpath, input,
  paste(nameset, input, "csv", sep = ".")
)
result <- lapply(fh, fread)
setattr(result, 'names', input)

# merge
result <- rbindlist(result)
input <- gsub(".APOE", "", input)

# thresholds of cutoffs, pval, padj, outer_decile, outer_quartile
result[, stat_pval := ifelse(pval < 0.05, `t stat`, NA)]
result[, stat_padj := ifelse(padj < 0.05, `t stat`, NA)]
lapply(input, function(i) {
  qs <- quantile(result[dx == i, `t stat`], c(0.1, 0.25, 0.75, 0.9))
  result[, stat_decile := ifelse(`t stat` < qs[[1]] | `t stat` > qs[[4]], `t stat`, NA)]
  result[, stat_quartile := ifelse(`t stat` < qs[[2]] | `t stat` > qs[[3]], `t stat`, NA)]
})

gene_sets <- c("stat_padj", "stat_decile")
# gene_sets <- c("stat_pval", "stat_padj", "stat_decile")
# stat_pval <- unstack(result[!is.na(stat_pval), .(Symbol, dx)])
lapply(gene_sets, function(i){
  assign(i, unstack(result[!is.na(get(i)), .(Symbol, dx)]), envir = .GlobalEnv)
})

# Sharma et al 2017 proteomic celltypes data ----------
load(file.path(ref, "Sharma_etal.EWCE.human.rda"))

# run them all in a loop at both levels
lapply(gene_sets, function(i) {
  plot_ewce(
    x = ctd, 
    gene_sets = get(i), 
    set_name = i, 
    ref_name = "Sharma_etal.cultured", 
    annot = 1,
    spec = "human"
  )
})


# # Allen Inst SEA-AD data --------------
# # get overall downsampled SCE converted to ctd format
# load(file.path(ref, "ewce_SEA-AD.Class.1kDS.rda"))

# # run them all in a loop at both levels
# lapply(gene_sets, function(i) {
  # lapply(seq(2), function(j) {
    # plot_ewce(
      # x = ctd, 
      # gene_sets = get(i), 
      # set_name = i, 
      # ref_name = "SEA-AD.Class", 
      # annot = j, 
      # spec = "human"
    # )
  # })
# })

# get overall downsampled SCE converted to ctd format
load(file.path(ref, "ewce_SEA-AD.Neighborhood.1kDS.rda"))

# run them all in a loop at both levels
lapply(gene_sets, function(i) {
  lapply(seq(2), function(j) {
    plot_ewce(
      x = ctd, 
      gene_sets = get(i), 
      set_name = i, 
      ref_name = "SEA-AD.Neighborhood", 
      annot = j, 
      spec = "human"
    )
  })
})

# # get GABA downsampled
# load(file.path(ref, "ewce_SEA-AD.GABA.500DS.rda"))

# # run them all in a loop at both levels
# lapply(gene_sets, function(i) {
  # lapply(seq(2), function(j) {
    # plot_ewce(
      # x = ctd, 
      # gene_sets = get(i), 
      # set_name = i, 
      # ref_name = "SEA-AD.GABA", 
      # annot = j, 
      # spec = "human"
    # )
  # })
# })

# # get Gluta downsampled
# load(file.path(ref, "ewce_SEA-AD.Glut.500DS.rda"))

# # run them all in a loop at both levels
# lapply(gene_sets, function(i) {
  # lapply(seq(2), function(j) {
    # plot_ewce(
      # x = ctd, 
      # gene_sets = get(i), 
      # set_name = i, 
      # ref_name = "SEA-AD.Glut", 
      # annot = j, 
      # spec = "human"
    # )
  # })
# })

# # get Non-neuron downsampled
# load(file.path(ref, "ewce_SEA-AD.Non-neuronal.500DS.rda"))

# # run them all in a loop at both levels
# lapply(gene_sets, function(i) {
  # lapply(seq(2), function(j) {
    # plot_ewce(
      # x = ctd, 
      # gene_sets = get(i), 
      # set_name = i, 
      # ref_name = "SEA-AD.Non-neuron", 
      # annot = j, 
      # spec = "human"
    # )
  # })
# })














