#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 05-EWCE_enrichment_WSS.R
##
## Version: 0.0.1
##
## Purpose of script: Show EWCE enrichments for WSS modules genesets
##
## Author: Robert B Butler III
##
## Date Created: 2026-03-12
##
## Copyright (c) 2026
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
    help = "WSS modules file to use"
  ),
  make_option(
    c("--inpath"),
    type = "character",
    default = "../module_scores/STRING_merged_scaled_scores.csv",
    help = "WSS modules file to use"
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
    c("-m", "--min-size"),
    dest = "min_size",
    default = 10,
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
  "min_size",
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

#' Filter Gene Sets Based on Reference Dataset
#'
#' Filters gene sets to retain only genes present in a reference dataset,
#' and removes modules that fall below a minimum size threshold after filtering.
#' Preserves the nested list structure of gene sets.
#'
#' @param gene_sets Character vector of gene set object names to filter.
#'   Each name should correspond to an object in the global environment
#'   containing a named list of gene modules (e.g., c("WSS_mod", "WSS_etc")).
#' @param ref_genes Character vector of gene symbols present in the reference
#'   dataset (e.g., rownames of expression matrix).
#' @param min_size Integer. Minimum number of genes required for a module to
#'   be retained after filtering. Default is 10.
#'
#' @return A named list where each element corresponds to a gene set from
#'   \code{gene_sets}. Each element contains a named list of filtered gene
#'   modules, where modules with fewer than \code{min_size} genes have been
#'   removed. The structure mirrors the input gene sets but with filtered
#'   gene lists.
#'
#' @details
#' The function performs the following steps for each gene set:
#' \enumerate{
#'   \item Retrieves the gene set object from the global environment
#'   \item Filters each module to retain only genes in \code{ref_genes}
#'   \item Removes modules with fewer than \code{min_size} genes
#'   \item Reports the number of modules retained
#' }
#'
#' Original gene set objects in the global environment are not modified.
#'
#' @examples
#' \dontrun{
#' # Load reference dataset
#' load("reference.rda")
#' ref_genes <- rownames(ctd[[1]]$mean_exp)
#'
#' # Create filtered gene sets
#' gene_sets <- c("WSS_mod", "stat_pval")
#' gene_sets_filt <- filter_gene_sets(gene_sets, ref_genes, min_size = 10)
#'
#' # Access filtered modules
#' gene_sets_filt$WSS_mod$Merged.M1
#' }
#'
#' @export
filter_gene_sets <- function(gene_sets, ref_genes, min_size = 10) {
  gene_sets_filt <- lapply(gene_sets, function(set_name) {
    # Get the current gene set from global environment
    current_set <- get(set_name, envir = .GlobalEnv)
    
    # Filter each module's genes to those in reference
    filtered_set <- lapply(current_set, function(module_genes) {
      module_genes[module_genes %in% ref_genes]
    })
    
    # Keep only modules with >= min_size genes after filtering
    filtered_set <- filtered_set[lengths(filtered_set) >= min_size]
    
    message(sprintf("%s: %d/%d modules retained after filtering (min_size=%d)", 
                    set_name, length(filtered_set), length(current_set), min_size))
    
    return(filtered_set)
  })
  
  # Name the list elements
  names(gene_sets_filt) <- gene_sets
  
  return(gene_sets_filt)
}


# Main ----------------------------------------------------------------
# DE file locations
result <- fread(inpath)

# filter and setup
# for genes in multiple modules, keep by top mean_alpha_specificity
result <- result[order(-mean_alpha_scaled)][, .SD[1], by = symbol]
result[, module := paste(Biomarker, cluster, sep = ".")]

# filter singles/doubles
keep <- result[, .N, by = module][N > 5, module]
result <- result[module %in% keep]

# just the ones specified above in input
result <- result[module %in% input]

# get module members as a named list
gene_sets <- c("WSS_mod")
# gene_sets <- c("stat_pval", "stat_padj", "stat_decile")
# stat_pval <- unstack(result[!is.na(stat_pval), .(Symbol, dx)])
lapply(gene_sets, function(i){
  gene_list <- unstack(result[, .(symbol, module)])
  assign(i, gene_list, envir = .GlobalEnv)
})

# Sharma et al 2017 proteomic celltypes data ----------
load(file.path(ref, "Sharma_etal.EWCE.human.rda"))

# Get genes filtered to reference genes
ref_genes <- rownames(ctd[[1]]$mean_exp)
gene_sets_filt <- filter_gene_sets(gene_sets, ref_genes, min_size)

# run them all in a loop at both levels
lapply(names(gene_sets_filt), function(i) {
  plot_ewce(
    x = ctd, 
    gene_sets = gene_sets_filt[[i]], 
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

# Get genes filtered to reference genes
ref_genes <- rownames(ctd[[1]]$mean_exp)
gene_sets_filt <- filter_gene_sets(gene_sets, ref_genes, min_size)

# run them all in a loop at both levels
lapply(names(gene_sets_filt), function(i) {
  lapply(seq(2), function(j) {
    plot_ewce(
      x = ctd, 
      gene_sets = gene_sets_filt[[i]], 
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














