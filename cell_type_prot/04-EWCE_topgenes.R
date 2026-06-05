#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 04-EWCE_topgenes.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Show top cell-type proteins in EWCE association
#'
#' Author: Robert B Butler III
#'
#' Date Created: 2025-10-02
#'
#' Copyright (c) 2025
#' Email: rrbutler@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'     
#'   Usage:
#'     sbatch -J MG --mem=5G -c 2 -t 01:00:00 -p interactive \
#'       -o %x/%A_biodomain_enrichment_%x.log \
#'       --wrap "ml R/4.0; biodomain_enrichment.R -n PS19_C31 -g PS19 -d C31 -r R6 -t MG"
#'
#'   or interactive session:
#'     sdev -m 5 -c 2 -t 01:00:00 -p interactive
#'
#' ---------------------------

#' load up the packages we will need:  (uncomment as required)

library(optparse)
library(EWCE)
library(data.table)
library(readxl)
library(ggrepel)


#' ----------------------------

# default variables

# read in args
option_list <- list(
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project "
  ),
  make_option(
    c("-c", "--cset"),
    type = "character",
    help = "cell type dataset nickname (with path)"
  ),
  make_option(
    c("-g", "--gset"),
    type = "character",
    help = "geneset tested: c(all_genes, causal_genes, stringent_genes)"
  ),
  make_option(
    c("-t", "--trait"),
    type = "character",
    help = "trait to extract"
  ),
  make_option(
    c("-s", "--spec"),
    type = "character",
    default = "human",
    help = "species of ctd dataset [default: %default]"
  ),
  make_option(
    c("-l", "--level"),
    default = 1,
    help = "annotation level (1 or 2) [default: %default]"
  ),
  make_option(
    c("-r", "--cref"),
    type = "character",
    help = "cell type dataset name (with path)"
  ),
  make_option(
    c("-m", "--max"),
    default = 25,
    help = "Max.overlap for plot labels with ggrepel [default: %default]"
  ),
  make_option(
    c("-o", "--outpath"),
    type = "character",
    default = "Bootstrap_plots",
    help = "File path to output directory [default: %default]"
  ),
  make_option(
    c("-i", "--inpath"),
    type = "character",
    default = ".",
    help = "File path to input data [default: %default]"
  )
)
args_list = parse_args(OptionParser(option_list = option_list))

# read in fields
fields <- c(
  "nameset",
  "cset",
  "gset",
  "trait",
  "spec",
  "level",
  "cref",
  "max",
  "outpath",
  "inpath"
)
for (field in fields) {
  if (!field %in% names(args_list)) {
    cat(paste0('Argument "', field, '" must be provided\n'), fill = TRUE)
    stop(parse_args(OptionParser(option_list=option_list), args='--help'), 
         call.=FALSE)
  }
  
  assign(field, args_list[[field]])
  print(paste0(field, ': ', get(field)))
}

# make outpath
if(!dir.exists(outpath)) dir.create(outpath)

options(ggrepel.max.overlaps = max)

# Prep gene lists -------------------------------------------------------
# input variables (hardcoded for now)
input <- c(
  "PlasmaGFAP", "PlasmaNfL", "PlasmaPTau181", "PlasmapTau217", 
  "PlasmaAB140P", "PlasmaAB142P", "PlasmaABRatio", "PlasmapTau217_AB42Ratio"
)

# DE file locations
fh <- file.path("..", "baseline_PCs_relatedness", input,
  paste("dataProt_SS-205063.PC.related", input, "csv", sep = ".")
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

gene_sets <- c("stat_pval", "stat_padj", "stat_decile", "stat_quartile")
# stat_pval <- unstack(result[!is.na(stat_pval), .(Symbol, dx)])
lapply(gene_sets, function(i){
  assign(i, unstack(result[!is.na(get(i)), .(Symbol, dx)]), envir = .GlobalEnv)
})

# Plot key genes per significant associations ---------------------------------
load(cref)

# load EWCE results
fh <- paste(nameset, gset, cset, paste0("lvl", level), sep = ".")
res <- readRDS(file.path(inpath, paste(fh, "rds", sep = ".")))

# find trait in results list object
trait_ind <- match(trait, names(get(gset)))

# plot
bt <- generate_bootstrap_plots(
  sct_data = ctd,
  hits = get(gset)[[trait]],
  sctSpecies = spec,
  genelistSpecies = "human",
  reps = 10000,
  annotLevel = level,
  full_results = res[[trait_ind]],
  save_dir = outpath,
  listFileName = paste(fh, trait, sep = "."), 
  show_plot = FALSE
)

# save gene data
saveRDS(bt, 
  file = file.path(outpath, paste(fh, trait, "bootstrap_result.rds", sep = "."))
)

# # cleanup
# system("squeue -j $SLURM_JOBID --Format=TimeUsed")
