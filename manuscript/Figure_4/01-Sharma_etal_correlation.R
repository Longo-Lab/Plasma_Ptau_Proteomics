#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 04-TREAT-AD_correlation.R
##
## Version: 0.0.2
##
## Purpose of script: Show Correlation of Sharma_etal proteomic celltype
##  profiles with PlasmaPtau181 geneset
##
## Author: Robert B Butler III, Michael P Brown
##
## Date Created: 2022-12-21
##
## Copyright (c) 2022
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##   Currently hardcodes the input files (good target for options funtion)
##   paper DOI: none, need to publish
##
##   Using T stat of Ptau with Pro_TE for TREAT-AD
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(data.table)
# library(corrplot)
library(stringr)
library(optparse)
library(Cairo)
library(RColorBrewer)
library(ggcorrplot)
library(ggh4x)
library(patchwork)

fdr_cutoff = 0.05

# Nameset helper function
nameset_default <- function(input){
  file_path_sans_ext(basename(input))
}

# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = ".",
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
    default = "../cell_type_prot",
    help = "File path to directory of reference data [default: %default]"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project [default: Nameset of input]"
  )
  # make_option(
    # c("-c", "--outcome"),
    # type = "character",
    # help = "String specifying the column in data containing the outcome variable [required]"     
  # )
)

args_list <- parse_args(OptionParser(option_list = option_list))

# Nameset default case
if(is.null(args_list$nameset)){
  args_list$nameset <- nameset_default(args_list$input)
}

# read in fields
fields <- c(
  "input",
  "output",
  "ref",
  "nameset"
  # "outcome"
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

# Functions -------------------------


# Main -------------------------
# DE file locations
fh <- file.path("../..", "baseline_PCs_relatedness", input,
  paste("dataProt_SS-205063.PC.related", input, "csv", sep = ".")
)
result <- lapply(fh, fread)
setattr(result, 'names', input)

# merge
result <- rbindlist(result)
# bottom <- round(min(result[["t stat"]], na.rm = TRUE) - 0.05, digits = 2)
# top <- round(max(result[["t stat"]], na.rm = TRUE) + 0.05, digits = 2)
# # drop APOE from input so it will match the results
input <- gsub(".APOE", "", input)

# thresholds of cutoffs, pval, padj, outer_decile, outer_quartile
result[, stat_pval := ifelse(pval < 0.05, `t stat`, NA)]
result[, stat_padj := ifelse(padj < 0.05, `t stat`, NA)]
lapply(input, function(i) {
  qs <- quantile(result[dx == i, `t stat`], c(0.1, 0.25, 0.75, 0.9))
  result[, stat_decile := ifelse(`t stat` < qs[[1]] | `t stat` > qs[[4]], `t stat`, NA)]
  result[, stat_quartile := ifelse(`t stat` < qs[[2]] | `t stat` > qs[[3]], `t stat`, NA)]
})

# check quantiles to see distributions
# lapply(input, function(i) {
  # qs <- quantile(result[dx == i, `t stat`], c(0.1, 0.25, 0.75, 0.9))
  # print(i)
  # print(qs)
# })

# read in biodomains
types = c("cultured", "isolated")
trt = lapply(types, function(i) {
  fread(paste0(ref, "/Sharma_etal_prot_celltypes.", i, ".scaled.csv"))
})

# merge tables
dt = lapply(trt, function(i) {
  merge(i, 
    result[, .(GENE, dx, `t stat`, stat_pval, stat_padj, stat_decile, stat_quartile)], 
    by = "GENE",
    allow.cartesian = TRUE
  )
})

# calc correlations (exclude all is.na(beta) )
plots = c("t stat", "stat_pval", "stat_padj", "stat_decile", "stat_quartile")
lapply(plots, function(k) {
  # checking the number of genes matched (not NA in both) per group to any mouse de
  counts <- lapply(dt, function(j) {
    j[celltype == "Neurons" & !is.na(get(k)) & !is.na(z.log2LFQ), .N, by = dx
      ][match(input, dx)]
  })
  
  corr_list = lapply(dt, function(j) lapply(input, function(i){
    j[dx == i & !is.na(get(k)) & !is.na(z.log2LFQ)
      ][, N := .N, by=celltype
        ][N >= 2
          ][, cor.test(get(k), z.log2LFQ, method="pearson")[-c(2, 9)], by=celltype]
  }))

  # label and collapse
  lapply(corr_list, function(i) {
    setattr(i, 'names', input)
  })
  corrs = lapply(corr_list, function(i) rbindlist(i, idcol="group"))

  # calculate FDR per group (normally across all groups if not significance thresholds)
  lapply(corrs, function(i) {
    i[, FDR := p.adjust(p.value, method = "fdr"), by = group]
  })

  # get matricies for heatmaps
  lapply(corrs, function(i) i[, group := factor(group, levels = input)])
  a = lapply(corrs, function(i) {
    as.matrix(
      x = dcast(i, group ~ celltype, value.var = "estimate"), 
      rownames = "group"
    )
  })
  b = lapply(corrs, function(i) {
    as.matrix(
      x = dcast(i, group ~ celltype, value.var = "p.value"), 
      rownames = "group"
    )
  })

  # plot cors
  p.list <- lapply(seq(corrs), function(i) {
    ggcorrplot(a[[i]], 
      method = "circle",
      ggtheme = ggplot2::theme_classic,
      p.mat = b[[i]], 
      sig.level = fdr_cutoff, 
      colors = c("#900D25", "#FEFEFE", "#15508E"),
      outline.color = "#f0f0f000",
      lab = TRUE,
      lab_col = '#4D4D4D', 
      tl.col = 'black', 
      tl.srt = 60
    ) + 
      coord_flip() + 
      scale_x_discrete(limits = rev(rownames(a[[i]]))) +
      scale_y_discrete(position = "right") +
      ggtitle(paste(str_to_title(types[i]), "cells")) +
      # secondary axis with gene counts
      guides( 
        y.sec = guide_axis_manual(
          breaks = rev(counts[[i]][, .I]), 
          labels = counts[[i]][, N],
          label_size = 10,
          angle = 0,
          position = 'right',
          title = "Number of genes"
        )
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 60, vjust = 0, hjust = 0)
      )
  })
  ggsave(wrap_plots(p.list) + plot_layout(guides = 'collect', axes = 'collect_y'),
    filename = file.path(output, 
      paste(nameset, k, "Sharma_etal_correlations.png", sep = ".")
    ), 
    width = 11, 
    height = 4 + 0.25,
    units = "in", 
    dpi = 144, 
    type = "cairo"
  )
  ggsave(wrap_plots(p.list),
    filename = file.path(output, 
      paste(nameset, k, "Sharma_etal_correlations.pdf", sep = ".")
    ), 
    width = 11, 
    height = 4,
    units = "in", 
    device = cairo_pdf
  )
})















