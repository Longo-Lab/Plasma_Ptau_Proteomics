#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 04-TREAT-AD_correlation.R
##
## Version: 0.0.3
##
## Purpose of script: Show Correlation of TREAT-AD pathway traced 
## modules with Ptau geneset
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
library(tools)
library(ggcorrplot)
library(ggh4x)
library(gtools)

# Global variables
# ensembl version
ens <- "113"
# ortholog conversion lookup
orth <- "one2one"

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
    default = "../biodomains",
    help = "File path to directory of reference data [default: %default]"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project [default: Nameset of input]"
  ),
  make_option(
    c("-c", "--outcome"),
    type = "character",
    help = "String specifying the column in data containing the outcome variable [required]"     
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
  "nameset",
  "outcome"
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

# Functions -------------------------


# Main -------------------------
## DE file locations
# load files (dataProt_SS-205063.pTau.baseline.csv)
# `bm`~ biomarker
bm = fread(input)

#thresholds of cutoff
bm[, stat_pval := ifelse(pval < 0.05, `t stat`, NA)]
bm[, stat_padj := ifelse(padj < 0.05, `t stat`, NA)]

# read in biodomains
trt = fread(file.path(ref,
  paste(modules, ens, orth, "Human.Genes.txt", sep = ".")
))

# table downloaded from synapse
# synapse get -q "SELECT * FROM syn22758536"
expr = fread(file.path(ref,
  paste(modules, ens, orth, "Human.Corr.txt", sep = ".")
))
# expr[, Pro_stat := -log10(Pro_fdr_CorPVal) * Pro_TE]

# merge tables
dt = merge(trt, expr[, .(GENE, Pro_stat)], all.x = TRUE, by = "GENE")
dt = merge(
  bm[, .(GENE, `t stat`, stat_pval, stat_padj)], 
  dt, 
  all.x = TRUE, 
  by = "GENE"
)
dt <- dt[!is.na(module)]

# checking the number of genes matched (not NA in both) per group to any mouse de
counts <- dt[
  !is.na(Pro_stat) & (!is.na(`t stat`) | !is.na(stat_pval) | !is.na(stat_padj)),
  .N, 
  by = module
]
counts$module = factor(counts$module, levels = mixedsort(unique(counts$module)))
setorder(counts, module)


# calc correlations (exclude all is.na(beta) ) ----------------------------
cors = c("t stat", "stat_pval", "stat_padj")
corr_list = lapply(cors, function(i){
  dt[!is.na(get(i)) & !is.na(Pro_stat)
     ][, N := .N, by=module
       ][N >= 2
         ][, cor.test(get(i), Pro_stat, method="spearman")[-2], by=module]
})
setattr(
  corr_list, 
  'names', 
  c(
    outcome, 
    sprintf("%s Pval < 0.05", outcome), 
    sprintf("%s Padj < 0.05", outcome)
  )
)
corrs = rbindlist(corr_list, idcol = "group")
corrs[, FDR := p.adjust(p.value, method = "fdr"), by = group]
corrs$module = factor(corrs$module, levels = mixedsort(unique(counts$module)))
corrs$group = factor(corrs$group, levels = unique(corrs$group))

# plot cors
plot_corr = function(corrs_data, outname, width=11, height=4, fdr_cutoff=0.05) {
  a <- as.matrix(
    dcast(corrs_data, group ~ module, value.var="estimate"), 
    rownames = "group"
  )
  b <- as.matrix(
    dcast(corrs_data, group ~ module, value.var="FDR"), 
    rownames = "group"
  )
  
  # legend title
  test_name <- unique(corrs_data$method)[1]
  if (test_name %like% 'Ken') {
    titl <- "Kendall's\n\U03C4"
  } else if (test_name %like% 'Spe') {
    titl <- "Spearman's\n\U03C1"
  } else if (test_name %like% 'Pea') {
    titl <- "Pearson's\n\U03C1"
  }
  
  p <- ggcorrplot(a, 
    method = "circle",
    ggtheme = ggplot2::theme_classic,
    p.mat = b, 
    sig.level = fdr_cutoff, 
    # insig = "blank", 
    # colors = brewer.pal(9, "RdBu")[ c(1, 5, 9)],
    colors = c("#900D25", "#FEFEFE", "#15508E"),
    outline.color = "#f0f0f000",
    legend.title = titl,
    lab = TRUE,
    lab_col = '#4D4D4D', 
    tl.col = 'black', 
    tl.srt = 60
  ) + 
    coord_flip() + 
    scale_x_discrete(limits = rev(rownames(a))) +
    scale_y_discrete(position = "right") +
    # secondary axis with gene counts
    guides( 
      x.sec = guide_axis_manual(
        breaks = counts[, .I], 
        labels = counts[, N],
        label_size = 10,
        angle = 0,
        position = 'bottom',
        title = "Number of genes"
      )
    ) +
    theme(
      axis.text.x = element_text(angle = 60, vjust = 0, hjust = 0),
      legend.title = element_text(face = 'bold', hjust = 0.5)
    )
  ggsave(p,
    filename = file.path(output, paste(outname, "png", sep=".")), 
    width = width, 
    height = height + 0.25,
    units = "in", 
    dpi = 144, 
    type = "cairo"
  )
  ggsave(p,
    filename = file.path(output, paste(outname, "pdf", sep=".")), 
    width = width, 
    height = height,
    units = "in", 
    device = cairo_pdf
  )

}
plot_corr(
  corrs,
  outname = paste(nameset, modules, "correlations", sep = "."),
  height = ifelse(modules == "TMTSoma-CSF", 5, 4),
  width = ifelse(modules == "TMTSoma-CSF", 20, 11)
)

# # plot cors
# a = as.matrix(dcast(corrs, group ~ module, value.var="estimate"), 
              # rownames="group")
# b = as.matrix(dcast(corrs, group ~ module, value.var="FDR"), 
              # rownames="group")
# pdf(
  # file = file.path(output, 
    # paste(nameset, modules, "correlations.pdf", sep = ".")
  # ), 
  # width = 11, 
  # height = 8
# )
# corrplot(a, 
  # p.mat = b, 
  # sig.level = 0.05, 
  # addCoef.col = 'dark grey', 
  # na.label = "-"
# )
# dev.off()
