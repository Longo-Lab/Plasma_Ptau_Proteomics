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
library(corrplot)
library(stringr)
library(optparse)

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
types = c("cultured", "isolated")
trt = lapply(types, function(i) {
  fread(paste0(ref, "/Sharma_etal_prot_celltypes.", i, ".scaled.csv"))
})

# merge tables
dt = lapply(trt, function(i) {
  merge(i, bm[, .(GENE, `t stat`, stat_pval, stat_padj)], all.x=T, by="GENE")
})

# calc correlations (exclude all is.na(beta) )
cors = c("t stat", "stat_pval", "stat_padj")
corr_list = lapply(dt, function(j) lapply(cors, function(i){
  j[!is.na(get(i)) & !is.na(z.log2LFQ)
    ][, N := .N, by=celltype
      ][N >= 2
        ][, cor.test(get(i), z.log2LFQ, method="pearson")[-c(2, 9)], by=celltype]
}))

# label and collapse
lapply(corr_list, function(i) {
  setattr(i, 'names', c(outcome, 
                        sprintf("%s Pval < 0.05", outcome), 
                        sprintf("%s Padj < 0.05", outcome)))
})
corrs = lapply(corr_list, function(i) rbindlist(i, idcol="group"))

# calculate FDR per group (normally across all groups if not significance thresholds)
lapply(corrs, function(i) i[, FDR := p.adjust(p.value, method="fdr"), by=group])

# get matricies for heatmaps
lapply(corrs, function(i) i[, group := factor(group, levels=names(corr_list[[1]]))])
a = lapply(corrs, function(i) {
  as.matrix(dcast(i, group ~ celltype, value.var="estimate"), rownames="group")
})
b = lapply(corrs, function(i) {
  as.matrix(dcast(i, group ~ celltype, value.var="p.value"), rownames="group")
})

# plot cors
pdf(file=file.path(output, paste(nameset, "Sharma_etal_correlations.pdf", sep=".")), width=14, height=6)
par(mfrow = c(2, 2))
lapply(seq(a), function(i) {
  corrplot(a[[i]], p.mat=b[[i]], sig.level=0.05, method="ellipse",
           addCoef.col='dark grey', na.label="-", mar=c(0,0,2,0), tl.cex=1.5,
           tl.col="black", tl.srt=45, cl.ratio=0.3,
           title=sprintf("%s celltypes", types[i]))$corrPos -> p1
  text(p1$x, p1$y, round(p1$corr, 2), col='dark grey')
})
dev.off()















