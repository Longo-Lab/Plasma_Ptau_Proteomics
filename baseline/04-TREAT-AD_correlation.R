#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 04-TREAT-AD_correlation.R
##
## Version: 0.0.2
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
library(corrplot)
library(stringr)
library(optparse)
library(tools)

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
trt = fread(paste0(ref, "/TREAT-AD.Human.Genes.txt"))

# table downloaded from synapse
# synapse get -q "SELECT * FROM syn22758536"
expr = fread(paste0(ref, "/SYNAPSE_TABLE_QUERY_96970654.csv"))
expr[, Pro_stat := -log10(Pro_fdr_CorPVal) * Pro_TE]

# merge tables
dt = merge(trt, expr[, .(ENSG, Pro_stat)], all.x = T, by.x = "GENE", by.y = "ENSG")
dt = merge(dt, bm[, .(GENE, `t stat`, stat_pval, stat_padj)], all.x=T, by="GENE")

# calc correlations (exclude all is.na(beta) )
cors = c("t stat", "stat_pval", "stat_padj")
corr_list = lapply(cors, function(i){
  dt[!is.na(get(i)) & !is.na(Pro_stat)
     ][, N := .N, by=module
       ][N >= 2
         ][, cor.test(get(i), Pro_stat, method="spearman")[-2], by=module]
})
setattr(corr_list, 'names', c(outcome, sprintf("%s Pval < 0.05", outcome), sprintf("%s Padj < 0.05", outcome)))
corrs = rbindlist(corr_list, idcol="group")
corrs[, FDR := p.adjust(p.value, method="fdr"), by=group]

# plot cors
a = as.matrix(dcast(corrs, group ~ module, value.var="estimate"), 
              rownames="group")
b = as.matrix(dcast(corrs, group ~ module, value.var="FDR"), 
              rownames="group")
pdf(file=file.path(output, paste(nameset, "TREAT-AD_correlations.pdf", sep=".")), width=11, height=8)
corrplot(a, p.mat=b, sig.level=0.05, insig="blank", addCoef.col='dark grey', 
         na.label="-")
dev.off()
