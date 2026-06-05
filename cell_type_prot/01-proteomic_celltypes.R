#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 01-proteomic_celltypes.R
##
## Version: 0.0.1
##
## Purpose of script: Generate proteomic celltype correlation references for 
## Sharma et al 2015
##
## Author: Robert B Butler III
##
## Date Created: 2023-09-14
##
## Copyright (c) 2023
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##
##   Operates on the Log2 LFQ intensity values for each xlsx, copied into csv.
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(data.table)
library(biomaRt)

# directory locations
rootdir = '/mnt/c/Users/rrbutler/Documents/SRCF_shutdown'
wk_dir = paste0(rootdir, '/plasma_ptau_proteomics/cell_type_prot')
setwd(wk_dir)
refdir = paste0(rootdir, '/reference/biodomains')

# variables
nameset = "Sharma_etal_prot_celltypes"
types = c("cultured", "isolated")
ensver = "105"

# Functions -----------------------

#' For a given column of Gene Symbols, replace them with Ensembl IDs
#' includes a filter for only protein coding genes
#'
#' @param dt data table to replace names
#' @param colname name of the column containing ensgs
#' @param keep.ensgs boolean, retain the ENSG column?
#'
#' @returns dt with a GENE column with or without the Symbol column
convertGeneSymbol = function(dt, colname, keep.symbols=F){
  require("biomaRt")
  print(paste(length(unique(dt[[colname]])), "genes to look up"))
  ensembl = useEnsembl("ensembl", dataset="mmusculus_gene_ensembl", version=ensver)
  # match with ensembl_gene_id across both
  x = data.table(getBM(attributes=c("ensembl_gene_id", "external_gene_name",
                                    "gene_biotype", "chromosome_name"),
                       filters="external_gene_name", 
                       values=unique(dt[[colname]]), 
                       mart=ensembl))
  print(paste(nrow(x), "total genes returned"))
  # restricting to only protein coding 
  keep_categ = c('IG_C_gene',
                 'IG_D_gene',
                 'IG_J_gene',
                 'IG_V_gene',
                 'protein_coding',
                 'TR_C_gene',
                 'TR_D_gene',
                 'TR_J_gene',
                 'TR_V_gene')
  x = x[gene_biotype %in% keep_categ]
  print(paste(nrow(x), "genes correct biotype"))
  # use only canonical chromosomes
  use_chrs = c(1:19, "X", "Y")
  x = x[chromosome_name %in% use_chrs]
  print(paste(nrow(x), "genes in chr 1:19, X, Y"))
  setnames(x, c("ensembl_gene_id", "external_gene_name"), c("GENE", "Symbol"))
  # replace gene names and write to file
  dt = merge(x[, .(GENE, Symbol)], dt, all.y=T, by.x="Symbol", by.y=colname)
  dt = dt[!is.na(GENE)]
  print(paste(length(unique(dt[["GENE"]])), "unique genes with ensembl ids"))
  if (keep.symbols == FALSE) dt[, Symbol := NULL]
  return(unique(dt))
}

#' For a given column of human ensIDs, replace them with mouse ensembl gene ids 
#' includes a filter for only protein coding genes
#'
#' @param dt data table to replace names
#' @param colname name of the column containing mouse ensgs
#' @param keep.mouse boolean, retain the mouse ENSG column?
#'
#' @returns dt with a GENE column with or without the mGENE column
convertMouseHuman = function(dt, colname, keep.mouse=F){
  require("biomaRt")
  print(paste(length(unique(dt[[colname]])), "genes to look up"))
  mouse = useEnsembl("ensembl", dataset="mmusculus_gene_ensembl", version=ensver)
  human = useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", version=ensver)
  # match with ensembl_gene_id across both
  x = data.table(getLDS(mart=mouse, attributes=c('ensembl_gene_id'), martL=human,
                        attributesL=c('ensembl_gene_id', 'chromosome_name', 
                                      'gene_biotype'), 
                        filters='ensembl_gene_id', values=unique(dt[[colname]]),
                        bmHeader=F))
  names(x) = make.unique(names(x))
  print(paste(nrow(x), "total genes returned"))
  # restricting to only protein coding 
  keep_categ = c('IG_C_gene',
                 'IG_D_gene',
                 'IG_J_gene',
                 'IG_V_gene',
                 'protein_coding',
                 'TR_C_gene',
                 'TR_D_gene',
                 'TR_J_gene',
                 'TR_V_gene')
  x = x[gene_biotype %in% keep_categ]
  print(paste(nrow(x), "genes correct biotype"))
  # use only canonical chromosomes
  use_chrs = c(1:22, "X", "Y", "MT")
  x = x[chromosome_name %in% use_chrs]
  print(paste(nrow(x), "genes in chr 1:22 X Y MT"))
  setnames(x, c("ensembl_gene_id", "ensembl_gene_id.1"), c("mGENE", "GENE"))
  # replace gene names and write to file
  dt = merge(x[, .(GENE, mGENE)], dt, all.y=T, by.x="mGENE", by.y=colname)
  dt = dt[!is.na(GENE)]
  print(paste(length(unique(dt[["GENE"]])), "unique genes with ensembl ids"))
  if (keep.mouse == FALSE) dt[, mGENE := NULL]
  return(unique(dt))
}

#' Applies scale to the rows of a data.table with double transposition
#' (t -> scale -> t)
#'
#' @param dt data table to replace names
#' @param key_col name of the key column, cannot be a combination of columns
#' @param scale_cols vector of column names included in the scaling
#'
#' @returns dt with scale_cols scaled
row_scale = function(dt, key_col, scale_cols){
  cols = c(key_col, scale_cols)
  x = transpose(dt[, ..cols], keep.names="rn", make.names=key_col)
  scale_rows = names(x)[-1]
  x[, (scale_rows) := lapply(.SD, scale), .SDcols=scale_rows]
  y = transpose(x, keep.names=key_col, make.names="rn")
  others = setdiff(names(dt), scale_cols)
  z = merge(dt[, ..others], y, by=key_col, all.x=T, sort=F)
  setcolorder(z, names(dt))
  return(z)
}

# Main -----------------------
# read in list of files 
fl = lapply(types, function(i) fread(paste(nameset, i, 'csv', sep=".")))
setattr(fl, 'names', types)
lapply(fl, function(i) setnames(i, "Gene names", "GENE"))

# explode prot names column by ";"
fl = lapply(fl, function(i){
  i[, .(GENE=tstrsplit(GENE, ";", fixed=TRUE)), by=c(names(i)[-1])]
  i[, GENE := unlist(GENE)]
})

# Convert to ENSMUSG, then human ENSG
fl = lapply(fl, function(i) convertGeneSymbol(i, "GENE")) 
fl = lapply(fl, function(i) convertMouseHuman(i, "GENE")) 

# there are some duplicated ENSGs, randomly select one
set.seed(123)
fl = lapply(fl, function(i) i[, .SD[sample(.N)[1]], by=GENE])

# row scale the genes 
fl = lapply(fl, function(i) row_scale(i, "GENE", setdiff(names(i), "GENE")))

# melt the table into celltype & z.log2LFQ
fm = lapply(fl, function(i) {
  melt(i, variable.name="celltype", value.name="z.log2LFQ", id.vars="GENE")
})

# save to file
lapply(seq(fm), function(i) {
  fwrite(fm[[i]], file=paste(nameset, types[i], "scaled.csv", sep="."))
})

# # Specificity Generation --------------------------------------------------
# fh <- fread("Sharma_etal.EWCE.human.lvl1.spec_matr.txt")
# fh <- 




