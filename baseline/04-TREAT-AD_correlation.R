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
## Author: Robert B Butler III
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

# directory locations
rootdir = '/mnt/box/plasma_ptau_proteomics'
wk_dir = paste0(rootdir, '/baseline')
setwd(wk_dir)
refdir = paste0(rootdir, '/biodomains')
prefix = "dataProt_SS-205063.pTau.baseline"
nameset = "TREAT-AD"

# Functions -------------------------


# Main -------------------------
## DE file locations
# load files
ptau = fread(paste(prefix, "csv", sep="."))

#thresholds of cutoff
ptau[, stat_pval := ifelse(pval < 0.05, `t stat`, NA)]
ptau[, stat_padj := ifelse(padj < 0.05, `t stat`, NA)]

# read in biodomains
trt = fread(paste0(refdir, "/TREAT-AD.Human.Genes.txt"))

# table downloaded from synapse
# synapse get -q "SELECT * FROM syn22758536"
expr = fread(paste0(refdir, "/SYNAPSE_TABLE_QUERY_96970654.csv"))

# merge tables
dt = merge(trt, expr[, .(ENSG, Pro_TE)], all.x=T, by.x="GENE", by.y="ENSG")
dt = merge(dt, ptau[, .(GENE, `t stat`, stat_pval, stat_padj)], all.x=T, by="GENE")

# calc correlations (exclude all is.na(beta) )
cors = c("t stat", "stat_pval", "stat_padj")
corr_list = lapply(cors, function(i){
  dt[!is.na(get(i)) & !is.na(Pro_TE)
     ][, N := .N, by=module
       ][N >= 2
         ][, cor.test(get(i), Pro_TE, method="spearman")[-2], by=module]
})
setattr(corr_list, 'names', c("PlasmaPTau181", "PlasmaPTau181 Pval < 0.05", "PlasmaPTau181 Padj < 0.05"))
corrs = rbindlist(corr_list, idcol="group")
corrs[, FDR := p.adjust(p.value, method="fdr"), by=group]

# plot cors
a = as.matrix(dcast(corrs, group ~ module, value.var="estimate"), 
              rownames="group")
b = as.matrix(dcast(corrs, group ~ module, value.var="FDR"), 
              rownames="group")
pdf(file=paste(prefix, nameset, "correlations.pdf", sep="."), width=11, height=8)
corrplot(a, p.mat=b, sig.level=0.05, insig="blank", addCoef.col='dark grey', 
         na.label="-")
dev.off()
