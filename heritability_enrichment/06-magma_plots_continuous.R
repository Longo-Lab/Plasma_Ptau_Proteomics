#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 09-magma_plots_continuous.R
##
## Version: 0.0.1
##
## Purpose of script: Generate heatmaps for continuous heritability enrichment analysis 
##
## Author: Robert B Butler III
##
## Date Created: 2024-01-17
##
## Copyright (c) 2024
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##   Need to fix ordering...
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(ComplexHeatmap)
library(circlize)
library(data.table)

# directory locations
# file naming info
today = format(Sys.Date(), '%Y%m%d')
nameset = 'continuous_property_magma'
f_names = c('betaLOGP', 'signLOGP', 't_stat')
f_prefix = 'AdultBrain'
gwas = c('ALZ', 'ALS', 'FTD', 'LBD', 'PD', 'rAOO', 'BMI')



# Functions -------------------------------------------------------------------
heat = function(pval, fdr, tone, name) {
  col_fun = colorRamp2(c(0, 5), c("#FFFFFF", tone))
  a = Heatmap(pval, 
    cluster_rows = FALSE, 
    cluster_columns = FALSE, 
    col = col_fun, 
    name = "-Log P", 
    column_title = name, 
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    column_names_rot = 60, 
    column_names_gp = gpar(fontface = "bold"),
    # column_split = c(ceiling(ncol(pval) / 2)),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.rect(
        x = x, 
        y = y, 
        width = width, 
        height = height, 
        gp = gpar(
          lwd = 2, 
          col = ifelse(
            fdr[i, j] < 0.05, 
            "black", 
            ifelse(pval[i, j] > -log10(0.05), "#969696", NA)
          )
        )
      )
    }
  )
  return(a)
}

# MAIN ------------------------------------------------------------------------
# reading in tables, combine and melt by GWAS
fn = as.vector(outer(gwas, f_names, paste, sep="."))
fl = lapply(paste0(sub("\\.", paste0("_", f_prefix, "_"), fn), ".gsa.out"),
            fread, skip="VARIABLE")
setattr(fl, 'names', fn)
fd = rbindlist(fl, use.names=T, idcol="X")

# generate columns for metadata, logp, fdr
fd[, c("GWAS", "DE") := tstrsplit(X, "\\.")]
# fd[, c("GENESET", "DIR") := tstrsplit(VARIABLE, "\\.")]
# fd[, DIR_NA := is.na(DIR)]
fd[, LOGP := -log10(P)]
fd[, FDR := p.adjust(P, method = "fdr"), by=c("DE")]
fd$GWAS <- factor(fd$GWAS, levels = gwas)
levels(fd$GWAS) <- gsub("_g1000", "", levels(fd$GWAS)) # drop hm3 label from plots


# all DE genes ----------------------------------------------------------------
# pvalue matrices
mats = lapply(f_names, function(x){
  a = as.matrix(dcast(fd[DE == x], GWAS ~ VARIABLE,
                value.var = "LOGP"), rownames="GWAS")
})

# fdr matrices
fdr_mats = lapply(f_names, function(x){
  a = as.matrix(dcast(fd[DE == x], GWAS ~ VARIABLE, 
                value.var = "FDR"), rownames="GWAS")
})

# Heatmap
ht_list = heat(mats[[1]], fdr_mats[[1]], "#1b7837", paste("DE continuous", f_names[1])) +
  heat(mats[[2]], fdr_mats[[2]], "#762a83", paste("DE continuous", f_names[2])) +
  heat(mats[[3]], fdr_mats[[3]], "#8c510a", paste("DE continuous", f_names[3]))
pdf(
  file = paste(today, nameset, "DE", "heatmap.pdf", sep="."), 
  width = 11, 
  height = 7,
)
draw(ht_list, ht_gap=unit(2, "cm"))
dev.off()

# ht_list = heat(mats[[4]], fdr_mats[[4]], "#8c510a", paste("Continuous", f_names[4]))
# pdf(file = paste(today, nameset, "specificity", "heatmap.pdf", sep="."))
# draw(ht_list)
# dev.off()


# save table ------------------------------------------------------------------
# fd = rbindlist(fl, use.names=T, idcol="ANALYSIS")
fd[, c("X", "TYPE") := NULL]
fwrite(fd, file=paste(today, nameset, "heatmap_stats.txt", sep="."), sep="\t", quote=F)


