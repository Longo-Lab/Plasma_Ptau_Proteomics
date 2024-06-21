#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 06-Milind_etal.R
##
## Version: 0.0.2
##
## Purpose of script: Show Fisher's Exact Enrichment of Submodules of 
## human/mouse co-expression with Ptau geneset
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
##   paper DOI: 10.1371/journal.pgen.1008775
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(ggplot2)
library(ggpubr)
library(data.table)
library(stringr)
library(plyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(GeneOverlap)
library(biomaRt)

# directory locations
rootdir = '/mnt/c/Users/rrbutler/Documents/SRCF_shutdown'
wk_dir = paste0(rootdir, '/plasma_ptau_proteomics/baseline')
setwd(wk_dir)
refdir = paste0(rootdir, '/reference/biodomains')
prefix = "dataProt_SS-205063.pTau.baseline"
nameset = "Milind_etal"

# Functions -------------------------------------
gommat_2_dt = function(go, dt, dtype=c("pval", "odds.ratio", 
                       "intersection", "union", "Jaccard")) {
  # add gom object to a data.table
  # gom.obj: GeneOverlap object
  # dt: data.table to add to
  # dtype: which GO measure to extract
  d = melt(data.table(getMatrix(go, name=dtype), keep.rownames="sample"), 
           variable.name="group", value.name=dtype, id.vars="sample")
  dt = merge(dt, d, by=c("sample", "group"))
  return(dt)
}

get_mcols = function(str, result, col) {
  # for a string of ens ids: explode, look up symbols, return symbols as string
  # string: comma separated ens_ids
  # result: data table with ens_ids in Gene_id
  # col: colname to grab
  require("stringr")
  y = result[[col]][match(unlist(str_split(str, ",")), result$GENE)]
  return(paste0(y, collapse=","))
}

gomtable = function(list1, list2, background, outname, width=7, 
                    height=7, ...) {
  # Takes two gene overlap lists and writes a table of stats, makes list object
  # list1: named list object
  # list2: named list object
  # outname: name prefix for file
  # width: width of pdf files to be plotted in inches*2
  # height: height of pdf files
  # returns a list object of table, r_names, c_names
  # enrichment basic
  gom.obj = newGOM(list1, list2, background)
  r_names = names(list1)
  c_names = names(list2)
  # # geneoverlap heatmap
  # pdf(file=paste(prefix, nameset, outname, "pval.pdf", sep="."), 
      # width=width*2, height=height)
  # drawHeatmap(gom.obj, grid.col="Blues", cutoff=1)
  # dev.off()
  # write DE vals to file
  b = melt(data.table(getMatrix(gom.obj, name="pval"), keep.rownames="sample"), 
           variable.name="group", value.name="P", id.vars="sample")
  # write dtypes to table
  b = gommat_2_dt(gom.obj, b, dtype="odds.ratio")
  b = gommat_2_dt(gom.obj, b, dtype="Jaccard")
  b = gommat_2_dt(gom.obj, b, dtype="intersection")
  b = gommat_2_dt(gom.obj, b, dtype="union")
  # add fdr
  b[, FDR := p.adjust(P, method = "fdr") ]
  b[, LOGP := -log10(P) ]
  b[, LOGFDR := -log10(FDR) ]
  # get intersection genes
  nest = getNestedList(gom.obj, name="intersection")
  nest2 = lapply(seq(nest), function(i) {
    lapply(nest[[i]], function(j) paste0(j, collapse=","))
  })
  nest3 = lapply(seq(nest2), function(i) stack(nest2[[i]]))
  setattr(nest3, 'names', names(nest))
  dt = rbindlist(nest3, idcol=T)
  setnames(dt, c("group", "ens_id", "sample"))
  dt$symbol = lapply(dt$ens_id, get_mcols, result, "EntrezGeneSymbol")
  b = merge(b, dt, by=c("sample", "group"), all.x=T, sort=F)
  # write and return
  fwrite(b, file=paste(prefix, nameset, outname, "txt", sep="."), sep="\t", 
         quote=F)
  return(list("b"=b, "r_names"=r_names, "c_names"=c_names))
}

enrich = function(list1, list2, gs, outname, width=7, height=7, fdr=TRUE, ...) {
  # function runs enrichment and saves to file for odds ratio and padj/p
  # takes two set lists, background gene count and output name and pdf width
  # DE enrich
  de_enrich = gomtable(list1, list2, gs, outname, width=width, height=height)

  # fetch Jaccard index
  f = as.matrix(dcast(de_enrich$b, sample ~ group, value.var="Jaccard"), 
                rownames="sample")
  f = f[de_enrich$r_names, de_enrich$c_names, drop=F]
  # log10p heatmap
  if (isTRUE(fdr)) {
    e = as.matrix(dcast(de_enrich$b, sample ~ group, value.var="LOGFDR"), 
                  rownames="sample")
    e = e[de_enrich$r_names, de_enrich$c_names, drop=F]
    # filename and heatmap
    pdf(file=paste(prefix, nameset, outname, "logfdr.pdf", sep="."), width=width, 
        height=height)
    heat(e, e, f, "-log FDR")
    dev.off()
  } else { 
    e = as.matrix(dcast(de_enrich$b, sample ~ group, value.var="LOGP"), 
                  rownames="sample")
    e = e[de_enrich$r_names, de_enrich$c_names, drop=F]
    # filename and heatmap
    pdf(file=paste(prefix, nameset, outname, "logp.pdf", sep="."), width=width, 
        height=height)
    heat(e, e, f, "-log P")
    dev.off()
  }
  # OR heatmap
  c = as.matrix(dcast(de_enrich$b, sample ~ group, value.var="odds.ratio"), 
                rownames="sample")
  c = c[de_enrich$r_names, de_enrich$c_names, drop=F]
  # filename and heatmap
  pdf(file=paste(prefix, nameset, outname, "logOR.pdf", sep="."), width=width, 
      height=height)
  heat(c, e, f, "Odds Ratio")
  dev.off()
}

heat = function(mat, pval, int, name) {
  # plot heatmap of overlaps colored by value, and sized by overlap (given matrices, name)
  # mat: matrix of heatmap values (fdr, p-val, or odds.ratio)
  # pval: matched matrix of fdr values, to bold significant circles
  # int: matched matrix of Jaccard indexes for intersection sizes
  # name: heatscale name for plot (fdr, p-val, or odds.ratio)
  # returns a draw graph call

  # resolve Inf to real number (round) and check for NAs
  mat[is.infinite(mat)] = -log10(2.225074e-308)
  pval[is.infinite(pval)] = -log10(2.225074e-308)
  if (any(is.na(int))) {
    print("Jaccard index contains NAs!")
    quit("no", 1)
  }

  # color palette
  col_fun = colorRamp2(seq(min(mat, na.rm=T), max(max(mat, na.rm=T), 3), length=9), 
                       brewer.pal(name="BuPu", n=9))
  col_fun2 = brewer.pal(name="Dark2", n=3)
  # max value in each column of matrix to order by 
  col_ord = data.table(max=sapply(colnames(mat), function(i) max(mat[,i])), module=colnames(mat))
  # study each column was from for annotation
  col_ord[['study']] = mapvalues(col_ord$module, from=tmt_anno$module, to=tmt_anno$Study)
  # pathway enrichment annotation for each module/column
  col_ord[['path']] = mapvalues(col_ord$module, from=tmt_anno$module, to=tmt_anno$Pathway_Enrichment)
  # reorder matrices for plotting by size within group
  col_ord = col_ord[order(-study, -max, module)]
  names(col_fun2) = unique(col_ord$study)
  mat = mat[, col_ord$module]
  pval = pval[, col_ord$module]
  int = int[, col_ord$module]
  ## Column annotations mapped from annotation table for Milind below
  ha = HeatmapAnnotation("Pathway"=anno_text(col_ord$path, rot=60, just="left", 
                                             location=unit(0, 'npc')), 
                         "Study"=col_ord$study,
                         col=list("Pathway"="black", "Study"=col_fun2))
  # scale max int (for circle size) to matrix "buffer"
  buffer = ((int/max(int))/4 + .25)
  dot_scalar = unit(1.5, "cm")
  # heatmap plot 
  a = Heatmap(mat, cluster_rows=F, cluster_columns=F, na_col="#FFFFFF", 
              col=col_fun, name=name, column_names_rot=60, row_names_side="left",
              top_annotation=ha, 
              column_split=factor(col_ord$study, levels=c("ROSMAP", "Mayo", "MSBB")), 
              column_gap=unit(5, "mm"), 
              cluster_column_slices=F,
              column_title_gp = gpar(fontsize = 20, fontface = "bold"),
              rect_gp=gpar(type="none"), 
              cell_fun=function(j, i, x, y, width, height, fill) { 
                grid.circle(x=x, y=y, r=(buffer[i, j] * dot_scalar), 
                            gp=gpar(fill=col_fun(mat[i, j]), 
                                    col=ifelse(pval[i, j] < -log10(0.05), NA, 
                                               "black")))
                if (pval[i, j] > -log10(0.05)) {
                  grid.text(sprintf("%.2f", int[i, j]), x, y, 
                            gp=gpar(fontsize=10, fontface="bold", 
                                    col=ifelse(mat[i, j] > quantile(mat, 0.95), 
                                               "white", "black")))
                }
              })
  # heatmap legend
  lgd_list = list(
    Legend(labels=c(rep("", 3)), 
           title = "Intersection Size",
           grid_height = unit(12, "mm"), 
           grid_width = unit(12, "mm"), 
           graphics = list(
              function(x, y, w, h) {
                grid.circle(x, y, r=min(buffer) * dot_scalar,
                            gp=gpar(fill="black", col="white"))
                grid.text(sprintf("%.2f", 0), x, y, 
                          gp=gpar(fontsize=10, fontface="bold", col="white"))
              },
              function(x, y, w, h) {
                grid.circle(x, y, r=mean(buffer) * dot_scalar,
                            gp=gpar(fill="black", col="white"))
                grid.text(sprintf("%.2f", mean(int)), x, y, 
                          gp=gpar(fontsize=10, fontface="bold", col="white"))
              },
              function(x, y, w, h) {
                grid.circle(x, y, r=max(buffer) * dot_scalar,
                            gp=gpar(fill="black", col="white"))
                grid.text(sprintf("%.2f", max(int)), x, y, 
                          gp=gpar(fontsize=10, fontface="bold", col="white"))
              }
           )
    )
  )
  draw(a, merge_legend=T, annotation_legend_list=lgd_list)
}


# Main ------------------------------------------
# DE file locations 
result = fread(paste(prefix, "csv", sep="."))
f_names = sub(".Milind_etal", "", nameset)
replacement = f_names

# Fisher exact testing for enrichment of DE sets -----------------
# generate list of gene sets for AD genes
# DE genes (all)
de_list = lapply(f_names, function(i) {
  result[padj < 0.05, GENE] # geneoverlap errors with only one list, so run it with the up/down
})
# names(de_list) = c(replacement, paste0(replacement, "_"))

# DE genes (split up/down)
up_list = lapply(f_names, function(i) {
  result[padj < 0.05 & beta > 0, GENE]
})
down_list = lapply(f_names, function(i) {
  result[padj < 0.05 & beta < 0, GENE]
})

# DE genes (split up/down)
qs = quantile(result[padj < 0.05, beta], c(0.1, 0.9))
top10_list = lapply(f_names, function(i) {
  result[padj < 0.05 & beta > qs[[2]], GENE]
})
bottom10_list = lapply(f_names, function(i) {
  result[padj < 0.05 & beta < qs[[1]], GENE]
})

# join and name
updown_list = c(top10_list, up_list, de_list, down_list, bottom10_list)
f_names_updown = c("Top 10%", "Up", "All", "Down", "Bottom 10%")
names(updown_list) = paste(replacement, f_names_updown)

# Milind_etal gene sets -------------------------
tmt = fread(paste0(refdir, "/Milind_etal.Human.Genes.txt"))

# limit to background of SOMALOGIC detectable proteins
tmt = tmt[GENE %in% result$GENE]

tmt_anno = fread(paste0(refdir, "/Milind_etal.group_annots.txt"))

tmt_list = sapply(unique(tmt$module), function(k) {
  k = list(unique(tmt[module == k , GENE]))
}, USE.NAMES=T)

##############
# background DE genes list
gs = length(unique(c(result[,GENE]))) # only the Somalogic set are detectable

# # enrich for all DE genes
# enrich(de_list, tmt_list, gs, "Modules_allDE", width=24, height=8)

# enrich for up/down genes
enrich(updown_list, tmt_list, gs, "Modules_Up-Down", width=24, height=8)





