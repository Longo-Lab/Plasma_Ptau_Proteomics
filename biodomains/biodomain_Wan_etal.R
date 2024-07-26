#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: biodomain_Wan_etal.R
#'
#' Version: 0.0.2
#'
#' Purpose of script: Functions using Wan, et al. modules
#'
#' Author: Robert B Butler III; Crystal Han
#'
#' Date Created: 2022-06-12
#'
#' Copyright (c) 2022
#' Email: rrbutler@stanford.edu; cyouh95@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'   v0.0.2 - updated for relative paths
#'
#' ---------------------------

#' load up the packages we will need:  (uncomment as required)

library(data.table)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

#' ---------------------------

heat = function(mat, pval, int, name, ...) {
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
  col_fun = colorRamp2(seq(min(mat, na.rm=T), max(max(mat, na.rm=T), 1), length=9), 
                       brewer.pal(name="BuPu", n=9))
  col_fun2 = brewer.pal(name="Dark2", n=7)
  # column annotations
  ad = get_ad_file(...)
  
  ad_anno = unique(ad, by=c("module", "Consensus_Cluster"))[, Consensus_Cluster]
  ad_tiss = unique(ad, by=c("module", "Brain_Region"))[, Brain_Region]
  names(col_fun2) = unique(ad_tiss)
  ha = HeatmapAnnotation("Consensus Cluster"=ad_anno, 
                         show_legend=F,
                         col=list("Consensus Cluster"=c("Astroglial"="#FAFF34", 
                                                        "Microglial"="#62B0B0", 
                                                        "Neuronal"="#7E2221", 
                                                        "Oligodendroglial"="#6CC835", 
                                                        "Proteostasis"="#280808")))
  hb = HeatmapAnnotation("Brain Region"=ad_tiss, col=list("Brain Region"=col_fun2))
  # scale max int (for circle size) to matrix "buffer"
  buffer = ((int/max(int))/4 + .25)
  buffer[is.nan(buffer)] <- 0
  dot_scalar = unit(1.5, "cm")
  # heatmap plot 
  a = Heatmap(mat, cluster_rows=F, cluster_columns=F, na_col="#FFFFFF", 
              col=col_fun, name=name, column_names_rot=60, row_names_side="left",
              top_annotation=ha, column_split=ad_anno, column_gap = unit(5, "mm"), 
              column_title_gp = gpar(fontsize = 20, fontface = "bold"),
              rect_gp=gpar(type="none"), bottom_annotation=hb, 
              cell_fun=function(j, i, x, y, width, height, fill) { 
                grid.circle(x=x, y=y, r=(buffer[i, j] * dot_scalar), 
                            gp=gpar(fill=col_fun(mat[i, j]), 
                                    col=ifelse(pval[i, j] < -log10(0.05), NA, 
                                               "black")))
                if (pval[i, j] > -log10(0.05)) {
                  grid.text(sprintf("%.2f", int[i, j]), x, y, 
                            gp=gpar(fontsize=10, fontface="bold", 
                                    col=ifelse(mat[i, j] > quantile(mat, 0.97), 
                                               "white", "black")))
                }
              })
  # heatmap legend
  lgd_list = list(
    Legend(labels=c(rep("", 3)), 
           title = "Jaccard Index",
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
                grid.circle(x, y, r=quantile(buffer, 0.9) * dot_scalar,
                            gp=gpar(fill="black", col="white"))
                grid.text(sprintf("%.2f", quantile(int, 0.9)), x, y, 
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
  draw(a, merge_legend=T, heatmap_legend_list=lgd_list, padding = unit(c(2, 2, 2, 12), "mm"))
}

get_ad_file <- function(path) {
  ad = fread(file.path(path, "Wan_etal.Human.Genes.txt"))
  
  ad[, Consensus_Cluster := fcase(
    Consensus_Cluster == 'A', 'Astroglial',
    Consensus_Cluster == 'B', 'Microglial',
    Consensus_Cluster == 'C', 'Neuronal',
    Consensus_Cluster == 'D', 'Oligodendroglial',
    Consensus_Cluster == 'E', 'Proteostasis',
    default = NA
  )]
  
  # order by Consensus Cluster (for figure coloring order)
  setorder(ad, Consensus_Cluster, module)
  
  ad
}

get_modules_list <- function(...) {
  ad = get_ad_file(...)
  
  sapply(unique(ad$module), function(k) {
    k = list(unique(ad[module == k , GENE]))
  }, USE.NAMES=T)
}

# dimensions for saving figures
w = 24
dh = 3
h1 = 8
