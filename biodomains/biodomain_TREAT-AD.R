#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: biodomain_TREAT-AD.R
#'
#' Version: 0.0.3
#'
#' Purpose of script: Functions using TREAT-AD modules
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
#'   v0.0.3 - adding symbol lookup table & title (scripts #8)
#'   v0.0.2 - removed unnecessary libraries
#'   v0.0.2 - updated for relative paths
#'
#' ---------------------------

#' load up the packages we will need:  (uncomment as required)

library(data.table)
# library(stringr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

#' ---------------------------
# global variables
target <- "TREAT-AD.113.one2one.Human.Genes.txt"

# dimensions for saving figures
w = 12
dh = 0
h1 = 6


#' ---------------------------

#' ComplexHeatmap template for this target
#' overlaps colored by value, and sized by overlap (given matrices, name)
#'
#' @param mat matrix of heatmap values (fdr, p-val, or odds.ratio)
#' @param pval matched matrix of fdr values, to bold significant circles
#' @param int matched matrix of Jaccard indexes for intersection sizes
#' @param name heatscale name for plot (fdr, p-val, or odds.ratio)
#'
#' @return a draw graph call
heat = function(mat, pval, int, name, title, ...) {
  # resolve Inf to real number (round) and check for NAs
  mat[is.infinite(mat)] = -log10(2.225074e-308)
  pval[is.infinite(pval)] = -log10(2.225074e-308)
  if (any(is.na(int))) {
    print("Jaccard index contains NAs!")
    quit("no", 1)
  }

  # color palette
  col_fun = colorRamp2(seq(min(mat, na.rm=T), 40, length=9), 
                       brewer.pal(name="BuPu", n=9))
  # col_fun2 = brewer.pal(name="Dark2", n=7)
  # scale max int (for circle size) to matrix "buffer"
  buffer = ((int/max(int))/4 + .25)
  buffer[is.nan(buffer)] <- 0
  dot_scalar = unit(1.5, "cm")
  # heatmap plot 
  a <- Heatmap(mat, 
    cluster_rows = FALSE, 
    cluster_columns = FALSE, 
    na_col = "#FFFFFF", 
    col = col_fun, 
    name = name, 
    column_names_rot = 60, 
    row_names_side = "left",
    column_order = sort(colnames(mat)), 
    # top_annotation=ha, bottom_annotation=hb, column_split=ad_anno, 
    # column_gap = unit(5, "mm"), 
    column_title = paste("Enrichment of", title),
    column_title_gp = gpar(fontsize = 20, fontface = "bold"),
    rect_gp = gpar(type = "none"), 
    cell_fun = function(j, i, x, y, width, height, fill) { 
      grid.circle(
        x = x, 
        y = y, 
        r = (buffer[i, j] * dot_scalar), 
        gp = gpar(
          fill = col_fun(mat[i, j]), 
          col = ifelse(pval[i, j] < -log10(0.05), NA, "black")
        )
      )
      if (pval[i, j] > -log10(0.05)) {
        grid.text(
          label = sprintf("%.2f", int[i, j]), 
          x = x, 
          y = y, 
          gp = gpar(
            fontsize = 10, 
            fontface = "bold", 
            col = ifelse(mat[i, j] > quantile(mat, 0.97), "white", "black")
          )
        )
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
  draw(a, merge_legend=T, annotation_legend_list=lgd_list)
}

#' load gene modules from reference in to a list of named lists
#'
#' @param path string path to reference file
#'
#' @return list of named lists of gene ids, sorted
get_modules_list <- function(path) {
  ad = fread(file.path(path, target))
  
  sapply(unique(ad$module), function(k) {
    k = list(unique(ad[module == k , GENE]))
  }, USE.NAMES=T)
}

#' load gene symbol lookup table from reference
#'
#' @param path string path to reference file
#'
#' @return list of named lists of gene ids, sorted
get_lookup <- function(path) {
  ad = fread(file.path(path, target))
  return(unique(ad[, .(GENE, symbol)]))
}

