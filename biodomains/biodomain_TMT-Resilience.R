#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: biodomain_TMT-Resilience.R
#'
#' Version: 0.0.3
#'
#' Purpose of script: Functions using TMT-Resilience modules
#'
#' Author: Robert B Butler III; Crystal Han
#'
#' Date Created: 2023-10-09
#'
#' Copyright (c) 2023
#' Email: rrbutler@stanford.edu; cyouh95@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'   v0.0.3 - adding symbol lookup table & title (scripts #8)
#'   v0.0.2 - updated for relative paths
#'   
#'
#' ---------------------------

#' load up the packages we will need:  (uncomment as required)

library(data.table)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

#' ---------------------------
# global variables
target <- "TMT-Resilience.113.one2one.Human.Genes.txt"

# dimensions for saving figures
w = 24
dh = 1
h1 = 7.5


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
  
  # column annotations
  ad = get_ad_file(...)
  
  # upper annot
  col = list(
    "Greater resilience" = "olivedrab", 
    "Less resilience" = "tomato3"
  )
  ad_anno <- unique(ad, by = c("module", "cluster"))[, cluster]
  ad_ind <- lapply(names(col), function(i) which(ad_anno %in% i))
  names(ad_ind) <- names(col)
  panel_fun = function(index, nm) {
    grid.rect(gp = gpar(fill = col[[nm]], col = col[[nm]]), height = unit(1, "cm"))
    grid.text(
      label = nm, 
      x = 0.5, 
      y = 0.5, 
      gp = gpar(col = "white", fontsize = 16, fontface = "bold")
    )
  }
  ha <- HeatmapAnnotation(
    "Cognitive resilience" = anno_block(
      align_to = ad_ind, 
      panel_fun = panel_fun #, 
      # show_name = TRUE
    )
  )

  # color palette
  col_fun = colorRamp2(seq(min(mat, na.rm=T), max(max(mat, na.rm=T), 1), length=9), 
                       brewer.pal(name="BuPu", n=9))
  # scale max int (for circle size) to matrix "buffer"
  buffer = ((int/max(int))/4 + .25)
  buffer[is.nan(buffer)] <- 0
  dot_scalar = unit(1.5, "cm")
  
  # adding padding for the title
  ht_opt$TITLE_PADDING = unit(c(0.75, 0, 0.75, 0), "cm")
  
  # heatmap plot 
  a <- Heatmap(mat, 
    cluster_rows = FALSE, 
    cluster_columns = FALSE, 
    na_col = "#FFFFFF", 
    col = col_fun, 
    name = name, 
    column_names_rot = 60, 
    row_names_side = "left",
    top_annotation = ha, 
    column_split = ad_anno, 
    column_gap = unit(5, "mm"), 
    column_title = paste("Enrichment of", title),
    column_title_gp = gpar(fontsize = 22, fontface = "bold"),
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
  draw(a, merge_legend=T, annotation_legend_list=lgd_list, padding = unit(c(2, 2, 2, 14), "mm"))
}

#' load gene modules from reference and format annotations
#'
#' @param path string path to reference file
#'
#' @return list of named lists of gene ids, sorted
get_ad_file <- function(path) {
  ad = fread(file.path(path, target))
  
  ad[, cluster := ' ']
  ad[module %in% c('M22 Synapse', 'M5 Synapse', 'M36 Exocytosis', 'M30 Mitochondria/ER'), cluster := 'Greater resilience']
  ad[module %in% c('M11 Proteosome', 'M15 MAPK Signaling', 'M32 GPCR Signaling', 'M16 Gluconeogensis'), cluster := 'Less resilience']
  
  ad$cluster = factor(ad$cluster, levels = c('Greater resilience', ' ', 'Less resilience'))
  ad$module = factor(ad$module, levels = gtools::mixedsort(unique(ad$module)))
  
  setorder(ad, cluster, module)
  
  ad
}

#' load gene modules from reference in to a list of named lists
#'
#' @param path string path to reference file
#'
#' @return list of named lists of gene ids, sorted
get_modules_list <- function(...) {
  ad = get_ad_file(...)
  
  modules <- unique(ad$module)
  names(modules) <- unique(ad$module)
  
  sapply(modules, function(k) {
    k = list(unique(ad[module == k, GENE]))
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

