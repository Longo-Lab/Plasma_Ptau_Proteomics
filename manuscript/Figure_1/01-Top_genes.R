#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 01-Top_genes.R
#'
#' Version: 0.0.4
#'
#' Purpose of script: Tops hits for biomarker associations in Fig 1
#'
#' Author: Robert B Butler III, Michael Peter Brown
#'
#' Date Created: 2024-09-11
#'
#' Copyright (c) 2024
#' Email: rrbutler@stanford.edu, mpb@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'
#'
#' ---------------------------
#' load up the packages we will need:  (uncomment as required)
library(optparse)
library(tools)
library(ggplot2)
library(RColorBrewer)
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
library(data.table)
library(patchwork)
# library(Cairo)


# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "List of input datasets comma separated"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = ".",
    help = "File path to output directory [default: %default]"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project [default: Nameset of input]"
  )
)

args_list <- parse_args(OptionParser(option_list = option_list))

# Nameset default case
if(is.null(args_list$nameset)){
  args_list$nameset <- file_path_sans_ext(basename(args_list$input))
}

# read in fields
fields <- c(
  "input",
  "output",
  "nameset"
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

# input variables
input <- strsplit(input, ",")[[1]]


# Functions ---------------------------
#' plot heatmap of overlaps colored by value, and sized by overlap (given matrices, name)
#'
#' @param mat matrix of heatmap values (fdr, p-val, or odds.ratio)
#' @param pval matched matrix of fdr values, to bold significant circles
#' @param count matched matrix of gene count frequencies for asterisks
#' @param name legend name for plot (fdr, p-val, or odds.ratio)
#' @param range vector of two numbers for upper/lower bounds of heatmap default: auto select
#'
#' @return a draw graph call
heat <- function(mat, pval, count, name, range = NULL, ...) {
  # color palette
  if (is.null(range)) {
    bottom <- floor(min(mat, na.rm = TRUE))
    top <- ceiling(max(mat, na.rm = TRUE))
  } else {
    bottom <- range[1]
    top <- range[2]
  }
  col_fun <- colorRamp2(c(bottom, 0, top), c("#2166ac", "#FFFFFF", "#b2182b"))

  # heatmap plot 
  a <- Heatmap(mat, 
    cluster_rows = FALSE, 
    cluster_columns = FALSE, 
    na_col = "#FFFFFF", 
    col = col_fun, 
    name = name, 
    show_column_names = FALSE, 
    heatmap_legend_param = list(at = c(bottom, 0, top)), 
    cell_fun = function(j, i, x, y, w = unit(1.5, "cm"), h, col) {
      # add asterisks text to each grid spot based on count freq matrix
      grid.text(paste(rep("*", count[i, j]), collapse = ""), x, y, 
        gp = gpar(fontsize = 14, fontface = "bold", col = "black")
      )
      # bold box if padj < 0.05
      grid.rect(x = x, y = y, width = w, height = h, 
        gp = gpar(
          lwd = 3, 
          fill = 'transparent', 
          col = ifelse(pval[i, j] < 0.05, "black", NA)
        )
      )
    }, 
    ...
  )
  draw(a)
}

#' gets matricies for complex heatmap to use
#'
#' @param dt data.table
#' @param colname1 string, colname for first matrix (often log2fc)
#' @param colname2 string, colname for second matrix (often padj), values sorted by this column
#' @param colname3 string, colname for third matrix (count frequency)
#'
#' @return a list of two matrices
get_mats <- function(dt, rownam, colname1, colname2, colname3){
  dt1 <- dt[get(colname1) > 0][order(get(colname2))][1:20]
  dt2 <- dt[get(colname1) < 0][order(get(colname2))][1:20][order(-get(colname2))]
  dt <- rbindlist(list(dt1, dt2))
  a <- as.matrix(dt[, c(mget(rownam), mget(colname1))], rownames = rownam)
  b <- as.matrix(dt[, c(mget(rownam), mget(colname2))], rownames = rownam)
  c <- as.matrix(dt[, c(mget(rownam), mget(colname3))], rownames = rownam)
  return(list(a, b, c))
}

#' Make a volcano with labels & cutoffs
#'
#' @param x data.table, results 
#' @param titl string, label for graph
#' @param thresh numeric, cutoff for adjusted p-val
#' @param lfc numeric, log2 fold change cutoff for significance
#'
#' @returns ggplot object
my_volcano <- function(x, titl = NULL, sub = NULL, thresh = 0.05, lfc = 0) {
  require("ggplot2")
  require("EnhancedVolcano")
  # limit to pval subset
  x <- x[!is.na(x[["padj"]]), ]
  x[padj < 1e-150, padj := 1e-150]
  # color by up/downregulated
  keyvals <- ifelse(
    x$beta < -lfc & x[["padj"]] < thresh, 
    "#2d004b",
    ifelse(x$beta > lfc & x[["padj"]] < thresh, "#7f3b08", "grey45")
  )
  names(keyvals)[keyvals == "#7f3b08"] <- "upregulated"
  names(keyvals)[keyvals == "grey45"] <- "not significant"
  names(keyvals)[keyvals == "#2d004b"] <- "downregulated"
  # plot with EnhancedVolcano
  p1 <- EnhancedVolcano(x,
    x = "beta",
    y = "pval",
    lab = x[["Target"]],
    selectLab = x[abs(beta) > lfc & padj < thresh, Target],
    labFace = "bold",
    drawConnectors = TRUE,
    widthConnectors = 0.25,
    max.overlaps = 40,
    arrowheads = FALSE, 
    # xlim = c(-9,9),
    xlab = "beta",
    # ylab = bquote(~-Log[10] ~ italic(P)),
    title = titl,
    subtitle = sub,
    pCutoff = NA,
    FCcutoff = lfc,
    cutoffLineCol = "light gray",
    pointSize = 1.5,
    labSize = 3.5,
    colCustom = keyvals,
    colAlpha = 1,
    legendPosition = 'none',
    legendLabSize = 11,
    legendIconSize = 5,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    border = "full",
    borderWidth = 0.5,
    borderColour = "black"
  ) + theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
  p1
}


# Main ------------------------------------------
# DE file locations
fh <- file.path("../..", "baseline_PCs_relatedness", input,
  paste("dataProt_SS-205063.PC.related", input, "csv", sep = ".")
)
result <- lapply(fh, fread)
setattr(result, 'names', input)

# merge
result <- rbindlist(result)
# count of significant occurences across biomarkers
result[, count := sum(.SD < 0.05), by = NCBI_ID, .SDcols = "padj"]

# drop APOE from input so it will match the results
input <- gsub(".APOE", "", input)

# plot colored by beta and t stat
lapply(c("beta", "t stat"), function(j) {
  # select range across all sets
  bottom <- round(min(result[[j]], na.rm = TRUE) - 0.05, digits = 2)
  top <- round(max(result[[j]], na.rm = TRUE) + 0.05, digits = 2)

  # heatmap
  m.list <- lapply(input, function(i) {
    get_mats(result[dx == i], "Symbol", j, "padj", "count")
  })
  p.list <- lapply(seq(input), function(i) {
    grid.grabExpr(
      heat(
        mat = m.list[[i]][[1]], 
        pval = m.list[[i]][[2]], 
        count = m.list[[i]][[3]], 
        name = j, 
        range = c(bottom, top),
        column_title = gt_render(input[i], padding = unit(c(0, 0, 0, 0.5), "in")),
        column_title_gp = gpar(fontsize = 16, fontface = "bold")
      )
    )
  })


  ggsave(
    filename = file.path(output, 
      paste(nameset, gsub(" ", "_", j), "top_genes.pdf", sep=".")
    ), 
    plot = wrap_plots(p.list) + plot_layout(nrow = 1), 
    device = cairo_pdf,
    width = 20, 
    height = 14
  )
})

# Volcano printing ------------------------------------
# prep for printing
# options(ggrepel.max.overlaps = 150)
# CairoPDF(
pdf(
  file = file.path(output, paste(nameset, "volcanos.pdf", sep = ".")),
  height = 11,
  width = 8.5
)

# if file exists, make a volcano
lapply(input, function(i) {
  x <- result[dx == i]
  my_volcano(x, titl = i, sub = paste("N =", unique(x[, nobs])))
})

graphics.off()


# if file exists, make a volcano
p.list <- lapply(input, function(i) {
  x <- result[dx == i]
  my_volcano(x, titl = i, sub = paste("N =", unique(x[, nobs])))
})

ggsave(
  filename = file.path(output, paste(nameset, "volcanos.png", sep = ".")), 
  plot = wrap_plots(p.list) + plot_layout(nrow = 2), 
  width = 20, 
  height = 14,
  dpi = 150
)
