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
library(GGally)
library(RColorBrewer)
library(data.table)
library(patchwork)
# library(UpSetR)
library(ComplexUpset)
# library(cowplot)

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


# Main ------------------------------------------
# DE file locations
fh <- file.path("../..", "baseline_PCs_relatedness", input,
  paste("dataProt_SS-205063.PC.related", input, "csv", sep = ".")
)
result <- lapply(fh, fread)
setattr(result, 'names', input)

# merge
result <- rbindlist(result)
# bottom <- round(min(result[["t stat"]], na.rm = TRUE) - 0.05, digits = 2)
# top <- round(max(result[["t stat"]], na.rm = TRUE) + 0.05, digits = 2)
# # drop APOE from input so it will match the results
input <- gsub(".APOE", "", input)
result[, upset := ifelse(padj < 0.05, 1, 0)]

# calculate genomic inflation
lambda.c <- function(pval) {
  median(qchisq(pval, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
}
group_data <- result[, .(lambda = lambda.c(pval)), by = dx]
group_data <- merge(group_data, result[padj < 0.05, .N, by = dx])

# Use n equally spaced breaks to assign each value to n-1 equal sized bins 
ii <- cut(group_data$lambda, 
  breaks = seq(min(group_data$lambda), max(group_data$lambda), len = 256), 
  include.lowest = TRUE
)
# Use bin indices, ii, to select color from vector of n-1 equally spaced colors
group_data$colors <- colorRampPalette(brewer.pal(9, "BuPu"))(255)[ii]


# plot
p <- upset(
  data = dcast(result, NCBI_ID ~ dx, value.var = "upset"), 
  intersect = input, 
  name = "overlapping proteins", 
  queries = lapply(seq(nrow(group_data)), function(i) {
    upset_query(
      set = group_data$dx[i],
      fill = group_data$colors[i]
    )
  }),
  width_ratio = 0.1,
  min_size = 3,
  min_degree = 1,
  keep_empty_groups = TRUE
)

ggsave(
  filename = file.path(output, paste(nameset, "overlaps.pdf", sep=".")), 
  plot = p, 
  device = cairo_pdf,
  width = 16, 
  height = 9
)


# Pairwise correlation plots -----------------------------------------------
upper_fn <- function(data, mapping, method = "s", use = "pairwise", ...){
  # grab data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)

  # calculate correlation
  corr <- cor(x, y, method = method, use = use)

  # calculate colour based on correlation value
  colFn <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)
  fill <- colFn[findInterval(corr, seq(-1, 1, length = 255))]

  ggally_cor(data = data, mapping = mapping, title = "\U03C1", size = 6, ...) + 
    theme(
      panel.background = element_rect(fill = fill, color = NA),
      panel.grid.major = element_blank()
    ) 
}

lower_fn <- function(data, mapping, method = "lm", ...) {
  ggplot(data = data, mapping = mapping, ...) +
    geom_point(colour = "#2166ac", size = 0.5) +
    geom_smooth(method = method, color = "black", ...) + 
    theme_classic()
}

p <- ggpairs(
  data = dcast(result, NCBI_ID ~ dx, value.var = "t stat"),
  upper = list(continuous = upper_fn),
  diag = list(continuous = function(...) ggally_densityDiag(...) + theme_linedraw()),
  lower = list(continuous = lower_fn),
  columns = 2:9
) 

ggsave(p,
  filename = file.path(output, paste(nameset, "correlations.pdf", sep=".")), 
  device = cairo_pdf,
  width = 12, 
  height = 12
)
ggsave(p,
  filename = file.path(output, paste(nameset, "correlations.png", sep=".")), 
  width = 12, 
  height = 12,
  units = "in", 
  dpi = 144, 
  type = "cairo"
)
