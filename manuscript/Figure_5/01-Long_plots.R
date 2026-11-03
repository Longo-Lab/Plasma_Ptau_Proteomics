#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 01-Long_plots.R
#'
#' Version: 0.0.2
#'
#' Purpose of script: Plots of proteomic score distributions & diagnosis
#'
#' Author: Robert B Butler III
#'
#' Date Created: 2025-10-20
#'
#' Copyright (c) 2025
#' Email: rrbutler@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'   v0.0.2 - incorporate splitting of plots for more than 15 scores
#'
#' ---------------------------
#' load up the packages we will need:  (uncomment as required)
library(optparse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(data.table)
library(patchwork)
library(ggpmisc)
library(ggridges)
library(paletteer)
library(stringr)

# read in args
option_list <- list(
  make_option( 
    c("-i", "--input"),
    type = "character",
    help = "List of input score columns comma separated"
  ),
  make_option(
    c("-s", "--score"),
    type = "character",
    help = "Scoretype to analyse"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project [default: Nameset of input]"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = ".",
    help = "File path to output directory [default: %default]"
  ),
  make_option(
    c("-w", "--box-width"),
    default = 14,
    dest = "box_width",
    help = "File width for box plots in inches [default: %default]"
  ),
  make_option(
    c("-b", "--batch-size"),
    default = 15,
    dest = "batch",
    help = "For many scores, batch size for splitting plots [default: %default]"
  ),
  make_option(
    c("-d", "--data-path"),
    type = "character",
    dest = "data_path",
    default = "../../module_scores",
    help = "File path to module scores [default: %default]"
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
  "score",
  "nameset",
  "output",
  "box_width",
  "batch",
  "data_path"
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
#' box plot with anova + tukey can facet 
#'
#' @param dt data.table 
#' @param x string, column name to plot on x
#' @param y string, column name to plot on y
#' @param facet string, column name for faceting
#' @param title main title
#' ...: additional objects passed to ggboxplot
#'
#' returns: ggobject
plot_box = function (dt, x, y, group, title = NULL, ...) {
  ggboxplot(dt, 
    x = x, 
    y = y, 
    # facet.by = facet, 
    add = "jitter", 
    col = group,
    add.params = list(), 
    palette = get_palette("BuPu", 4)[-1],
    na.rm = TRUE,
    # scales = "free_y", 
    # nrow = 1,
    ...
  ) + 
  geom_pwc(
    aes(group = .data[[group]]), 
    label = "p.format",
    p.adjust.method = "bonferroni",
    tip.length = 0,
    bracket.nudge.y = -0.3
  ) +
  ggtitle(title) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
}

#' Plot Density Ridges for Score Distribution
#'
#' Creates a ggplot2 density ridges plot showing the distribution of scores
#' across different categories defined by the x variable. If there are more
#' than 15 categories, splits into batches and combines horizontally.
#'
#' @param df A data frame containing the data to be plotted.
#' @param x A character string specifying the column name in df to be used for
#'        the x-axis values (continuous distribution).
#' @param cols Optional. A vector of custom colors to add on. If NULL, the default
#'        palette will be used.
#' @param pal A character string specifying the paletteer palette to use. 
#'        Default is "colorBlindness::paletteMartin".
#' @param batch_size Integer. Maximum number of categories per plot. Default is 15.
#' @param score_labels Optional. A named vector of custom labels for the score categories.
#'        Names should match values in the Score column, values are the display labels.
#'
#' @return A ggplot2 object or patchwork object representing density ridges plot(s).
#'
#' @import ggplot2
#' @import ggridges
#' @import paletteer
#' @import patchwork
#' @import dplyr
#' @import stringr
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' plot_comp(my_dataframe, "CellType")
#'
#' # With custom labels
#' labels <- c("T_cell" = "T Cell", "B_cell" = "B Cell", "NK_cell" = "NK Cell")
#' plot_comp(my_dataframe, "CellType", score_labels = labels)
#' }
#'
#' @export
plot_comp <- function(df, x, cols = NULL, pal = "colorBlindness::paletteMartin", 
                      batch_size = 15, score_labels = NULL) {
  
  # Apply custom labels if provided
  if (!is.null(score_labels)) {
    df <- df %>%
      dplyr::mutate(Score = dplyr::recode(.data[["Score"]], !!!score_labels))
  }
  
  # Get unique values of the Score column (categories to be batched)
  unique_scores <- unique(df[["Score"]])
  n_scores <- length(unique_scores)
  
  # If <= batch_size, create single plot
  if (n_scores <= batch_size) {
    return(create_single_ridge_plot(df, x, cols, pal, show_labels = TRUE))
  }
  
  # Split into batches
  batches <- split(unique_scores, ceiling(seq_along(unique_scores) / batch_size))
  
  # Create plots with consistent dimensions
  plots <- list()
  for (i in seq_along(batches)) {
    batch_data <- df[df[["Score"]] %in% batches[[i]], ]
    plots[[i]] <- create_single_ridge_plot(batch_data, x, cols, pal, show_labels = FALSE) +
      theme(plot.margin = margin(t = 10, r = 5, b = 10, l = 5, unit = "pt"))
  }
  
  # Simple patchwork combination
  require(patchwork)
  combined_plot <- wrap_plots(plots, 
    nrow = 1, 
    widths = rep(1, length(plots)),
    heights = c(rep(batch_size, length(plots) - 1), 3)
  )
  
  final_plot <- combined_plot + 
    plot_annotation(
      title = "Score Distribution by Type",
      caption = paste("Values:", x),
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.caption = element_text(hjust = 0.5, size = 12)
      )
    )
  
  return(final_plot)
}
#' Helper function to create individual ridge plots
#' @keywords internal
create_single_ridge_plot <- function(df, x, cols = NULL, pal = "colorBlindness::paletteMartin", show_labels = TRUE) {
  
  # Create the plot
  p <- ggplot(df, aes(x = .data[[x]], y = Score, fill = Score, color = Score)) +
    geom_density_ridges(alpha = 0.7, scale = 0.9, rel_min_height = 0.01) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 26)) +
    theme_ridges() +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Add labels only if show_labels is TRUE
  if (show_labels) {
    p <- p + labs(
      title = "Score Distribution by Type",
      x = paste("Values:", x),
      y = "Score Type"
    )
  } else {
    p <- p + labs(
      title = NULL,
      x = NULL,
      y = "Score Type"
    )
  }
  
  # Get colors from paletteer
  require(paletteer)
  
  # Add colors
  if (!is.null(cols)) {
    base_colors <- paletteer_d(pal)
    all_colors <- c(base_colors, cols)
    p <- p + scale_color_manual(values = all_colors)
    p <- p + scale_fill_manual(values = all_colors)
  } else {
    p <- p + scale_color_paletteer_d(pal)
    p <- p + scale_fill_paletteer_d(pal)
  }
  
  return(p)
}

#' Scatterplot with a trendline and stats, can facet
#'
#' @param dt data.table 
#' @param x string, column name to plot on x
#' @param y string, column name to plot on y
#' @param facet string, column name for faceting
#' @param group string, column name for color groups
#' @param title main title
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param trendline_type string, either "group" for separate trendlines per group
#'   or "overall" for a single trendline per facet (default: "group")
#' @param ... additional objects passed to ggplot layers
#'
#' @return ggplot object
plot_corr = function(dt, x, y, facet, group, title = NULL, xlab = NULL, ylab = NULL,
                     trendline_type = "group") {
  
  # Validate trendline_type parameter
  if(!trendline_type %in% c("group", "overall")) {
    stop("trendline_type must be either 'group' or 'overall'")
  }
  
  # Create the base plot
  p <- ggplot(dt, aes_string(x = x, y = y, fill = group, color = group)) +
    geom_point() +
    scale_color_manual(
      values = get_palette("BuPu", 4)[-1], 
      aesthetics = c("color", "fill")
    ) +
    labs(title = title, x = xlab, y = ylab) +
    facet_wrap(c(facet), nrow = 2, scales = "free_y") +
    theme_pubr() +
    # Add custom legend with larger points
    guides(
      color = guide_legend(override.aes = list(size = 5)),
      fill = guide_legend(override.aes = list(size = 5))
    )
  
  # Add trendlines and statistics based on the trendline_type
  if(trendline_type == "group") {
    # Group-specific trendlines (keeps the grouping aesthetic for stat_smooth)
    p <- p + 
      stat_smooth(method = "lm", na.rm = TRUE, show.legend = FALSE) + 
      stat_poly_eq(aes(label = paste(..adj.rr.label.., ..p.value.label.., sep = "~~~~")))
  } else {
    # Overall trendline (ignores grouping for stat_smooth)
    p <- p + 
      stat_smooth(
        method = "lm", 
        na.rm = TRUE, 
        show.legend = FALSE,
        aes(color = NULL, fill = NULL), # Remove group aesthetic
        color = "#525252" # Set a neutral color
      ) +  
      stat_poly_eq(
        aes(color = NULL, fill = NULL, # Remove group aesthetic
            label = paste(..adj.rr.label.., ..p.value.label.., sep = "~~~~")),
        color = "#525252"  # Set a neutral color
      )
  }
  
  return(p)
}


# Main ------------------------------------------
# load overall base file
result <- fread("../../dataProt_SS-205063.merged.QC.csv")
a <- result[, 1:27]

# Map in PCs assumes your primary data table exists 
x <- fread("../../genomic/mypcair.PCs.csv", colClasses = c(PIDN = "character"))
a <- merge(a, x, by = "PIDN")

# setup for eigengenes vs NES scores
if (score %like% "eigengenes") {
  tail <- "csv" 
  score_val <- "Eigengene"
  scores <- fread(file.path(data_path, "STRING_merged_scaled_scores.csv"))
  scores <- unique(scores[, .(Biomarker, cluster, term)])
  scores[, name := paste0("ME", Biomarker, "-", cluster)]
  scor <- scores[, term]
  names(scor) <- scores[, name]
} else if (score %like% "WSS") {
  tail <- "csv" 
  score_val <- "Weighted Sum of Scores"
  scores <- fread(file.path(data_path, "STRING_merged_scaled_scores.csv"))
  scores <- unique(scores[, .(Biomarker, cluster, term)])
  scores[, name := paste0(Biomarker, ".", cluster)]
  scor <- scores[, term]
  names(scor) <- scores[, name]
} else {
  tail <- "normalized_scores.csv"
  score_val <- "NES"
  scor <- NULL
}

# Map in PCs assumes your primary data table exists 
x <- fread(
  file.path(data_path, 
    paste("dataProt_SS-205063", score, tail, sep = ".")
  ), 
  colClasses = c(SampleID = "character")
)
a <- merge(a, x, by = "SampleID")

# drop outlier samples
drop = a[abs(ConnectivityZscore) > 10, SampleID]
a = a[!SampleID %in% drop]

# only APOE
a <- a[APOE != ""]

# Add in phenotypes for individuals -----------------
x <- data.table(
  readxl::read_xlsx(
    "../../../ADRC traits/Z1XADRCNACC_DATA_LABELS_2021_10_18_matchMetadata_allSamples.xlsx",
    sheet = "Z1XADRCNACC_DATA_LABELS_2021_10"
  )
)

a <- merge(a, 
  x[, .(Sample_Id, CDRSUM, MOCATOTS)], 
  by.x = "SampleID", 
  by.y = "Sample_Id", 
  all.x = TRUE, 
  sort = FALSE
)
# Three individuals had a MOCATOTS of 88, which is possibly a missing value
a[MOCATOTS > 30, MOCATOTS := NA]

# Setup for baseline visit only analysis ---------------------------
# get the first visit for each score, with Diagnosis_group, Age, Gender, APOE, HybControlNormScale
b <- lapply(input, function(i) {
  x <- a[!is.na(get(i))][, value := get(i)]
  x <- x[, .(PIDN, DS1V, value, Diagnosis_group, Age, Gender, APOE, HybControlNormScale, CDRSUM, MOCATOTS)]
  # keep only first visit where plasma available
  x <- x[order(DS1V)][, head(.SD, 1), by = "PIDN"]
})
setattr(b, 'names', input)
b <- rbindlist(b, idcol = 'Score')
# b[, .N, by = Score]
b[, Score := factor(Score, levels = sort(unique(Score)))]
setnames(b, "value", score_val)

# for this plot, reclassify MCI-AD as MCI
b[, Diagnosis_group := gsub("MCI-AD", "MCI", Diagnosis_group)]
b$Diagnosis_group <- factor(b$Diagnosis_group, levels = c("HC", "MCI", "AD"))


# plot vs diagnosis boxplot ----------------------------------------
p1 = plot_box(b, 
  x = "Score", 
  y = score_val, 
  group = "Diagnosis_group",
  xlab=paste("Diagnosis"),
  ylab=paste(score_val)
)
ggsave(
  filename = file.path(output, paste(nameset, score, "Diag_box.pdf", sep=".")), 
  plot = p1, 
  width = box_width, 
  height = box_width, 
  device = cairo_pdf
)
ggsave(
  filename = file.path(output, paste(nameset, score, "Diag_box.png", sep=".")), 
  plot = p1, 
  width = box_width, 
  height = box_width, 
  dpi = 144, 
  type = cairo
)


# maybe distribution plots ----------------------------------------
# plot_comp(b, "NES")
extra_colors = paletteer_d("ggsci::planetexpress_futurama") 
p1 <- plot_comp(
  df = b, 
  x = score_val, 
  cols = extra_colors, 
  score_labels = scor,
  batch_size = batch
)

ggsave(
  filename = file.path(output, paste(nameset, score, "Dist_joy.pdf", sep=".")), 
  plot = p1, 
  width = box_width, 
  height = 8, 
  device = cairo_pdf
)
ggsave(
  filename = file.path(output, paste(nameset, score, "Dist_joy.png", sep=".")), 
  plot = p1, 
  width = box_width, 
  height = 10, 
  dpi = 144, 
  type = cairo
)

# plot vs Age (consider MOCA or memory composite) scatter w/ overall trend -----
# scatterplot of L2FCs
lapply(c("Age", "CDRSUM", "MOCATOTS"), function(i) {
  p1 <- plot_corr(b, 
    x = i, 
    y = score_val, 
    facet = "Score",
    group = "Diagnosis_group",
    xlab = paste(i, "at visit"),
    ylab = paste("Score", score_val),
    trendline_type = "overall"
  )
  ggsave(
    filename = file.path(
      output, paste(nameset, score,  i, "overall_reg.pdf", sep=".")
    ), 
    plot = p1, 
    width = box_width, 
    height = 7, 
    device = cairo_pdf
  )
  ggsave(
    filename = file.path(
      output, paste(nameset, score, i, "overall_reg.png", sep=".")
    ), 
    plot = p1, 
    width = box_width, 
    height = 7, 
    dpi = 144, 
    type = cairo
  )
})
