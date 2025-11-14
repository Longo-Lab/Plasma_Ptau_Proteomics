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
#'   v0.0.2 - 
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
    c("-d", "--data-path"),
    type = "character",
    dest = "data_path",
    default = "../../celltype_scores",
    help = "File path to celltype scores [default: %default]"
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
#' across different categories defined by the x variable.
#'
#' @param df A data frame containing the data to be plotted.
#' @param x A character string specifying the column name in df to be used for
#'        grouping on the y-axis of the ridge plot.
#' @param cols Optional. A vector of custom colors to add on. If NULL, the default
#'        palette will be used.
#' @param pal A character string specifying the paletteer palette to use. 
#'        Default is "colorBlindness::paletteMartin".
#' @param custom_color Optional. A single color to append to the palette.
#'
#' @return A ggplot2 object representing a density ridges plot.
#'
#' @import ggplot2
#' @import ggridges
#' @import paletteer
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' plot_comp(my_dataframe, "CellType")
#'
#' # With custom palette
#' plot_comp(my_dataframe, "CellType", pal = "viridis::viridis")
#'
#' # With custom color added to palette
#' plot_comp(my_dataframe, "CellType", custom_color = "#FF5733")
#' }
#'
#' @export
plot_comp <- function(df, x, cols = NULL, pal = "colorBlindness::paletteMartin") {
  # plot
  p <- ggplot(df, aes(x = .data[[x]], y = Score, fill = Score, color = Score)) +
    geom_density_ridges(alpha = 0.7, scale = 0.9, rel_min_height = 0.01) +
    theme_ridges() +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = "NES Distribution by Cell Type",
      x = paste("Specificity", x),
      y = "Cell Type"
    )
  
  # Get colors from paletteer
  require(paletteer)
  base_colors <- paletteer_d(pal)
  
  # Add custom color if provided
  if (!is.null(cols)) {
    all_colors <- c(base_colors, cols)
    p <- p + scale_color_manual(values = all_colors)
    p <- p + scale_fill_manual(values = all_colors)
  } else {
    p <- p + scale_color_paletteer(palette_name)
    p <- p + scale_fill_paletteer(palette_name)
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

# Map in PCs assumes your primary data table exists 
if (score %like% "eigengenes") tail <- "csv" else tail <- "normalized_scores.csv"
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
setnames(b, "value", "NES")

# for this plot, reclassify MCI-AD as MCI
b[, Diagnosis_group := gsub("MCI-AD", "MCI", Diagnosis_group)]
b$Diagnosis_group <- factor(b$Diagnosis_group, levels = c("HC", "MCI", "AD"))


# plot vs diagnosis boxplot ----------------------------------------
p1 = plot_box(b, 
  x = "Score", 
  y = "NES", 
  group = "Diagnosis_group",
  xlab=paste("Diagnosis"),
  ylab=paste("Normalized Enrichment Score")
)
ggsave(
  filename = file.path(output, paste(nameset, score, "Diag_box.pdf", sep=".")), 
  plot = p1, 
  width = box_width, 
  height = 7, 
  device = cairo_pdf
)
ggsave(
  filename = file.path(output, paste(nameset, score, "Diag_box.png", sep=".")), 
  plot = p1, 
  width = box_width, 
  height = 7, 
  dpi = 144, 
  type = cairo
)


# maybe distribution plots ----------------------------------------
# plot_comp(b, "NES")
extra_colors = paletteer_d("ggsci::planetexpress_futurama") 
ggsave(
  filename = file.path(output, paste(nameset, score, "Dist_joy.pdf", sep=".")), 
  plot = plot_comp(b, "NES", cols = extra_colors), 
  width = 6, 
  height = box_width, 
  device = cairo_pdf
)
ggsave(
  filename = file.path(output, paste(nameset, score, "Dist_joy.png", sep=".")), 
  plot = plot_comp(b, "NES", cols = extra_colors), 
  width = 6, 
  height = box_width, 
  dpi = 144, 
  type = cairo
)

# plot vs Age (consider MOCA or memory composite) scatter w/ overall trend -----
# scatterplot of L2FCs
lapply(c("Age", "CDRSUM", "MOCATOTS"), function(i) {
  p1 <- plot_corr(b, 
    x = i, 
    y = "NES", 
    facet = "Score",
    group = "Diagnosis_group",
    xlab = paste(i, "at visit"),
    ylab = paste("Score NES"),
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
