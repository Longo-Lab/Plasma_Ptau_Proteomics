#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 03-Long_plots.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Plots of biomarker characteristics in Fig 1
#'
#' Author: Robert B Butler III
#'
#' Date Created: 2025-02-24
#'
#' Copyright (c) 2025
#' Email: rrbutler@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'
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
plot_box = function (dt, x, y, facet, title = NULL, ...) {
  ggboxplot(dt, 
    x = x, 
    y = y, 
    facet.by = facet, 
    add = "jitter", 
    col = "Diagnosis_group",
    add.params = list(), 
    palette = get_palette("BuPu", 4)[-1],
    na.rm = TRUE,
    scales = "free_y", 
    nrow = 2,
    ...
  ) + 
  geom_pwc(
    aes(group = Diagnosis_group), 
    tip.length = 0,
    method = "tukey_hsd",
    bracket.nudge.y = -0.3
  ) +
  ggtitle(title) ### add a graying for close to y=x
}

# Main ------------------------------------------
# load overall base file
result <- fread("../../dataProt_SS-205063.merged.QC.csv")
a <- result[, 1:27]

# Map in PCs assumes your primary data table exists 
x <- fread("../../genomic/mypcair.PCs.csv", colClasses = c(PIDN = "character"))
a <- merge(a, x, by = "PIDN")

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
# get the first visit for each biomarker, with Diagnosis_group, Age, Gender, APOE, HybControlNormScale
b <- lapply(input, function(i) {
  x <- a[!is.na(get(i))][, value := get(i)]
  x <- x[, .(PIDN, DS1V, value, Diagnosis_group, Age, Gender, APOE, HybControlNormScale, CDRSUM, MOCATOTS)]
  # keep only first visit where plasma available
  x <- x[order(DS1V)][, head(.SD, 1), by = "PIDN"]
})
setattr(b, 'names', input)
b <- rbindlist(b, idcol = 'Biomarker')
# b[, .N, by = Biomarker]

# for this plot, reclassify MCI-AD as MCI
b[, Diagnosis_group := gsub("MCI-AD", "MCI", Diagnosis_group)]
b$Diagnosis_group <- factor(b$Diagnosis_group, levels = c("HC", "MCI", "AD"))

# plot vs Age (consider MOCA or memory composite) scatter w/trends ---------
# scatterplot of L2FCs
lapply(c("Age", "CDRSUM", "MOCATOTS"), function(i) {
  p1 <- plot_corr(b, 
    x = i, 
    y = "value", 
    facet = "Biomarker",
    group = "Diagnosis_group",
    xlab = paste(i, "at visit"),
    ylab = paste("Biomarker level")
  )
  ggsave(
    filename = file.path(output, paste(nameset, i, "reg.pdf", sep=".")), 
    plot = p1, 
    width = 10, 
    height = 7, 
    device = cairo_pdf
  )
  ggsave(
    filename = file.path(output, paste(nameset, i, "reg.png", sep=".")), 
    plot = p1, 
    width = 10, 
    height = 7, 
    dpi = 144, 
    type = cairo
  )
})

# plot vs Age (consider MOCA or memory composite) scatter w/ overall trend -----
# scatterplot of L2FCs
lapply(c("Age", "CDRSUM", "MOCATOTS"), function(i) {
  p1 <- plot_corr(b, 
    x = i, 
    y = "value", 
    facet = "Biomarker",
    group = "Diagnosis_group",
    xlab = paste(i, "at visit"),
    ylab = paste("Biomarker level"),
    trendline_type = "overall"
  )
  ggsave(
    filename = file.path(output, paste(nameset, i, "overall_reg.pdf", sep=".")), 
    plot = p1, 
    width = 10, 
    height = 7, 
    device = cairo_pdf
  )
  ggsave(
    filename = file.path(output, paste(nameset, i, "overall_reg.png", sep=".")), 
    plot = p1, 
    width = 10, 
    height = 7, 
    dpi = 144, 
    type = cairo
  )
})

######plot vs diagnosis boxplot
p1 = plot_box(b, 
  x = "Diagnosis_group", 
  y = "value", 
  facet = "Biomarker",
  xlab=paste("Diagnosis"),
  ylab=paste("Biomarker level")
)
ggsave(
  filename = file.path(output, paste(nameset, "Diag_box.pdf", sep=".")), 
  plot = p1, 
  width = 10, 
  height = 7, 
  device = cairo_pdf
)
ggsave(
  filename = file.path(output, paste(nameset, "Diag_box.png", sep=".")), 
  plot = p1, 
  width = 10, 
  height = 7, 
  dpi = 144, 
  type = cairo
)




#####

# p <- ggpairs(
  # data = dcast(result, NCBI_ID ~ dx, value.var = "t stat"),
  # upper = list(continuous = upper_fn),
  # diag = list(continuous = function(...) ggally_densityDiag(...) + theme_linedraw()),
  # lower = list(continuous = lower_fn),
  # columns = 2:9
# ) 

# ggsave(p,
  # filename = file.path(output, paste(nameset, "correlations.pdf", sep=".")), 
  # device = cairo_pdf,
  # width = 12, 
  # height = 12
# )
# ggsave(p,
  # filename = file.path(output, paste(nameset, "correlations.png", sep=".")), 
  # width = 12, 
  # height = 12,
  # units = "in", 
  # dpi = 144, 
  # type = "cairo"
# )
