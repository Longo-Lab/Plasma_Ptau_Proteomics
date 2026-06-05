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

# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "List of input markers comma separated"
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
#' scatterplot with a trendline and stats, can facet
#'
#' @param dt data.table 
#' @param x string, column name to plot on x
#' @param y string, column name to plot on y
#' @param facet string, column name for faceting
#' @param title main title
#' ...: additional objects passed to ggscatter
#'
#' returns: ggobject
plot_corr = function (dt, x, y, facet, title = NULL, ...) {
  ggscatter(dt, 
    x = x, 
    y = y, 
    facet.by = facet, 
    add = "reg.line", 
    conf.int = TRUE, 
    col = "Diagnosis_group",
    add.params = list(), 
    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(
      method = "kendall", 
      label.x.npc = "left", 
      label.y.npc = "top"
      # label.sep = "\n"
    ),
    palette = get_palette("BuPu", 4)[-1],
    na.rm = TRUE,
    scales = "free_y", 
    nrow = 2,
    ...
  ) + 
  ggtitle(title) ### add a graying for close to y=x
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
a <- result
# a <- result[, 1:27]

# Map in PCs assumes your primary data table exists 
x <- fread("../../genomic/mypcair.PCs.csv", colClasses = c(PIDN = "character"))
a <- merge(a, x, by = "PIDN")

# drop outlier samples
drop = a[abs(ConnectivityZscore) > 10, SampleID]
a = a[!SampleID %in% drop]

# only APOE
a <- a[APOE != ""]

# Setup for baseline visit only analysis ---------------------------
# get the first visit for each biomarker, with Diagnosis_group, Age, Gender, APOE, HybControlNormScale
b <- lapply(input, function(i) {
  x <- a[!is.na(get(i))][, value := get(i)]
  x[, value := log2(value + 1)]
  x <- x[, .(PIDN, DS1V, value, Diagnosis_group, Age, Gender, APOE, HybControlNormScale)]
  # keep only first visit where plasma available
  x <- x[order(DS1V)][, head(.SD, 1), by = "PIDN"]
})
setattr(b, 'names', input)
b <- rbindlist(b, idcol = 'Probe')
# b[, .N, by = Biomarker]

# for this plot, reclassify MCI-AD as MCI
b[, Diagnosis_group := gsub("MCI-AD", "MCI", Diagnosis_group)]
b$Diagnosis_group <- factor(b$Diagnosis_group, levels = c("HC", "MCI", "AD"))

# plot vs Age (consider MOCA or memory composite) scatter w/trends ---------
# scatterplot of L2FCs
p1 = plot_corr(b, 
  x = "Age", 
  y = "value", 
  facet = "Probe",
  xlab=paste("Age at visit"),
  ylab=paste("Probe Log2RFU")
)
ggsave(
  filename = file.path(output, paste(nameset, "Age_reg.pdf", sep=".")), 
  plot = p1, 
  width = 10, 
  height = 7, 
  device = cairo_pdf
)
ggsave(
  filename = file.path(output, paste(nameset, "Age_reg.png", sep=".")), 
  plot = p1, 
  width = 10, 
  height = 7, 
  dpi = 144, 
  type = cairo
)


######plot vs diagnosis boxplot
p1 = plot_box(b, 
  x = "Diagnosis_group", 
  y = "value", 
  facet = "Probe",
  xlab=paste("Diagnosis"),
  ylab=paste("Probe Log2RFU")
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

######plot vs APOE boxplot
p1 = plot_box(b, 
  x = "APOE", 
  y = "value", 
  facet = "Probe",
  xlab=paste("APOE status"),
  ylab=paste("Probe Log2RFU")
)
ggsave(
  filename = file.path(output, paste(nameset, "APOE_box.pdf", sep=".")), 
  plot = p1, 
  width = 10, 
  height = 7, 
  device = cairo_pdf
)
ggsave(
  filename = file.path(output, paste(nameset, "APOE_box.png", sep=".")), 
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
