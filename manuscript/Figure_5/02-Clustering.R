#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 02-Clustering.R
#'
#' Version: 0.0.3
#'
#' Purpose of script: heatmaps and clustering of individuals by scores
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
#'   v0.0.3 - adding in eigengene flexibility
#'
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
library(ComplexHeatmap)
library(stringr)
library(NbClust)
# library(pheatmap)

# read in args
option_list <- list(
  make_option( 
    c("-i", "--input"),
    type = "character",
    help = "List of input score columns comma separated, if multiple scoretypes, separate with |"
  ),
  make_option(
    c("-s", "--score"),
    type = "character",
    help = "Scoretype to analyse, multiple comma separated"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project [default: Nameset of input]"
  ),
  make_option(
    c("-p", "--dat-prefix"),
    dest = "dat_prefix",
    type = "character",
    default = "dataProt_SS-205063",
    help = "Prefix of project dataset for filenames [default: Nameset of input]"
  ),
  make_option(
    c("-k", "--kmeans"),
    type = "integer",
    default = NULL,
    help = "Apply kmeans clustering (NbClust) either k-value or 0 for best [default: %default]"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = ".",
    help = "File path to output directory [default: %default]"
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
  "dat_prefix",
  "output",
  "data_path",
  "kmeans"
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

# split input variables and add suffixes if multiple datasets
score <- strsplit(score, ",")[[1]]
if (length(score) > 1) {
  # split scores and values into list of lists
  input <- strsplit(input, "\\|")[[1]]
  input <- lapply(input, function(i) strsplit(i, ",")[[1]])
  # define breaks in heatmap by rep score names for lengths of lists
  breaks <- sapply(input, length) 
  breaks <- unlist(
    mapply(
      FUN = function(x, y) rep(x, y), 
      x = str_extract(score, "^[^.]+"), 
      y = breaks, 
      SIMPLIFY = FALSE, 
      USE.NAMES = FALSE
    )
  )
  # collapse for rownames
  input <- unlist(
    lapply(seq(input), function(i) paste(input[[i]], score[i], sep = "."))
  )
} else {
  input <- strsplit(input, ",")[[1]]
  breaks = NULL
}


# Functions ---------------------------
#' Merge Two Data Tables with Custom Column Suffixes
#'
#' This function merges two data tables while adding custom suffixes to column names
#' to avoid naming conflicts. Suffixes are added before merging and then removed from
#' the original data tables after merging.
#'
#' @param x A data.table, the first table to merge
#' @param y A data.table, the second table to merge
#' @param by Character vector of column names to merge on
#' @param suffixes A list with two elements specifying the suffixes for x and y columns
#' @param ... Additional arguments passed to merge()
#' 
#' @return A merged data.table with columns from both input tables
#' @export
#'
#' @references \url{https://stackoverflow.com/a/19143137/8034603}
#'
#' @examples
#' require(data.table)
#' dt1 <- data.table(id = 1:3, value = c("a", "b", "c"))
#' dt2 <- data.table(id = 2:4, value = c("x", "y", "z"))
#' merged <- mergeWithSuffix(dt1, dt2, by = "id", suffixes = list("dt1", "dt2"))
mergeWithSuffix <- function(x, y, by, suffixes = NULL, ...) {

  # Add Suffixes
  .mkSuffix(x, suffixes[[1]], merge.col = by)
  .mkSuffix(y, suffixes[[2]], merge.col = by)

  # Merge
  ret <- merge(x, y, by = by, suffixes = NULL, ...)

  # Remove Suffixes
  .undoSuffix(x, suffixes[[1]], merge.col = by)
  .undoSuffix(y, suffixes[[2]], merge.col = by)
  return(ret)
}

#' Add Suffix to Column Names
#'
#' Helper function that adds a suffix to column names in a data.table,
#' excluding the specified merge columns.
#'
#' @param x A data.table whose column names will be modified
#' @param sfx The suffix to add to column names
#' @param sep The separator between column name and suffix, default is "."
#' @param merge.col Character vector of column names to exclude from renaming
#'
#' @return None, modifies the input data.table by reference
#' @keywords internal
.mkSuffix <- function(x, sfx, sep = ".", merge.col=NULL)  {
  if (sfx != "") {
    nms <- setdiff(names(x), merge.col)
    setnames(x, nms, paste(nms, sfx, sep = "."))
  }
}

#' Remove Suffix from Column Names
#'
#' Helper function that removes a suffix from column names in a data.table,
#' excluding the specified merge columns.
#'
#' @param x A data.table whose column names will be modified
#' @param sfx The suffix to remove from column names
#' @param sep The separator between column name and suffix, default is "."
#' @param merge.col Character vector of column names to exclude from renaming
#'
#' @return None, modifies the input data.table by reference
#' @keywords internal
.undoSuffix <- function(x, sfx, sep = ".", merge.col = NULL) {
  if (sfx != "") {
    nms <- setdiff(names(x), merge.col)
    setnames(x, nms, sub(paste0(get("sep"), sfx, "$"), "", nms))
  }
}



# Main ------------------------------------------
# load overall base file
result <- fread(paste0("../../", dat_prefix, ".merged.QC.csv"))
a <- result[, 1:27]

# Map in PCs assumes your primary data table exists 
x <- fread("../../genomic/mypcair.PCs.csv", colClasses = c(PIDN = "character"))
a <- merge(a, x, by = "PIDN")

# eigen_import <- function(dat_prefix, score, tail) {
  
# }

# Map in Scores assumes your primary data table exists 
tail <- "normalized_scores.csv"
score_val <- "NES"
scor <- NULL
if (length(score) > 1) {
  lapply(score, function(i){
    x <- fread(
      file = file.path(data_path, paste(dat_prefix, i, tail, sep = ".")), 
      colClasses = c(SampleID = "character")
    )
    assign(
      x = "a", 
      value = mergeWithSuffix(a, x, by = "SampleID", suffixes = c("", i)), 
      envir = .GlobalEnv
    )
  })
} else {
  x <- fread(
    file = file.path(data_path, paste(dat_prefix, score, tail, sep = ".")), 
    colClasses = c(SampleID = "character")
  )
  a <- merge(a, x, by = "SampleID")
}

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
  x <- x[, .(PIDN, DS1V, value, Diagnosis_group, Age, Gender, APOE, HybControlNormScale, CDRSUM, MOCATOTS, Subarray, Study, PlateId)]
  # keep only first visit where plasma available
  x <- x[order(DS1V)][, head(.SD, 1), by = "PIDN"]
})
setattr(b, 'names', input)
b <- rbindlist(b, idcol = 'Score')
# b[, .N, by = Score]
if (length(score) > 1) {
  b[, Score := factor(Score, levels = input)]
} else {
  b[, Score := factor(Score, levels = sort(unique(Score)))]
}
setnames(b, "value", score_val)

# for this plot, reclassify MCI-AD as MCI
b[, Diagnosis_group := gsub("MCI-AD", "MCI", Diagnosis_group)]
b$Diagnosis_group <- factor(b$Diagnosis_group, levels = c("HC", "MCI", "AD"))

# irritating annotation frames
annots <- data.frame(
  unique(b[, .(PIDN, Subarray, Study, PlateId, Gender, APOE, Age, CDRSUM, MOCATOTS, Diagnosis_group)]),
  row.names = 'PIDN'
)

# keep just the complete.cases (drop cols)
mat <- as.matrix(dcast(b, Score ~ PIDN, value.var = score_val), rownames = 'Score')
mat <- mat[, colSums(is.na(mat)) == 0]
mat <- t(scale(t(mat)))  # Transpose twice to scale rows instead of columns
mat <- mat[, match(rownames(annots), colnames(mat))]


# dendrogram cutoff --------------------------------------------
# Nbclust
if (kmeans == 0) {
  png(
    file = file.path(output, 
      paste(nameset, ifelse(length(score) > 1, "all", score), "heat_k.png", sep=".")
    ), 
    type = "cairo"
  )
  res <- NbClust(t(mat), method = "kmeans")
  dev.off()
  # Get the best number of clusters according to the majority rule
  kval <- as.numeric(names(sort(table(res$Best.nc[1,]), decreasing = TRUE)[1]))
  print(paste("Optimal number of column clusters:", kval))
} else if (kmeans > 0) {
  kval = as.integer(kmeans)
} else {
  kval = NULL
}


# heatmap w/ dendrogram ----------------------------------------
set.seed(3)
# Convert your pheatmap code to ComplexHeatmap
hm <- Heatmap(
  mat,
  name = "Score clustering across cohort",
  cluster_columns = TRUE,
  column_dend_height = unit(4, "cm"), 
  column_km = kval,
  cluster_rows = FALSE,
  show_column_names = FALSE,
  row_split = breaks,
  row_labels = str_extract(rownames(mat), "^[^.]+"),
  top_annotation = HeatmapAnnotation(
    df = annots,
    col = list(
      Diagnosis_group = setNames(get_palette("BuPu", 4)[-1], levels(b$Diagnosis_group))
    )
  ),
  # Set legend columns parameter
  heatmap_legend_param = list(
    title = paste("Scaled", score_val), 
    legend_width = unit(4, "cm"),
    nrow = 2  # Number of rows in the legend (creates columns)
  )
)
# Draw the heatmap and return as a grid object
drawn <- draw(hm)

set.seed(3)
ggsave(
  filename = file.path(output, 
    paste(nameset, ifelse(length(score) > 1, "all", score), "heat_clust.pdf", sep=".")
  ), 
  plot = grid.grabExpr(draw(hm)), 
  width = 18, 
  height = 4 + 2 * length(score), 
  device = cairo_pdf
)
set.seed(3)
ggsave(
  filename = file.path(output, 
    paste(nameset, ifelse(length(score) > 1, "all", score), "heat_clust.png", sep=".")
  ), 
  plot = grid.grabExpr(draw(hm)), 
  width = 18, 
  height = 4 + 2 * length(score), 
  dpi = 144, 
  type = cairo
)


# Save the clustering groups
fwrite(
  x = data.table(stack(column_order(drawn))),
  file = file.path(output, 
    paste(nameset, ifelse(length(score) > 1, "all", score), "heat_k.csv", sep=".")
  )
)


# z-summary score for representation across datasets? ----------


