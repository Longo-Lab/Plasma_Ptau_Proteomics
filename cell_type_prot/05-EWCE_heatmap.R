#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 04-EWCE_heatmap.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: replot EWCE associations as a heatmap
#'
#' Author: Robert B Butler III
#'
#' Date Created: 2026-04-27
#'
#' Copyright (c) 2026
#' Email: rrbutler@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'     
#'   Usage:
#'     sbatch -J MG --mem=5G -c 2 -t 01:00:00 -p interactive \
#'       -o %x/%A_biodomain_enrichment_%x.log \
#'       --wrap "ml R/4.3.3; 04-EWCE_heatmap.R -i bio1,bio2,bio3 -n all_bios -gset stat_padj -cset Sharma_etal.cultured"
#'
#'   or interactive session:
#'     sdev -m 5 -c 2 -t 01:00:00 -p interactive
#'
#' ---------------------------

#' load up the packages we will need:  (uncomment as required)

library(ggplot2)
library(dplyr)
library(tidyr)
library(optparse)
library(data.table)
library(extrafont)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(Cairo)
library(stringr)


#' ----------------------------

# default variables

# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Biomarker or group names list to analyze"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project "
  ),
  make_option(
    c("-c", "--cset"),
    type = "character",
    help = "cell type dataset nickname (with path)"
  ),
  make_option(
    c("-g", "--gset"),
    type = "character",
    default = "stat_padj",
    help = "geneset tested: [default %default]"
  ),
  make_option(
    c("-w", "--width"),
    dest = "wid",
    default = 5,
    help = "width of the plot, more for more cell types [default: %default]"
  ),
  make_option(
    c("--height"),
    dest = "hgt",
    default = 4,
    help = "height of the plot, more for more cell types [default: %default]"
  ),
  make_option(
    c("-l", "--level"),
    default = 1,
    help = "annotation level (1 or 2) [default: %default]"
  ),
  make_option(
    c("-r", "--ref"),
    default = "../biomarker_names.txt",
    help = "reference for biomarker names/groups [default: %default]"
  ),
  make_option(
    c("-o", "--outpath"),
    type = "character",
    default = ".",
    help = "File path to output directory [default: %default]"
  ),
  make_option(
    c("--inpath"),
    type = "character",
    default = ".",
    help = "File path to input data [default: %default]"
  )
)
args_list = parse_args(OptionParser(option_list = option_list))

# read in fields
fields <- c(
  "input",
  "nameset",
  "cset",
  "gset",
  "wid",
  "hgt",
  "level",
  "ref",
  "outpath",
  "inpath"
)
for (field in fields) {
  if (!field %in% names(args_list)) {
    cat(paste0('Argument "', field, '" must be provided\n'), fill = TRUE)
    stop(parse_args(OptionParser(option_list=option_list), args='--help'), 
         call.=FALSE)
  }
  
  assign(field, args_list[[field]])
  print(paste0(field, ': ', get(field)))
}

# make outpath
if(!dir.exists(outpath)) dir.create(outpath)

# input variables
input <- strsplit(input, ",")[[1]]



# Prep gene lists -------------------------------------------------------
# DE file locations
fh <- file.path(inpath, 
  paste(nameset, gset, cset, paste0("lvl", level), "rds", sep = ".")
)
print(sprintf("Reading in enrichment data from %s", fh))
result <- readRDS(fh)

# merge
heatmap_data <- lapply(result, function(i) data.table(i$result))
heatmap_data <- rbindlist(heatmap_data, idcol = 'Biomarker')

# issues with rounded off p-values at 4 digits, calculate from sd
heatmap_data[p == 0 , p:= 2 * pnorm(-abs(sd_from_mean))]
heatmap_data[, log10p := -log10(p)]

# Heatmap -------------------------------------------------------------
# load term names 
names <- fread(ref)
if (!ref %like% "biomarker_names.txt") {
  names[, group := tstrsplit(name, "\\.", keep = c(1))]
  names[, group := { 
    group %>%
      str_remove_all(., '_|Plasma') %>%
      str_replace_all(., 'AB140P', 'A\u03b240') %>%
      str_replace_all(., 'AB142P', 'A\u03b242') %>%
      str_replace_all(., 'PTau', 'pTau') %>%
      str_replace_all(., 'ABRatio', 'A\u03b242/A\u03b240') %>%
      str_replace_all(., 'AB42Ratio', '/A\u03b242')
  }]
  g_names <- c("Merged", sort(unique(names[group != "Merged", group])))
  names[, group := factor(group, levels = g_names)]
  names <- names[match(input, name)]
} else {
  names[, group := factor(group, levels = c("A\u03b2", "Tau", "Ratio", "NDegen"))]
  setorder(names, group, term)
}

# add to scores
scores <- merge(heatmap_data, names, by.x = 'Biomarker', by.y = 'name', all.x = TRUE)
# scores$BioDomain <- factor(scores$BioDomain, levels = b_names)
scores$term <- factor(scores$term, levels = unique(scores[order(group, term), term]))
setorder(scores, group, term)

# report for syncing up heatmaps
print(
  do.call(sprintf, 
    c("-log10p min: %.02f, max %.02f", as.list(range(scores$log10p)))
  )
) # maybe could make a combined heatmap option in the future


# Convert to matrix
score_mat <- as.matrix(
  x = dcast(scores, term ~ CellType, value.var = "log10p", fill = 0),
  rownames = "term"
)
sig_mat <- as.matrix(
  x = dcast(scores, term ~ CellType, value.var = "q", fill = 0),
  rownames = "term"
)

# Create color function
col_fun <- colorRamp2(
  c(0, 8), # manually set to sync heatmaps from datasets
  c("white", "#85070C")
)

# Create right-side row annotation with GO term counts
row_ha <- rowAnnotation(
  Biomarker = names$group[match(rownames(score_mat), names$term)],
  col = list(
    Biomarker = setNames(
      object = brewer.pal(length(levels(names$group)), "Dark2"), 
      nm = levels(names$group))
  ),
  show_legend = TRUE,
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial")
)

# Create heatmap
p <- Heatmap(
  score_mat,
  name = "-log10p",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10, fontface = 'bold', fontfamily = 'Arial'),
  row_names_side = 'left',
  column_names_gp = gpar(fontsize = 9, fontface = 'bold', fontfamily = 'Arial'),
  column_names_side = 'bottom',
  column_names_rot = 45,
  column_title = paste("Core Biomarker Cell Type Enrichment", cset),
  left_annotation = row_ha,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(sig_mat[i, j] < 0.05) {
      grid.rect(x = x, y = y, width = width, height = height, 
        gp = gpar(col = "black", fill = NA, lwd = 1.5)
      )
    }
  },
  heatmap_legend_param = list(
    # title = "Score\nmeanNES * prop(FDR < 0.05)",
    # direction = "horizontal",
    title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial"),
    labels_gp = gpar(fontsize = 9, fontfamily = "Arial")
  )
)

CairoPDF(
  file = file.path(outpath, 
    paste(nameset, gset, cset, paste0("lvl", level), "heatmap.pdf", sep = ".")
  ), 
  width = wid, 
  height = hgt
)
draw(p)
dev.off()



# # Create the heatmap with gene sets as rows, cell types as columns
# p <- ggplot(heatmap_data, aes(x = CellType, y = Biomarker, fill = -log10(p))) +
  # geom_tile(color = "grey90", linewidth = 0) +
  # # Add bold black borders for significant cells
  # geom_tile(data = filter(heatmap_data, q < 0.05),
            # aes(x = CellType, y = Biomarker),
            # color = "black", 
            # linewidth = 1,
            # fill = NA) +
  # scale_fill_distiller(
    # palette = "Purples",
    # direction = 1,
    # name = "-log10(p)",
    # na.value = "grey50",
    # limits = 
  # ) +
  # theme_classic() +
  # scale_y_discrete(limits = rev) +  # Flip y-axis
  # theme(
    # text = element_text(face = 'bold', family = "Arial"),
    # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    # axis.text.y = element_text(size = 10),
    # panel.grid = element_blank(),
    # axis.ticks = element_blank(), 
    # axis.line = element_blank() 
  # ) +
  # labs(
    # x = sprintf("%s Cell Type", cset), 
    # y = "Core Biomarker",
    # title = "EWCE Core Biomarker Cell Type Enrichment",
  # )
# # p

# # Save
# ggsave(p,
  # filename = file.path(outpath, 
    # paste(nameset, gset, cset, paste0("lvl", level), "heatmap.pdf", sep = ".")
  # ), 
  # width = wid, 
  # height = 4
# )

# # cleanup
# system("squeue -j $SLURM_JOBID --Format=TimeUsed")
