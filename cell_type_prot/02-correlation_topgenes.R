#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 02-correlation_topgenes.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Show top cell-type proteins with biomarker association
#'
#' Author: Robert B Butler III
#'
#' Date Created: 2025-03-18
#'
#' Copyright (c) 2025
#' Email: rrbutler@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'
#'   Chooses the top genes for each TREAT-AD biodomain overlap by the 
#'     abs(sum(DE L2FC & TREAT-AD meta-analysis Ratio of means))
#'
#'     
#'   Usage:
#'     sbatch -J MG --mem=5G -c 2 -t 01:00:00 -p interactive \
#'       -o %x/%A_biodomain_enrichment_%x.log \
#'       --wrap "ml R/4.0; biodomain_enrichment.R -n PS19_C31 -g PS19 -d C31 -r R6 -t MG"
#'
#'   or interactive session:
#'     sdev -m 5 -c 2 -t 01:00:00 -p interactive
#'
#' ---------------------------

#' load up the packages we will need:  (uncomment as required)

library(data.table)
library(corrplot)
library(stringr)
library(stringi)
library(optparse)
library(Cairo)
library(RColorBrewer)
library(ggcorrplot)
library(ggrepel)
library(ggh4x)
library(patchwork)





#' ----------------------------

# default variables

# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = ".",
    help = "File path to DE input file [default: %default]"
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
  ),
  make_option(
    c("-r", "--ref"),
    type = "character",
    default = "../cell_type_prot",
    help = "File path to directory of reference data [default: %default]"
  ),
  make_option(
    c("--rtype"),
    type = "character",
    default = "isolated",
    help = "Reference type isolated or cultured [default: %default]"
  ),
  make_option(
    c("--lab-n"),
    dest = "lab_n",
    default = 20,
    help = "Number of top corr genes to label [default: %default]"
  ),
  make_option(
    c("--thresh"),
    type = "character",
    default = "stat_padj",
    help = "Significance threshold ('t stat', 'stat_pval', 'stat_padj', 'stat_decile', 'stat_quartile') [default: %default]"
  ),
  make_option(
    c("-b", "--biomarker"),
    type = "character",
    help = "String specifying the biomarker to plot [required]"     
  )
)
args_list = parse_args(OptionParser(option_list = option_list))

# read in fields
fields = c(
  "input",
  "output",
  "nameset",
  "ref",
  "rtype",
  "lab_n",
  "thresh",
  "biomarker"
)
for (field in fields) {
  if (!field %in% names(args_list)) {
    cat(str_c('Argument "', field, '" must be provided\n'), fill=T)
    stop(parse_args(OptionParser(option_list=option_list), args='--help'), 
         call.=FALSE)
  }
  
  assign(field, args_list[[field]])
  print(paste0(field, ': ', get(field)))
}


# Prep DE data -------------------------------------------------------
# DE file locations
fh <- file.path(input, biomarker,
  paste("dataProt_SS-205063.PC.related", biomarker, "csv", sep = ".")
)
result <- fread(fh)

# # to get the value in the dx column if APOE
# bio_no_apoe <- gsub(".APOE", "", biomarker)

# thresholds of cutoffs, pval, padj, outer_decile, outer_quartile ### maybe unnecessary
result[, stat_pval := ifelse(pval < 0.05, `t stat`, NA)]
result[, stat_padj := ifelse(padj < 0.05, `t stat`, NA)]
qs <- quantile(result[, `t stat`], c(0.1, 0.25, 0.75, 0.9))
result[, stat_decile := ifelse(`t stat` < qs[[1]] | `t stat` > qs[[4]], `t stat`, NA)]
result[, stat_quartile := ifelse(`t stat` < qs[[2]] | `t stat` > qs[[3]], `t stat`, NA)]

# Prep Reference data -----------------------------------------------------
if (rtype == 'specificity') {
  trt <- fread(paste0(ref, "/Sharma_etal.EWCE.human.lvl1.spec_matr.txt"))
  setnames(trt, 
    old = c("GENE", "Neuron", "Astrocyte", "Oligodendrocyte"), 
    new = c("Symbol", "Neurons", "Astrocytes", "Oligodendrocytes")
  )
  trt <- melt(trt, variable.name = "celltype", value.name = "z.log2LFQ", id.vars = "Symbol")
  trt <- unique(merge(trt, result[, .(Symbol, GENE)]))
  trt[, Symbol := NULL]
} else {
  trt <- fread(paste0(ref, "/Sharma_etal_prot_celltypes.", rtype, ".scaled.csv"))
}

# merge tables
dt <- merge(trt, 
  result[, .(GENE, Symbol, `t stat`, stat_pval, stat_padj, stat_decile, stat_quartile)], 
  by = "GENE"
)



# Plot key genes per category --------------------------------------------
# checking the number of genes matched (not NA in both) per group to any mouse de
counts <- dt[
  !is.na(get(thresh)) & !is.na(z.log2LFQ), 
  .(.N, r = cor(z.log2LFQ, get(thresh), use = "complete.obs")), 
  by = celltype
][order(celltype)]

# grab only the top five in each section
dt[, comb := sum(abs(z.log2LFQ), abs(get(thresh))), by=.I]
dt[, comb_sign := ifelse(sign(get(thresh)) == sign(z.log2LFQ), 1, -1)]
# c <- dt[order(-comb)][head(.SD, 10), by=c("celltype", "comb_sign")]

# Plot
plist <- lapply(seq(nrow(counts)), function(i) {
  ggplot(dt[celltype == counts[i, celltype]], 
    aes(x = z.log2LFQ, y = get(thresh))) +
  geom_point(size = 4, col = "#2166ac") + 
  geom_text_repel(
    data = dt[celltype == counts[i, celltype] & comb_sign == sign(counts[i, r])
              ][order(-comb)][1:as.integer(lab_n)],
    aes(label = Symbol), 
    fontface = "bold",
    max.overlaps = 15
  ) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_minimal() +
  labs(
    x = ifelse(rtype == "specificity", "Specificity", "z-scaled log2LFQ"), 
    y = "Biomarker T stat"
  ) +
  ggtitle(counts[i, celltype]) +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(face = "bold", size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
})
# wrap_plots(plist) + plot_layout(nrow = 1, guides = "collect")

ggsave(
  plot = wrap_plots(plist) + 
    plot_layout(nrow = 1, guides = "collect") +
    plot_annotation(
      title = paste(biomarker, "vs", rtype, "cell type protein abundance"), 
      theme = theme(plot.title = element_text(size = 20, face = "bold")
    )
  ),
  filename = file.path(output, 
    paste(
      nameset, 
      biomarker,
      rtype, 
      thresh, 
      "correlations", 
      paste0("top", lab_n), 
      "pdf", 
      sep="."
    )
  ), 
  width = 16, 
  height = 5,
  units = "in", 
  device = cairo_pdf
)




# # cleanup
# system("squeue -j $SLURM_JOBID --Format=TimeUsed")
