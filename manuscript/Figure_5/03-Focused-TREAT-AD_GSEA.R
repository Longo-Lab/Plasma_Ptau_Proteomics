#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 03-Focused-TREAT-AD_GSEA.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Focused analysis of top GSEA BioDomain enrichments for specified modules
#'
#' Author: Robert B Butler III
#'
#' Date Created: 2026-06-17
#'
#' Copyright (c) 2026
#' Email: rrbutler@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'   Creates violin plots and summary tables for specified modules showing
#'   top BioDomain enrichments with labeled GO terms
#'
#' ---------------------------
#' load up the packages we will need:
library(optparse)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(patchwork)
library(data.table)
library(stringr)
library(dplyr)
library(gt)
library(RColorBrewer)
library(Cairo)
library(extrafont)

# read in args
option_list <- list(
  make_option(
    c("-m", "--modules"),
    type = "character",
    help = "Comma-separated list of modules to analyze (e.g., 'Merged.M8,Merged.M1')"
  ),
  make_option(
    c("--inpath"),
    type = "character",
    default = "../../WSS_prot",
    help = "Path to input datasets folder [default: %default]"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = ".",
    help = "File path to output directory [default: %default]"
  ),
  make_option(
    c("-r", "--ref"),
    type = "character",
    default = "../../module_scores/STRING.INT.WSS.terms.Zsum.csv",
    help = "File path to module term names [default: %default]"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    default = "dataProt_SS-205063.baseline_WSS",
    help = "filename prefix of input project [default: %default]"
  ),
  make_option(
    c("-s", "--suffix"),
    type = "character",
    default = "all_WSS",
    help = "output suffix to add to nameset [default: %default]"
  ),
  make_option(
    c("--n_top"),
    type = "integer",
    default = 15,
    help = "Number of top terms to show in table and label in plot [default: %default]"
  ),
  make_option(
    c("--exclude_biodomains"),
    type = "character",
    default = "none",
    help = "Comma-separated biodomains to exclude [default: %default]"
  )
)

args_list <- parse_args(OptionParser(option_list = option_list))

# read in fields
fields <- c(
  "modules",
  "inpath",
  "output",
  "ref",
  "nameset",
  "suffix",
  "n_top",
  "exclude_biodomains"
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

# Parse input variables
modules <- strsplit(modules, ",")[[1]]
exclude_biodomains <- strsplit(exclude_biodomains, ",")[[1]]

# Main ------------------------------------------
# Read all module GSEA results
all_modules <- modules
fh <- file.path(inpath, all_modules,
  paste(nameset, all_modules, "GSEA_Biodomain.rds", sep = ".")
)

# Check if files exist
if (!all(file.exists(fh))) {
  missing <- fh[!file.exists(fh)]
  stop("Missing GSEA result files:\n", paste(missing, collapse = "\n"))
}

result <- lapply(fh, readRDS)
setattr(result, 'names', all_modules)

# Extract relevant columns from each gseaResult object and convert to data.table
result_dt <- lapply(names(result), function(m) {
  
  # Extract the result data.frame from the gseaResult object
  df <- result[[m]]@result
  
  # Convert to data.table and select relevant columns
  dt <- as.data.table(df)[, .(
    Module = m,
    BioDomain = label,
    abbr = abbr,
    color = color,
    ONTOLOGY = ONTOLOGY,
    ID = ID,
    Description = Description,
    NES = NES,
    padj = p.adjust
  )]
  
  return(dt)
})
result_dt <- rbindlist(result_dt)

# Load term names for module annotations
names_dt <- fread(ref)
names_dt[, group := tstrsplit(name, "\\.", keep = c(1))]
names_dt[, group := { 
  group %>%
    str_remove_all(., '_|Plasma') %>%
    str_replace_all(., 'AB140P', 'A\u03b240') %>%
    str_replace_all(., 'AB142P', 'A\u03b242') %>%
    str_replace_all(., 'PTau', 'pTau') %>%
    str_replace_all(., 'ABRatio', 'A\u03b242/A\u03b240') %>%
    str_replace_all(., 'AB42Ratio', '/A\u03b242')
}]

# Process each module
for (m in modules) {
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Processing module:", m, "\n")
  cat(rep("=", 60), "\n\n", sep = "")
  
  # Get module annotation
  module_term <- names_dt[name == m, term]
  if (length(module_term) == 0) {
    module_term <- m
  }
  
  # Driver domains/proteins table ---------------------------------------------
  top_terms_table <- unique(
    result_dt[Module == m 
      & !BioDomain %in% exclude_biodomains
      & padj < 0.05
    ][order(-NES)]
  )[1:n_top, .(ID, Description, BioDomain, NES)]
  
  if (nrow(top_terms_table) > 0) {
    p_table <- gt(top_terms_table) |> 
      tab_header(title = md(paste0("**", m, " - ", module_term, "**"))) |>
      cols_align(
        align = "right",
        columns = c(BioDomain, NES)
      ) |>
      fmt_number(columns = NES, decimals = 3)
    
    ggsave(
      plot = wrap_table(p_table), 
      file = file.path(output, 
        paste(nameset, suffix, m, "BioDomain_Top_Terms.png", sep = ".")
      ),
      width = 10,
      height = max(4, nrow(top_terms_table) * 0.3)
    )
    
    cat("Saved top terms table for", m, "\n")
  } else {
    cat("No significant terms found for", m, "\n")
    next
  }
  
  # Focused plot for top biodomains -------------------------------------------
  
  # Get top biodomains for this module
  top_biodomains <- unique(
    result_dt[Module == m 
      & !BioDomain %in% exclude_biodomains
      & padj < 0.05
      & NES > 0
    ][order(-NES)][1:n_top, BioDomain]
  )
  
  if (length(top_biodomains) == 0) {
    cat("No upregulated biodomains found for", m, "\n")
    next
  }
  
  # Filter data for plotting
  plot_data <- unique(
    result_dt[
      Module == m 
      & BioDomain %in% top_biodomains
      & NES > 0
    ]
  )
  
  if (nrow(plot_data) == 0) {
    cat("No data to plot for", m, "\n")
    next
  }
  
  # Get top terms to label
  top_terms <- unique(
    result_dt[Module == m 
      & !BioDomain %in% exclude_biodomains
      & padj < 0.05
      & NES > 0
    ][order(-NES)]
  )[1:n_top, ID]
  
  # Add label column for top terms only
  plot_data[, label_text := ifelse(ID %in% top_terms, Description, "")]
  
  # # Check for duplicate descriptions within the labeled set
  # dup_descriptions <- plot_data[label_text != "", .N, by = label_text][N > 1, label_text]
  
  # # For duplicate descriptions, append the BioDomain abbreviation
  # plot_data[label_text %in% dup_descriptions & label_text != "", 
            # label_text := paste0(label_text, " (", abbr, ")")]
  
  # Reorder BioDomain factor by top NES
  plot_data[, BioDomain := factor(BioDomain, levels = rev(top_biodomains))]
  
  # Create the plot
  p_violin <- ggplot(plot_data, aes(NES, BioDomain)) + 
    theme_minimal() +
    scale_color_identity(guide = 'none') + 
    scale_fill_identity(guide = 'none') +
    scale_y_discrete(drop = FALSE) +
    scale_size_continuous(
      name = "-log10(FDR)",
      limits = range(-log10(plot_data$padj), na.rm = TRUE),
      range = c(0.5, 5)
    ) +
    # Non-significant points
    geom_jitter(
      data = subset(plot_data, padj > 0.05),
      color = 'grey85', size = 0.3, height = 0.2
    ) +
    # Violin plots for significant terms
    geom_violin(
      data = subset(plot_data, padj <= 0.05),
      scale = 'width', 
      aes(color = color), 
      lwd = 0.5,
      fill = NA
    ) +
    # Significant points
    geom_jitter(
      data = subset(plot_data, padj <= 0.05),
      aes(fill = color, color = color, size = -log10(padj)),
      color = 'grey20', 
      alpha = 0.4, 
      shape = 21,
      height = 0.2
    ) +
    # Label top terms
    geom_text_repel(
      data = subset(plot_data, label_text != ""),
      aes(label = label_text, color = color),
      size = 2.5,
      max.iter = 50000,
      max.overlaps = Inf,
      segment.size = 0.2,
      fontface = "bold",
      family = "Arial",
      box.padding = 0.5,
      point.padding = 0.3
    ) +
    labs(
      y = '',
      x = 'Normalized Enrichment Score (NES)',
      title = paste0(m, ' - ', module_term),
      subtitle = 'Upregulated BioDomain Enrichments (NES > 0, FDR < 0.05)'
    ) +
    theme(
      legend.position = 'bottom',
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12, family = "Arial"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, family = "Arial"),
      axis.text = element_text(size = 9, family = "Arial", face = "bold"),
      axis.title = element_text(size = 10, family = "Arial", face = "bold"),
      legend.text = element_text(size = 9, family = "Arial"),
      legend.title = element_text(size = 10, family = "Arial", face = "bold")
    ) +
    guides(size = guide_legend(nrow = 1))
  
  # Save the plot
  ggsave(
    filename = file.path(output, 
      paste(nameset, suffix, m, "Top_BioDomains_NES.pdf", sep = ".")
    ),
    plot = p_violin,
    width = 6,
    height = max(7, length(unique(plot_data$BioDomain)) * 0.6 + 2),
    device = cairo_pdf
  )
  
  cat("Saved violin plot for", m, "\n")
  
  # Create a summary table of the plotted data
  summary_table <- plot_data[padj <= 0.05, .(
    N_Terms = .N,
    Median_NES = median(NES),
    Mean_NES = mean(NES),
    Max_NES = max(NES),
    Min_padj = min(padj)
  ), by = BioDomain][order(-Median_NES)]
  
  fwrite(summary_table,
    file = file.path(output,
      paste(nameset, suffix, m, "BioDomain_Summary.csv", sep = ".")
    )
  )
  
  cat("Saved summary table for", m, "\n\n")
}

cat("\n", rep("=", 60), "\n", sep = "")
cat("Analysis complete!\n")
cat(rep("=", 60), "\n", sep = "")

