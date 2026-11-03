#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 04-WSS_validation.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Compare Weighted Sum of Scores for STRING network modules
#'  
#'
#' Author: Robert B Butler III
#'
#' Date Created: 2025-11-12
#'
#' Copyright (c) 2025
#' Email: rrbutler@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'   v0.0.1 - initial build
#'   
#'   
#'   
#'
#' ---------------------------

#' load up the packages we will need:  (uncomment as required)
library(data.table)
library(optparse)
library(WGCNA)
library(RNOmni)
library(ggplot2)
library(ggpubr)
library(patchwork)


# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "celltype_scores",
    help = "Data path for eigengenes"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project [default: Nameset of input]"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "validation",
    help = "File path to output directory [default: %default]"
  ),
  make_option(
    c("-r", "--ref"),
    type = "character",
    default = "dataProt_SS-205063",
    help = "Reference expression data prefix to be validated [default: %default]"
   ),
  make_option(
    c("-v", "--val"),
    type = "character",
    default = "GNPC_table",
    help = "Validation dataset prefix[default: %default]"
   ),
  make_option(
    c("-p", "--prots"),
    type = "character",
    default = "dataProt_SS-205063.probe_list.csv",
    help = "Path to probe/GENE table for SOMA [default: %default]"
  ),
  make_option(
    c("-t", "--threads"),
    default = 0,
    help = "num threads [default: parallelly::availableCores()]"
  )
)

args_list <- parse_args(OptionParser(option_list = option_list))

# Nameset default case
if(is.null(args_list$nameset)){
  args_list$nameset <- nameset_default(args_list$input)
}

# read in fields
fields <- c(
  "input",
  "nameset",
  "output",
  "threads",
  "val",
  "ref"
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

# # testing vals
# input <- "celltype_scores"
# nameset <- "ADRCvGNPC.01"
# output <- "validation"
# threads <- 2
# val <- "GNPC_01"
# ref <- "dataProt_SS-205063"
# prots <- "dataProt_SS-205063.probe_list.csv"

# check the current active plan
if (threads == 0) {
  threads <- as.integer(parallelly::availableCores()) - 2
} else if (is.numeric(threads) & threads != 0) {
  threads <- as.integer(threads)
} else {
  cat("threads must be an integer value")
  stop(parse_args(OptionParser(option_list = option_list), args = "--help"),
    call. = FALSE
  )
}
print(sprintf("Setting up %d threads for samples", threads))

# check outpath status
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}


# Functions ------------------------------------------------------------------
# Much faster alternative using eigengenes and network statistics
quick_preservation_check <- function(multiExpr, colorh, nPermutations = 1000) {
  discovery_data <- multiExpr[[1]]$data
  validation_data <- multiExpr[[2]]$data
  
  if (nrow(discovery_data) != nrow(validation_data)) {
    min_samp <- min(nrow(discovery_data), nrow(validation_data))
    discovery_data <- discovery_data[1:min_samp, , drop = FALSE]
    validation_data <- validation_data[1:min_samp, , drop = FALSE]
  }
  modules <- unique(colorh)
  modules <- modules[modules != "MEgrey"]

  results <- data.frame(
    Module = character(),
    Size = numeric(),
    Eigengene_Cor = numeric(),
    Eigengene_Zscore = numeric(),
    Density_Cor = numeric(),
    Connectivity_Cor = numeric(),
    stringsAsFactors = FALSE
  )

  for(module in modules) {
    cat("Processing module:", module, "\n")

    module_genes <- names(colorh)[colorh == module]

    if(length(module_genes) < 5) next

    # Extract module data
    disc_module <- discovery_data[, module_genes]
    val_module <- validation_data[, module_genes]

    # 1. Eigengene preservation
    if(ncol(disc_module) > 1) {
      disc_eigen <- prcomp(disc_module)$x[,1]
      val_eigen <- prcomp(val_module)$x[,1]
      eigen_cor <- cor(disc_eigen, val_eigen, use = "complete.obs")

      # Permutation test for eigengene
      perm_cors <- replicate(nPermutations, {
        perm_idx <- sample(nrow(val_module))
        cor(disc_eigen, prcomp(val_module[perm_idx, ])$x[,1], use = "complete.obs")
      })
      
      eigen_zscore <- (eigen_cor - mean(perm_cors, na.rm = TRUE)) / sd(perm_cors, na.rm = TRUE)
    } else {
      eigen_cor <- cor(disc_module[,1], val_module[,1], use = "complete.obs")
      eigen_zscore <- NA
    }

    # 2. Density preservation (correlation of correlation matrices)
    if(ncol(disc_module) > 2) {
      disc_cor_mat <- cor(disc_module, use = "complete.obs")
      val_cor_mat <- cor(val_module, use = "complete.obs")

      # Upper triangle correlations
      upper_tri <- upper.tri(disc_cor_mat)
      density_cor <- cor(disc_cor_mat[upper_tri], val_cor_mat[upper_tri], use = "complete.obs")

      # 3. Connectivity preservation
      disc_connectivity <- rowSums(abs(disc_cor_mat)) - 1  # -1 to exclude self
      val_connectivity <- rowSums(abs(val_cor_mat)) - 1
      connectivity_cor <- cor(disc_connectivity, val_connectivity, use = "complete.obs")
    } else {
      density_cor <- NA
      connectivity_cor <- NA
    }

    # Store results
    results <- rbind(results, data.frame(
      Module = module,
      Size = length(module_genes),
      Eigengene_Cor = eigen_cor,
      Eigengene_Zscore = eigen_zscore,
      Density_Cor = density_cor,
      Connectivity_Cor = connectivity_cor
    ))
  }
  return(results)
}




# Main -----------------------------------------------------------------------
# Prep discovery cohort -------------
disc <- fread(file.path(paste(ref, "merged.QC.csv", sep = ".")))

# split string and explode piped terms into separate rows
prot <- fread(file.path(paste(ref, "probe_list.csv", sep = ".")))

# subset baseline samples
samp <- fread(file.path(paste(ref, "baseline_pidn.csv", sep = ".")),
  colClasses = list(character = c("SampleID", "PIDN"))
)
disc <- merge(samp, disc, by = c("SampleID", "PIDN"))
pheno <- disc[, .(PIDN, SampleID, Gender, Age, Diagnosis_group)] 

# Post filter simple table
cols <- c(2, 28:ncol(disc))
b <- merge(
  x = prot[, .(ix, GENE)], 
  y = transpose(disc[, ..cols], make.names = "PIDN", keep.names = "ix"), 
  all.x = TRUE
)
disc_mat <- as.matrix(
  x = transpose(b[, -1], make.names = "GENE", keep.names = "PIDN"), 
  rownames = "PIDN"
)

# Inverse rank normalization (INT) to z-scale
disc_mat <- apply(disc_mat, 2, RankNorm)


# Prep validation cohort -----------
vali <- fread(file.path(paste0(val, ".csv")))

# Convert -1 to pseudocount
vali[vali == -1] <- 0.01

# expand complexes & make matrix
samp <- names(vali)[3:ncol(vali)]
vali_mat <- as.matrix(
  x = transpose(vali[, -1], make.names = "GENE", keep.names = "SampleID"), 
  rownames = "SampleID"
)
vali_mat <- vali_mat[, colnames(disc_mat)]

# Inverse rank normalization (INT) to z-scale
vali_mat <- apply(vali_mat, 2, RankNorm)

# save some mem on GNPC's tiny instances
remove(vali, b, prot, disc)
gc()


# Module GENE membership -------------------------
eig <- fread(file.path("cell_type_prot", "STRING_merged_scaled_scores.csv"))
# for genes in multiple modules, keep by top mean_alpha_specificity
eig <- eig[order(-mean_alpha_scaled)][, .SD[1], by = symbol]
eig[, module := paste(Biomarker, cluster, sep = ".")]

# filter singles/doubles
keep <- eig[, .N, by = module][N > 2, module]
eig <- eig[module %in% keep]

# make color list for WGCNA
terms <- unique(eig[, .(Biomarker, module, term)])
terms$color <- labels2colors(terms$module)
lookup <- setNames(terms$color, terms$module)
colorh <- setNames(rep("grey", ncol(vali_mat)), colnames(vali_mat))
genes_in_modules <- intersect(names(colorh), eig$symbol)
colorh[genes_in_modules] <- lookup[eig$module[match(genes_in_modules, eig$symbol)]]

# Flip negative mean_beta signed proteins for signed analysis--------------
neg_proteins <- eig[mean_beta < 0, .(symbol, module)]

# Create flipped matrices
disc_mat_flipped <- disc_mat
vali_mat_flipped <- vali_mat

for(module in unique(neg_proteins$module)) {
  genes_to_flip <- neg_proteins[module == module, symbol]
  genes_to_flip <- genes_to_flip[genes_to_flip %in% colnames(disc_mat)]
  
  disc_mat_flipped[, genes_to_flip] <- -disc_mat_flipped[, genes_to_flip]
  vali_mat_flipped[, genes_to_flip] <- -vali_mat_flipped[, genes_to_flip]
}

remove(eig)
gc()


# Module preservation stats -----------------------------------------------
enableWGCNAThreads()
# Assuming you have:
# - multiExpr: list with discovery and validation expression data
# - colorh: module colors from discovery dataset
# - WSS_discovery: your WSS scores from ADRC
# Calculate module preservation statistics
mp <- modulePreservation(
  multiData = list(
    Discovery = list(data = disc_mat_flipped), 
    Validation = list(data = vali_mat_flipped)
  ),
  multiColor = list(Discovery = colorh, Validation = colorh),
  referenceNetworks = 1,     # Discovery dataset as reference
  testNetworks = 2,          # Validation dataset
  networkType = "signed",
  greyName = "grey",
  randomSeed = 12345,
  parallelCalculation = TRUE,
  verbose = 3
)

# Extract Zsummary scores
pres_stats <- as.data.table(
  x = mp$preservation$Z$ref.Discovery$inColumnsAlsoPresentIn.Validation, 
  keep.rownames = 'color'
)
pres_stats <- merge(pres_stats, terms, all.x = TRUE)

# Create preservation summary
pres_stats[, Preservation := ifelse(Zsummary.pres > 10, "Strong", ifelse(Zsummary.pres > 2, "Moderate", "Weak"))]
fwrite(pres_stats, 
  file = file.path(output, paste(nameset, "preservation_stats", "csv", sep = "."))
)
saveRDS(mp, 
  file = file.path(output, paste(nameset, "preservation_stats", "rds", sep = "."))
)

# Visualize
library(ggplot2)
library(ggrepel)

p <- ggplot(pres_stats[!is.na(module)], 
    aes(x = moduleSize, y = Zsummary.pres, color = Biomarker, label = term)) +
  geom_point(size = 3) +
  geom_text_repel(max.overlaps = 15) +
  geom_hline(yintercept = c(2, 10), linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(title = "Module Preservation Statistics",
       subtitle = "Zsummary > 10: Strong preservation, 2-10: Moderate, <2: Weak") +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(
    # legend.position = 'none'
  )

ggsave(
  plot = p,
  filename = file.path(output,
    paste(nameset, "preservation_stats", "pdf", sep = ".")
  ),
  width = 9
)


# Eigengene correlation -------------------------------------------------------
# # Calculate and visualize
# bootstrap_results <- bootstrap_eigengene_correlations(disc_eig, vali_eig)
# print(bootstrap_results)
# significance_tests <- test_preservation_significance(bootstrap_results)
# print(significance_tests)
# robustness_metrics <- advanced_robustness_analysis(bootstrap_results)
# print(robustness_metrics)
# 
# # Create plots
# plots <- plot_bootstrap_results(bootstrap_results, robustness_metrics)
# print(plots$combined)






