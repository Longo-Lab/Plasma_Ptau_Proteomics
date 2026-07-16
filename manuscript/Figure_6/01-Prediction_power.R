#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 01-Prediction_power.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Pairwise GSEA correlations of TREAT-AD in proteomic associations 
#'
#' Author: Robert B Butler III, Michael Peter Brown
#'
#' Date Created: 2024-11-18
#'
#' Copyright (c) 2024
#' Email: rrbutler@stanford.edu, mpb@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'
#'    If effects are attenuated by factor k < 1
#'    You need MORE samples to detect the same associations
#'
#'    Power ∝ n × effect_size²
#'    So if effect reduced by k, need n increased by 1/k²
#'
#' ---------------------------
#' load up the packages we will need:  (uncomment as required)
library(optparse)
library(ggplot2)
library(ggpubr)
library(data.table)
library(stringr)
library(plyr)
# library(ComplexHeatmap)
# library(circlize)
# library(RColorBrewer)
# library(GeneOverlap)
library(tools)
library(extrafont)

# Nameset helper function
nameset_default <- function(input){
  file_path_sans_ext(basename(input))
}

# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "List of input datasets comma separated"
  ),
  make_option(
    c("--inpath"),
    type = "character",
    default = "../../WSS_predict",
    help = "Directory of prediction outputs"
  ),
  make_option(
    c("-s", "--suffix"),
    type = "character",
    default = "top10",
    help = "File path to output directory [default: %default]"
  ),
  make_option(
    c("-c", "--corrs"),
    type = "character",
    default = "../../prediction_models/STRING_WSS/all_model/model_info/joinet/INT/joinet_model_correlations_STRING_WSS_INT.csv",
    help = "File path to prediction stats [default: %default]"
  ),
  make_option(
    c("-r", "--ref"),
    type = "character",
    default = "../../module_scores/STRING.INT.WSS.terms.csv",
    help = "File path to reference terms [default: %default]"
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
    default = "dataProt_SS-205063.follow-up",
    help = "Basename of project [default: Nameset of input]"
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
  "inpath",
  "suffix",
  "corrs",
  "ref",
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

#' for a string of ens ids: explode, look up symbols, return symbols as string
#'
#' @param str string of comma separated ens_ids
#' @param result data.table with ens_ids in GENE col
#' @param col colname to grab (i.e. Symbol)
#'
#' @return concatenated string of lookup terms 
get_mcols <- function(str, result, col) {
  require("stringr")
  y <- result[[col]][match(unlist(str_split(str, ",")), result$GENE)]
  return(paste0(y, collapse = ","))
}

#' Calculate power comparison metrics (continuous analysis)
#'
#' @param res_actual data.table, results from actual module score
#' @param res_predicted data.table, results from predicted module score
#' @param protein_col Character, protein identifier column. Default: "Symbol"
#' @param alpha Numeric, significance threshold. Default: 0.05
#'
#' @return Named vector of comparison metrics
power_metrics_continuous <- function(res_actual, res_predicted, 
                                    protein_col = "Symbol",
                                    alpha = 0.05) {
  
  # Significant proteins
  sig_actual <- res_actual[padj < alpha, get(protein_col)]
  sig_predicted <- res_predicted[padj < alpha, get(protein_col)]
  overlap <- intersect(sig_actual, sig_predicted)
  
  # Merge
  merged <- merge(
    res_actual[, c(protein_col, "beta", "se", "pval"), with = FALSE],
    res_predicted[, c(protein_col, "beta", "se", "pval"), with = FALSE],
    by = protein_col,
    suffixes = c("_actual", "_predicted")
  )
  
  # Correlations
  cor_spearman <- cor(merged$beta_actual, merged$beta_predicted, 
                     use = "complete.obs", method = "spearman")
  cor_pearson <- cor(merged$beta_actual, merged$beta_predicted,
                    use = "complete.obs", method = "pearson")
  
  # Attenuation (same direction only)
  merged[, same_dir := sign(beta_actual) == sign(beta_predicted)]
  attenuation <- merged[same_dir == TRUE, 
                       median(abs(beta_predicted) / abs(beta_actual), na.rm = TRUE)]
  
  # Other metrics
  var_ratio <- sd(merged$beta_predicted, na.rm = TRUE) / 
               sd(merged$beta_actual, na.rm = TRUE)
  se_ratio <- median(merged$se_predicted / merged$se_actual, na.rm = TRUE)
  direction_concordance <- 100 * mean(merged$same_dir, na.rm = TRUE)
  
  # Sample size multiplier
  n_multiplier <- ifelse(attenuation > 0 & attenuation < 1,
                        1 / (attenuation^2), NA)
  
  # Return as named vector
  c(
    n_actual = res_actual$nobs[1],
    n_predicted = res_predicted$nobs[1],
    n_proteins_tested = nrow(merged),
    n_sig_actual = length(sig_actual),
    n_sig_predicted = length(sig_predicted),
    n_overlap = length(overlap),
    pct_overlap = 100 * length(overlap) / length(sig_actual),
    cor_spearman = cor_spearman,
    cor_pearson = cor_pearson,
    direction_concordance = direction_concordance,
    attenuation = attenuation,
    var_ratio = var_ratio,
    se_ratio = se_ratio,
    n_multiplier = n_multiplier
  )
}


#' Generate diagnostic plots for power comparison
#'
#' @param res_actual data.table, results from actual module score
#' @param res_predicted data.table, results from predicted module score
#' @param protein_col Character, protein identifier column. Default: "Symbol"
#' @param alpha Numeric, significance threshold. Default: 0.05
#' @param title Character, plot title prefix. Default: NULL
#' @param filename Character, file path to save plot. Default: NULL (display only)
#'
#' @return Invisible merged data.table
plot_power_comparison <- function(res_actual, res_predicted,
                                 protein_col = "Symbol",
                                 alpha = 0.05,
                                 title = NULL,
                                 filename = NULL) {
  
  # Merge
  merged <- merge(
    res_actual[, c(protein_col, "beta", "se", "pval", "padj"), with = FALSE],
    res_predicted[, c(protein_col, "beta", "se", "pval", "padj"), with = FALSE],
    by = protein_col,
    suffixes = c("_actual", "_predicted")
  )
  
  # Calculate metrics for annotations
  cor_val <- cor(merged$beta_actual, merged$beta_predicted, 
                use = "complete.obs", method = "spearman")
  
  merged[, same_dir := sign(beta_actual) == sign(beta_predicted)]
  attenuation <- merged[same_dir == TRUE, 
                       median(abs(beta_predicted) / abs(beta_actual), na.rm = TRUE)]
  
  se_ratio <- median(merged$se_predicted / merged$se_actual, na.rm = TRUE)
  
  # Set up output
  if (!is.null(filename)) {
    png(filename, width = 10, height = 10, units = "in", res = 300)
  }
  
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  
  # 1. Beta comparison
  plot(merged$beta_actual, merged$beta_predicted,
       xlab = "Beta (actual module score)",
       ylab = "Beta (predicted module score)",
       main = if (!is.null(title)) paste(title, "- Effect sizes") else "Effect sizes",
       pch = 16, col = rgb(0, 0, 0, 0.3), cex = 0.8)
  abline(0, 1, col = "red", lty = 2, lwd = 2)
  if (!is.na(attenuation) && attenuation < 10) {
    abline(0, attenuation, col = "blue", lwd = 2)
  }
  legend("topleft", 
         c(sprintf("r = %.2f", cor_val),
           sprintf("Atten = %.2f", attenuation)),
         bty = "n", cex = 0.9)
  
  # 2. SE comparison
  plot(merged$se_actual, merged$se_predicted,
       xlab = "SE (actual)", ylab = "SE (predicted)",
       main = if (!is.null(title)) paste(title, "- Standard errors") else "Standard errors",
       pch = 16, col = rgb(0, 0, 0, 0.3), cex = 0.8)
  abline(0, 1, col = "red", lty = 2, lwd = 2)
  legend("topleft", sprintf("Ratio = %.2f", se_ratio), bty = "n", cex = 0.9)
  
  # 3. P-value comparison
  merged[, sig_group := fcase(
    padj_actual < alpha & padj_predicted < alpha, "Both sig",
    padj_actual < alpha & padj_predicted >= alpha, "Actual only",
    padj_actual >= alpha & padj_predicted < alpha, "Predicted only",
    default = "Neither"
  )]
  
  colors <- c("Both sig" = "red", "Actual only" = "blue", 
              "Predicted only" = "green", "Neither" = "gray")
  
  plot(-log10(merged$pval_actual), -log10(merged$pval_predicted),
       xlab = "-log10(p) actual", ylab = "-log10(p) predicted",
       main = if (!is.null(title)) paste(title, "- P-values") else "P-values",
       pch = 16, col = colors[merged$sig_group], cex = 0.8)
  abline(0, 1, col = "black", lty = 2, lwd = 1)
  abline(h = -log10(alpha), v = -log10(alpha), col = "gray", lty = 3)
  legend("topleft", names(colors), col = colors, pch = 16, cex = 0.8, bty = "n")
  
  # 4. Bland-Altman
  merged[, mean_beta := (beta_actual + beta_predicted) / 2]
  merged[, diff_beta := beta_predicted - beta_actual]
  mean_diff <- mean(merged$diff_beta, na.rm = TRUE)
  sd_diff <- sd(merged$diff_beta, na.rm = TRUE)
  bias_pct <- (abs(mean_diff) / mean(abs(merged$beta_actual), na.rm = TRUE)) * 100

  # Calculate percentage within LOA
  loa_lower <- mean_diff - 1.96 * sd_diff
  loa_upper <- mean_diff + 1.96 * sd_diff
  merged[, within_loa := diff_beta >= loa_lower & diff_beta <= loa_upper]
  in_loa <- 100 * mean(merged$within_loa, na.rm = TRUE)
  
  plot(merged$mean_beta, merged$diff_beta,
       xlab = "Mean beta", ylab = "Difference (predicted - actual)",
       main = if (!is.null(title)) paste(title, "- Bland-Altman") else "Bland-Altman",
       pch = 16, col = rgb(0, 0, 0, 0.3), cex = 0.8)
  abline(h = 0, col = "red", lty = 2, lwd = 2)
  abline(h = mean_diff, col = "blue", lwd = 2)
  abline(h = loa_upper, col = "blue", lty = 3)
  abline(h = loa_lower, col = "blue", lty = 3)
  legend("topleft", 
    sprintf("Bias = %.3f\nPercent of mean: %.1f%%", mean_diff, bias_pct), 
    bty = "n", 
    cex = 0.9
  )
  legend("topright", 
    sprintf("Points in LOA = %.1f%%", in_loa), 
    bty = "n", 
    cex = 0.9
  )
  
  if (!is.null(filename)) {
    dev.off()
    cat(sprintf("Plot saved to: %s\n", filename))
  }
  
  invisible(merged)
}


#' Plot required sample sizes with multiplier labels (ggplot2 version)
#'
#' @param power_results data.table with columns: module, n_predicted, n_multiplier
#' @param filename Character, file path to save. Default: NULL
#' @param width Numeric, plot width in inches. Default: 8
#' @param height Numeric, plot height in inches. Default: 6
#'
#' @return ggplot object
plot_sample_size_requirements_gg <- function(power_results,
                                            filename = NULL,
                                            width = 8,
                                            height = 6) {
  require("ggplot2")
  
  # Calculate n_required
  dt <- copy(power_results)
  dt[, n_required := round(n_predicted * n_multiplier)]
  dt[, multiplier_label := paste0(round(n_multiplier), "x")]
  
  # Order by n_required
  # dt <- dt[order(n_required)]

  # Wrap module names
  dt[, term_wrap := str_wrap(gsub("\\|", " \\| ", term), width = 30)]
  dt[, term_wrap := factor(term_wrap, levels = rev(term_wrap))]
  
  # Create plot
  p <- ggplot(dt, aes(x = term_wrap, y = n_required)) +
    geom_col(fill = "#8C1515", width = 0.7) +
    geom_text(
      aes(label = multiplier_label),
      color = "white",
      fontface = "bold",
      size = 6,
      family = "Arial",
      vjust = 0.5,
      hjust = 1.1
    ) +
    coord_flip() +
    labs(
      x = "WSS Module",
      y = "Required sample size",
      title = "Sample size needed to match actual WSS module power"
    ) +
    theme_minimal() +
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Save if filename provided
  if (!is.null(filename)) {
    ggsave(filename, p, width = width, height = height)
    cat(sprintf("Plot saved to: %s\n", filename))
  }
  
  return(p)
}


#' Plot correlation circles for WSS modules
#'
#' @param power_results data.table with columns: module, Pearson_Correlation, Correlation_P_Value
#' @param filename Character, file path to save. Default: NULL
#' @param width Numeric, plot width in inches. Default: 5
#' @param height Numeric, plot height in inches. Default: 6
#'
#' @return ggplot object
plot_correlation_circles <- function(power_results,
                                    filename = NULL,
                                    width = 5,
                                    height = 6) {
  require("ggplot2")
  require("stringr")
  
  # Prepare data
  dt <- copy(power_results)
  
  # Wrap module names (same as barplot)
  dt[, term_wrap := str_wrap(gsub("\\|", " \\| ", term), width = 30)]
  dt[, term_wrap := factor(term_wrap, levels = rev(term_wrap))]
  
  # Format correlation and p-value labels
  dt[, cor_label := sprintf("%.2f", Pearson_Correlation)]
  dt[, pval_label := sprintf("%.1e", Correlation_P_Value)]
  
  # Circle size based on absolute correlation
  dt[, circle_size := abs(Pearson_Correlation)]
  
  # Create plot
  p <- ggplot(dt, aes(x = 1, y = term_wrap)) +
    # Add circles
    geom_point(
      aes(size = circle_size, color = Pearson_Correlation),
      shape = 16
      # alpha = 0.8
    ) +
    # Add correlation labels inside circles
    geom_text(
      aes(label = cor_label),
      size = 3,
      fontface = "bold",
      family = "Arial",
      color = "white"
    ) +
    # Add p-value labels on the right
    geom_text(
      aes(x = 1.35, label = pval_label),
      size = 3,
      family = "Arial",
      hjust = 0
    ) +
    # Add "P-Value" header
    annotate(
      "text",
      x = 1.35,
      y = length(levels(dt$term_wrap)) + 0.7,
      label = "P-Value",
      size = 4,
      fontface = "bold",
      family = "Arial",
      hjust = 0
    ) +
    scale_size_continuous(
      range = c(3, 12),  # Reduced from c(8, 20)
      limits = c(0, 1),
      guide = "none"
    ) +
    scale_color_gradient2(
      low = "#4B92DB",
      mid = "white",
      high = "#8C1515",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Pearson r"
    ) +
    scale_x_continuous(limits = c(0.7, 1.6), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(
      x = NULL,
      y = "WSS Module",
      title = "Actual vs Predicted Correlation"
    ) +
    theme_minimal() +
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, lineheight = 0.9),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 12, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      plot.margin = margin(5, 40, 5, 5)  # Extra right margin for p-values
    )
  
  # Save if filename provided
  if (!is.null(filename)) {
    ggsave(filename, p, width = width, height = height)
    cat(sprintf("Plot saved to: %s\n", filename))
  }
  
  return(p)
}



# Main ------------------------------------------
# DE file locations
fh <- file.path(inpath, input,
  paste(nameset, input, "csv", sep = ".")
)
result_act <- lapply(fh, fread)
setattr(result_act, 'names', input)

fp <- file.path(inpath, paste0(input, ".pred"),
  paste(nameset, input, "pred.csv", sep = ".")
)
result_pred <- lapply(fp, fread)
setattr(result_pred, 'names', input)

# make table of comparisons ----------------
# Calculate metrics for all
power_comparison_results <- lapply(input, function(ep) {
  
  # Get results for this module score
  # (assuming you have them stored somehow)
  res_actual <- result_act[[ep]]
  res_predicted <- result_pred[[ep]]
  
  # Calculate metrics
  metrics <- power_metrics_continuous(
    res_actual, 
    res_predicted,
    protein_col = "Symbol"
  )
  
  # calculate n necessary
  dt <- setDT(as.list(metrics))
  dt[, n_required := n_actual * n_multiplier]
  # add in required for median effect size.....
  return(dt)
})
setattr(power_comparison_results, 'names', input)
power_comparison_results <- rbindlist(power_comparison_results, idcol = 'module')

# add names
terms <- fread(ref)
power_comparison_results <- merge(
  x = power_comparison_results, 
  y = terms, 
  by.x = 'module',
  by.y = 'name',
  sort = FALSE
)

fwrite(power_comparison_results,
  file = file.path(output, paste(nameset, suffix, "power_comparison.csv", sep = "."))
)

# plot 
p <- plot_sample_size_requirements_gg(
  power_comparison_results,
  filename = file.path(output, 
    paste(nameset, suffix, "power_comparison.pdf", sep = ".")
  )
)

# individual pairwise plots ------------------------------------
# Generate plots for each
for (ep in input) {
  plot_power_comparison(
    res_actual = result_act[[ep]],
    res_predicted = result_pred[[ep]],
    title = ep,
    filename = file.path(output, 
      paste(nameset, suffix, ep, "power_comparison.png", sep = ".")
    )
  )
}


# Possibly some summary corrplot? ---------------------------------------
samp_corr <- fread(corrs)
power_corr <- merge(
  x = power_comparison_results, 
  y = samp_corr[variant == "Biomarkers_Age_Sex", .(
    outcome, Pearson_Correlation, Correlation_P_Value
  )],
  by.x = 'module',
  by.y = 'outcome',
  sort = FALSE
)

p_cor <- plot_correlation_circles(power_corr,
  filename = file.path(output, 
    paste(nameset, suffix, "correlation_comparison.pdf", sep = ".")
  )
)


# Add pairwise actual vs prediction heatmaps for BioDomain ----------------




