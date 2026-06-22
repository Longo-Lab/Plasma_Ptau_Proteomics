#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 01-sila_age_onset_estimation.R
##
## Version: 0.0.1
##
## Purpose of script: Prepare PlasmaPTau217/AB42 sample data for Age of onset of AD prediction
##
## Author: Robert B Butler III
##
## Date Created: 2026-03-03
##
## Copyright (c) 2026
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##   As noted in 01a-Examining_SomaLogic_QC.R, it is really important to use the
##   data that is scale to the external reference. When comparing a protein 
##   across samples, it is also important to use the inverse rank normalized
##   data, as that is then scaled approriately to actually compare the lmm betas
##
##   REMEMBER, WHEN COMMENTING IN CODE, DO NOT REFERENCE INDIVIDUAL INFORMATION:
##   SAMPLES/PIDNS, DATES OR ANY INDIVIDUAL LEVEL INFO
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(data.table)
library(silaR)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(dplyr)
library(tibble)
library(optparse)
library(stringr)
library(RColorBrewer)
library(Cairo)
library(emmeans)
library(extrafont)

# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "../dataProt_SS-205063.merged.QC.csv",
    help = "File path to DE input file [default: %default]"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = ".",
    help = "File path to output directory [default: %default]"
  ),
  make_option(
    c("-b", "--biomarker"),
    type = "character",
    default = "PlasmapTau217_AB42Ratio",
    help = "File path to output directory [default: %default]"
  ),
  make_option(
    c("-e", "--eigengenes"),
    type = "character",
    default = "../module_scores/dataProt_SS-205063.STRING.INT.WSS.csv",
    help = "File path to eigengenes csv [default: %default]"
  ),
  make_option(
    c("-m", "--pcs"),
    type = "character",
    default = "../genomic/mypcair.PCs.csv",
    help = "File path to PCs csv [default: %default]"
  ),
  make_option(
    c("-r", "--ref"),
    type = "character",
    default = "../module_scores/STRING.INT.WSS.terms.Zsum.csv",
    help = "File path to module term names [default: %default]"
  ),
  make_option(
    c("-f", "--outliers"),
    type = "character",
    default = "iqr",
    help = "Filter the dataset by Cook's distance (> 4/n) or z-score (|z| > 3) instead of IQR (* 1.5) c('iqr', 'cooksd', 'zscore')"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    default = "dataProt_SS-205063.SILA",
    help = "Name of project [default: Nameset of input]"
  )
)

args_list <- parse_args(OptionParser(option_list = option_list))


# read in fields
fields <- c(
  "input",
  "output",
  "biomarker",
  "eigengenes",
  "pcs",
  "ref",
  "outliers",
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


# set env variables for data.table
envs <- list(aa = biomarker)
pTau_thresh <- 0.006312

# Functions --------------------------------------------------------------------
#' Plot Longitudinal Biomarker Data with Outlier Thresholds
#'
#' Creates a spaghetti plot showing individual trajectories of a biomarker over
#' age, with horizontal reference lines indicating outlier detection thresholds
#' (1.5×IQR and 3×IQR) and a clinical threshold.
#'
#' @param data A data.table or data.frame containing the biomarker data
#' @param biomarker Character string specifying the column name of the biomarker
#'   to plot
#' @param age_var Character string specifying the column name for age variable.
#'   Default is "AGE"
#' @param id_var Character string specifying the column name for participant ID.
#'   Default is "PIDN"
#' @param color_var Character string specifying the column name for the grouping
#'   variable used for coloring lines. Default is "ever_over"
#' @param threshold Numeric value for the clinical threshold to display as a
#'   horizontal line. Default is NULL (no threshold line)
#' @param iqr_values Numeric vector of length 3 containing pre-computed c(Q1, Q3, IQR)
#'   values. If NULL (default), values are calculated from the data. Use this to
#'   maintain consistent reference lines across filtered datasets.
#' @param plot_title Character string for the plot title. If NULL (default),
#'   creates a title using the biomarker name
#' @param show_iqr_labels Logical indicating whether to show IQR threshold labels.
#'   Default is TRUE
#' @param point_size Numeric value for the size of points. Default is 2
#' @param line_alpha Numeric value between 0 and 1 for IQR line transparency.
#'   Default is 0.6
#'
#' @return A ggplot2 object
#'
#' @details
#' The function calculates Q1, Q3, and IQR from the biomarker data (or uses
#' pre-specified values) and plots:
#' \itemize{
#'   \item Individual trajectories colored by the grouping variable
#'   \item 1.5×IQR outlier thresholds (transparent red dashed lines)
#'   \item 3×IQR extreme outlier thresholds (more transparent red dashed lines)
#'   \item Optional clinical threshold (gray dashed line)
#'   \item Labeled reference lines on the right side of the plot
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage - calculate IQR from data
#' p1 <- plot_biomarker_outliers(
#'   data = my_data,
#'   biomarker = "PlasmapTau217_AB42Ratio",
#'   threshold = 0.025
#' )
#'
#' # Save IQR values for use with filtered data
#' Q1 <- quantile(my_data$PlasmapTau217_AB42Ratio, 0.25)
#' Q3 <- quantile(my_data$PlasmapTau217_AB42Ratio, 0.75)
#' iqr <- IQR(my_data$PlasmapTau217_AB42Ratio)
#' iqr_vals <- c(Q1, Q3, iqr)
#'
#' # Use same IQR values with filtered data
#' p2 <- plot_biomarker_outliers(
#'   data = filtered_data,
#'   biomarker = "PlasmapTau217_AB42Ratio",
#'   threshold = 0.025,
#'   iqr_values = iqr_vals,
#'   plot_title = "Filtered Data"
#' )
#'
#' # Custom variables and no labels
#' plot_biomarker_outliers(
#'   data = my_data,
#'   biomarker = "PlasmapTau217_AB42Ratio",
#'   age_var = "age_at_visit",
#'   id_var = "subject_id",
#'   color_var = "diagnosis_group",
#'   threshold = 0.025,
#'   show_iqr_labels = FALSE
#' )
#' }
#'
#' @import ggplot2
#' @importFrom stats quantile IQR
#' @export
plot_biomarker_outliers <- function(data,
                                     biomarker,
                                     age_var = "AGE",
                                     id_var = "PIDN",
                                     color_var = "ever_over",
                                     threshold = NULL,
                                     iqr_values = NULL,
                                     plot_title = NULL,
                                     show_iqr_labels = TRUE,
                                     point_size = 2,
                                     line_alpha = 0.6) {
  
  # Input validation
  if (!biomarker %in% names(data)) {
    stop(paste("Biomarker column", biomarker, "not found in data"))
  }
  if (!age_var %in% names(data)) {
    stop(paste("Age variable", age_var, "not found in data"))
  }
  if (!id_var %in% names(data)) {
    stop(paste("ID variable", id_var, "not found in data"))
  }
  if (!color_var %in% names(data)) {
    stop(paste("Color variable", color_var, "not found in data"))
  }
  
  # Calculate or use provided quartiles and IQR
  if (is.null(iqr_values)) {
    Q1 <- quantile(data[[biomarker]], 0.25, na.rm = TRUE)
    Q3 <- quantile(data[[biomarker]], 0.75, na.rm = TRUE)
    iqr <- IQR(data[[biomarker]], na.rm = TRUE)
  } else {
    # Validate iqr_values
    if (!is.numeric(iqr_values) || length(iqr_values) != 3) {
      stop("iqr_values must be a numeric vector of length 3: c(Q1, Q3, IQR)")
    }
    Q1 <- iqr_values[1]
    Q3 <- iqr_values[2]
    iqr <- iqr_values[3]
  }
  
  # Create plot title if not provided
  if (is.null(plot_title)) {
    plot_title <- paste(biomarker, "Unfiltered Input Data")
  }
  
  # Initialize plot
  p <- ggplot(data, aes(x = .data[[age_var]], 
                        y = .data[[biomarker]], 
                        group = .data[[id_var]], 
                        color = .data[[color_var]])) +
    geom_point(size = point_size) +
    geom_line()
  
  # Add threshold line if provided
  if (!is.null(threshold)) {
    p <- p + 
      geom_hline(yintercept = threshold, color = "dark gray", linetype = "dashed")
    
    if (show_iqr_labels) {
      p <- p +
        annotate("label", x = Inf, y = threshold, label = "Threshold", 
                 hjust = 1, vjust = -0.2, color = "dark gray", size = 3, 
                 fill = "white", alpha = 0.8)
    }
  }
  
  # Add IQR reference lines
  p <- p +
    geom_hline(
      yintercept = c(Q1 - 1.5 * iqr, Q3 + 1.5 * iqr), 
      color = "red", 
      linetype = "dashed",
      alpha = line_alpha
    ) +
    geom_hline(
      yintercept = c(Q1 - 3 * iqr, Q3 + 3 * iqr), 
      color = "red", 
      alpha = max(line_alpha - 0.2, 0.1),  # Ensure alpha doesn't go below 0.1
      linetype = "dashed"
    )
  
  # Add IQR labels if requested
  if (show_iqr_labels) {
    p <- p +
      annotate("label", x = Inf, y = Q3 + 1.5 * iqr, label = "1.5×IQR", 
               hjust = 1, vjust = -0.2, color = "red", size = 3, 
               fill = "white", alpha = line_alpha) +
      annotate("label", x = Inf, y = Q1 - 1.5 * iqr, label = "1.5×IQR", 
               hjust = 1, vjust = 1.2, color = "red", size = 3, 
               fill = "white", alpha = line_alpha) +
      annotate("label", x = Inf, y = Q3 + 3 * iqr, label = "3×IQR", 
               hjust = 1, vjust = -0.2, color = "red", size = 3, 
               fill = "white", alpha = max(line_alpha - 0.2, 0.1)) +
      annotate("label", x = Inf, y = Q1 - 3 * iqr, label = "3×IQR", 
               hjust = 1, vjust = 1.2, color = "red", size = 3, 
               fill = "white", alpha = max(line_alpha - 0.2, 0.1))  # Fixed typo here
  }
  
  # Add labels and theme
  p <- p +
    xlab("Age (years)") +
    ylab(biomarker) +
    ggtitle(plot_title) +
    theme_classic()
  
  return(p)
}



#' Extract Curve Features from Biomarker Trajectories
#'
#' Fits loess smoothed curves to biomarker data and extracts key characteristics
#' including early/late slopes, curvature, and overall trend direction. These
#' features can be used for clustering biomarkers by trajectory shape.
#'
#' @param dt A data.table containing biomarker trajectory data with columns:
#'   \itemize{
#'     \item Biomarker: factor or character identifying the biomarker
#'     \item estdtt0: numeric, estimated time to threshold (x-axis)
#'     \item Score: numeric, biomarker z-score or value (y-axis)
#'   }
#' @param span Numeric, the smoothing span for loess fitting. Larger values 
#'   produce smoother curves. Default is 0.75.
#' @param n_grid Integer, number of points in the time grid for predictions.
#'   More points provide finer resolution. Default is 50.
#'
#' @return A data.table with one row per biomarker containing:
#'   \itemize{
#'     \item Biomarker: biomarker identifier
#'     \item mean_value: average predicted value across trajectory
#'     \item slope_early: linear slope in first third of trajectory
#'     \item slope_late: linear slope in last third of trajectory
#'     \item curvature: mean second derivative (acceleration/deceleration)
#'     \item range: total variation (max - min) in predicted values
#'     \item trend_direction: overall linear trend across full trajectory
#'     \item n_points: number of original data points for this biomarker
#'   }
#'
#' @details
#' The function operates in several steps:
#' \enumerate{
#'   \item For each biomarker, creates a time grid spanning its data range
#'   \item Fits a loess smooth curve to the Score ~ estdtt0 relationship
#'   \item Predicts scores across the time grid
#'   \item Divides the trajectory into thirds to calculate phase-specific slopes
#'   \item Computes curvature as the mean of the second derivative
#'   \item Returns NA for biomarkers with insufficient data (< 10 valid predictions)
#' }
#'
#' @note Biomarkers with fewer than 10 valid predictions or loess fitting errors
#'   will have NA values for all features except n_points.
#'
#' @examples
#' \dontrun{
#' # Assuming a_long is a properly formatted data.table
#' features <- extract_curve_features(a_long, span = 0.75, n_grid = 50)
#' 
#' # Use more smoothing for noisy data
#' features_smooth <- extract_curve_features(a_long, span = 0.9)
#' }
#'
#' @export
extract_curve_features <- function(dt, span = 0.75, n_grid = 50) {
  
  curve_features <- dt[, {
    # Create biomarker-specific time grid
    # This ensures we only predict where we have data
    time_range <- range(estdtt0, na.rm = TRUE)
    time_grid_local <- seq(time_range[1], time_range[2], length.out = n_grid)
    
    tryCatch({
      # Fit loess smooth curve to the data
      fit <- loess(Score ~ estdtt0, data = .SD, span = span)
      
      # Predict scores across the time grid
      pred <- predict(fit, newdata = data.frame(estdtt0 = time_grid_local))
      
      # Remove any NA predictions (rare but possible at edges)
      valid_idx <- !is.na(pred)
      pred_clean <- pred[valid_idx]
      time_clean <- time_grid_local[valid_idx]
      
      # Need at least 10 points for reliable feature extraction
      if(length(pred_clean) < 10) {
        list(
          mean_value = NA, slope_early = NA, slope_late = NA, 
          curvature = NA, range = NA, trend_direction = NA,
          n_points = .N
        )
      } else {
        # Divide curve into thirds for phase-specific slopes
        n <- length(pred_clean)
        third <- floor(n/3)
        
        # EARLY SLOPE: Linear trend in first third of trajectory
        # Captures initial rate of change
        early_idx <- 1:max(2, third)
        slope_early <- if(length(early_idx) > 1) {
          coef(lm(pred_clean[early_idx] ~ time_clean[early_idx]))[2]
        } else NA
        
        # LATE SLOPE: Linear trend in last third of trajectory
        # Captures final rate of change
        late_idx <- max(n-third+1, n-1):n
        slope_late <- if(length(late_idx) > 1) {
          coef(lm(pred_clean[late_idx] ~ time_clean[late_idx]))[2]
        } else NA
        
        # OVERALL TREND: Average rate of change across entire trajectory
        trend_direction <- coef(lm(pred_clean ~ time_clean))[2]
        
        # CURVATURE: Second derivative (measures acceleration/deceleration)
        # Positive = accelerating upward, Negative = decelerating
        if(length(pred_clean) > 2) {
          curvature <- mean(diff(diff(pred_clean)), na.rm = TRUE)
        } else {
          curvature <- NA
        }
        
        # Compile all features
        list(
          mean_value = mean(pred_clean, na.rm = TRUE),  # Average level
          slope_early = slope_early,
          slope_late = slope_late,
          curvature = curvature,
          range = max(pred_clean, na.rm = TRUE) - min(pred_clean, na.rm = TRUE),  # Total variation
          trend_direction = trend_direction,
          n_points = .N  # Number of original data points
        )
      }
    }, error = function(e) {
      # Return NAs if loess fitting fails
      list(
        mean_value = NA, slope_early = NA, slope_late = NA, 
        curvature = NA, range = NA, trend_direction = NA,
        n_points = .N
      )
    })
  }, by = Biomarker]
  
  return(curve_features)
}


#' Test for trajectory differences between hi/lo groups across time bins using LRT
#'
#' @param dat Data.table with the original data
#' @param outcome String specifying the outcome variable (e.g., "AGE")
#' @param exposure String specifying the exposure variable (e.g., "hi.lo")
#' @param time_var String specifying the time variable (e.g., "bin")
#' @param covs Vector of column names containing covariates
#'
#' @returns A named list with LRT results and the alternative model
get_trajectory_lrt <- function(dat, outcome, exposure, time_var, covs = c()) {
  
  # Null model: parallel trajectories (no interaction)
  formula_null <- paste(outcome, "~", paste(c(exposure, time_var, covs), collapse = " + "))
  model_null <- glm(formula_null, data = dat, family = gaussian())
  
  # Alternative model: different trajectories (with interaction)
  formula_alt <- paste(outcome, "~", paste(c(exposure, time_var, 
                                             paste0(exposure, ":", time_var), 
                                             covs), collapse = " + "))
  model_alt <- glm(formula_alt, data = dat, family = gaussian())
  
  # Likelihood ratio test
  lrt <- anova(model_null, model_alt, test = "LRT")
  
  return(list(
    deviance = lrt$Deviance[2],
    df = lrt$Df[2],
    pval = lrt$`Pr(>Chi)`[2],
    nobs = nobs(model_alt),
    model = model_alt  # Return for emmeans if needed
  ))
}


#' Generalized linear model based on a flexible formula
#'
#' @param dat Data.table with the original data
#' @param outcome String specifying the column in dat containing the outcome variable 
#' @param exposure String specifying the column in dat containing the exposure variable
#' @param covs Vector of column names containing covariates to be included in the formula
#' @param return_emmeans Logical, if TRUE also calculates and returns emmeans
#'
#' @returns A named list: 'beta', 'se', 't stat', 'pval', 'nobs', and optionally 'emmeans' (as a list column)
get_linreg = function(dat, outcome, exposure, covs=c(), return_emmeans=FALSE){
  model = glm(paste(outcome, "~", paste(c(exposure, covs), collapse=" + ")),
              na.action=na.exclude, data=dat)
  a = summary(model)
  
  result <- as.list(c(a$coefficients[2, ], nobs=nobs(model)))
  
  if(return_emmeans) {
    # Calculate emmeans and store as a list
    emm <- emmeans::emmeans(model, specs = exposure)
    result$emmeans <- list(as.data.frame(emm))
  }
  
  return(result)
}

#' function to make a barplot with swarm points given a fitted linear model
#' 
#' @params dt data.table with data
#' @params emm data.table with estimated marginal means from emmeans 
#' @params x string, x column data tables specifying group
#' @params y string, y column for dt only, emm depends on fixed emmeans values  
#' @params col string, grouping column
#' @params titl string, plot title
#' @params pal vector, list of color palette values for scale_color_manual
#' 
#' @return ggplot object
bar_swarm <- function(dt, emm, x, y, col, titl, 
                      pal = RColorBrewer::brewer.pal("Dark2", n = 7)) {
  # Convert emmeans to data.frame if needed
  if("emmGrid" %in% class(emm)) {
    emm_df <- as.data.frame(emm)
  } else {
    emm_df <- emm
  }

  p <- ggplot(dt, aes(x = .data[[x]], y = .data[[y]], color = .data[[col]], fill = .data[[col]])) +
    # Horizontal line at emmean using geom_errorbar
    geom_jitter(
      data = dt,
      aes(x = .data[[x]], y = .data[[y]]),
      position = position_jitter(width = 0.2, seed = 123), 
      size = 2,
      color = "gray"
    ) +
    geom_errorbar(
      data = emm_df,
      aes(x = .data[[x]], ymin = emmean, ymax = emmean, color = .data[[x]]),
      width = 0.6,
      linewidth = 2,
      inherit.aes = FALSE
    ) +
    # Error bars using emmean ± SE
    geom_errorbar(
      data = emm_df,
      aes(x = .data[[x]], y = emmean, ymin = emmean - SE, ymax = emmean + SE, color = .data[[x]]),
      width = 0.3,
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    theme_classic() +
    ggtitle(titl) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    theme(
      axis.text.x = element_text(
        angle = 45, 
        hjust = 1, 
        vjust = 1, 
        family = 'Arial',
        face = 'bold',
        size = 12
      ),
      axis.text.y = element_text(
        family = 'Arial',
        face = 'bold',
        size = 8
      ),
      axis.title.y = element_text(
        face = 'bold', 
        family = 'Arial', 
        size = 12
      ), 
      plot.title = element_text(
        hjust = 0.5, 
        face = 'bold', 
        family = 'Arial', 
        size = 12
      ), 
      plot.subtitle = element_text(
        hjust = 0.5, 
        face = 'bold', 
        family = 'Arial', 
        size = 8
      ),
      legend.position = 'none',
      axis.title.x = element_blank()
    )
  return(p)
}


# Main -------------------------------------------------------------------------
a <- fread(input, colClasses = c(PIDN = "character"))


# Shape data for analysis ---------
a <- a[, .(PIDN, SampleID, aa, Age, Gender, DS1V), env = envs
       ][order(PIDN, DS1V)]

# re-age the groups to make sure the decimal years (since we have days since first visit)
a[, AGE := Age[which.min(Age)] + DS1V / 365.2425, by = PIDN]
a[, ever_over := aa[which.max(aa)] > pTau_thresh, by = PIDN, env = envs]

# drop missing plasma values
a <- a[!is.na(aa), env = envs]

# spaghetti plot of value vs. age for simulated data
Q1 <- quantile(a[[biomarker]], 0.25)
Q3 <- quantile(a[[biomarker]], 0.75)
iqr <- IQR(a[[biomarker]])
p1 <- plot_biomarker_outliers(
  data = a,
  biomarker = biomarker,
  threshold = pTau_thresh,
  plot_title = paste(biomarker, "Unfiltered Input Data")
)

# plot z zcores
a[, zscore := as.vector(scale(aa)), env = envs]
p2 <- ggplot(a, aes(x = AGE, y = zscore, group = PIDN, color = ever_over)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = c(-3, 3), color = "dark gray", linetype = "dashed") +
  xlab("Age (years)") +
  ggtitle(paste(biomarker, "Z score deviation")) +
  theme_classic() 

if (outliers == "cooksd") {
  # looking at cooks distance to examine outliers
  mod <- lm(as.formula(paste(biomarker, "AGE", sep = " ~ ")), data = a)
  a[, cooksd := cooks.distance(mod)]

  # Plot the Cook's Distance using the traditional 4/n criterion
  p2 <- ggplot(a, aes(x = AGE, y = cooksd, group = PIDN, color = ever_over)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 4/nrow(a), color = "dark gray", linetype = "dashed") +
    xlab("Age (years)") +
    ggtitle(paste(biomarker, "Cooks distance from model"),
      sub = deparse(formula(mod))
    ) +
    theme_classic() 
  a <- a[cooksd < 4/nrow(a)]
} else if (outliers == "zscore"){
  # drop outlier plasma values
  a <- a[abs(zscore) <= 3]
} else {
  a <- a[aa > Q1 - 1.5 * iqr & aa < Q3 + 1.5 * iqr, env = envs]
}

ggsave(
  filename = file.path(output, 
    paste(nameset, biomarker, "unfiltered_input.png", sep = ".")
  ),
  plot = p1 + p2 + plot_layout(guides = 'collect'),
  width = 12,
  height = 6
)


# Modeling the data ---------------
res <- sila(
  df = tibble(subid = a[['PIDN']], age = a[['AGE']], val = a[[biomarker]]),
  dt = 0.25, 
  val0 = pTau_thresh, 
  maxi = 200
)

resfit <- sila_estimate(
  tsila = res$tsila, 
  df = tibble(subid = a[['PIDN']], age = a[['AGE']], val = a[[biomarker]]),
  align_event = "last"
)

a <- merge(
  x = a, 
  y = as.data.table(resfit), 
  by.x = c("PIDN", "AGE", biomarker),
  by.y = c("subid", "age", "val")
)

fwrite(a, 
  file = file.path(output, 
    paste(nameset, biomarker, "fitted_ages.csv", sep = ".")
  )
)


# Plotting  ----------------------
# plots showing the output from discrete rate sampling (i.e., rate vs. value)
# and modeled value vs. time data.
p1 <- ggplot(res$tdrs, aes(x = val, y = rate)) +
  geom_ribbon(aes(ymin = rate - ci, ymax = rate + ci), alpha = .3) +
  geom_line() +
  xlab("Value") +
  ylab(expression(Delta * "Value per Year")) +
  ggtitle("Discrete Rate Sampling Curve") +
  theme_bw()

p2 <- ggplot(res$tsila, aes(x = adtime, y = val)) +
  geom_line(aes(color = "Modeled curve")) +
  geom_hline(aes(yintercept = pTau_thresh, color = "threshold"), linetype = "dashed") +
  xlab("Time from Threshold") +
  ylab("Value") +
  ggtitle(
    expression(paste("SILA Modeled ", italic("Value vs. Time"), " Curve"))
  ) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(.15, .85)
  )

ggsave(
  filename = file.path(output, 
    paste(nameset, biomarker, "model_output.png", sep = ".")
  ),
  plot = p1 + p2,
  width = 12,
  height = 6
)

# spaghetti plot of value vs. age for simulated data
p1 <- plot_biomarker_outliers(
  data = a,
  biomarker = biomarker,
  threshold = pTau_thresh,
  plot_title = paste(biomarker, "Filtered Input Data"),
  iqr_values = c(Q1, Q3, iqr)
)

# value vs. time for all subjects
p2 <- ggplot(a, 
    aes(x = estdtt0, y = .data[[biomarker]], group = PIDN, color = ever_over)
  ) +
  geom_point() +
  geom_line() +
  geom_line(
    data = res$tsila,
    aes(x = adtime, y = val, group = NULL, color = NULL)
  ) +
  geom_hline(yintercept = pTau_thresh, color = "dark gray", linetype = "dashed") +
  xlab("Estimated time to threshold (years)") +
  ylab(biomarker) +
  ggtitle(paste(biomarker, "Aligned by Estimated Time to Threshold")) +
  theme_classic()

ggsave(
  filename = file.path(output, 
    paste(nameset, biomarker, "fitted_ages.png", sep = ".")
  ),
  plot = p1 + p2 + plot_layout(guides = 'collect'),
  width = 12,
  height = 6
)

ggsave(
  filename = file.path(output, 
    paste(nameset, biomarker, "fitted_plot.pdf", sep = ".")
  ),
  plot = p2,
  width = 12,
  height = 6,
  device = cairo_pdf
)

# relationship to eigengenes -----------------------------------------------
# Load eigengenes
x <- fread(eigengenes)
a <- merge(a, x, by = "SampleID")

# load term names 
names <- fread(ref)
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

# Define biomarker columns (all the module columns)
biomarker_cols <- names[mean.Zsum > 2.05][order(group, -mean.Zsum), name]

# Melt the data from wide to long
a_long <- melt(a[ever_over == TRUE], 
  id.vars = c("SampleID", "PIDN", "estdtt0", "estaget0", "AGE", "Gender"),
  measure.vars = biomarker_cols,
  variable.name = "Biomarker",
  value.name = "Score"
)
# # filter to outer quartiles
# a_long[, quartile := dplyr::ntile(Score, 4), by = Biomarker]
# a_long <- a_long[quartile %in% c(1, 4)]



# Extract curve features using the function
curve_features <- extract_curve_features(a_long, span = 0.75, n_grid = 50)

# Check for NAs
print("NA counts:")
print(colSums(is.na(curve_features)))

# Keep only biomarkers with complete features
curve_features_complete <- curve_features[complete.cases(
  curve_features[, .(slope_early, slope_late, curvature, trend_direction)]
)]

print(paste("Biomarkers with complete features:", 
            nrow(curve_features_complete), "out of", nrow(curve_features)))

# Prepare feature matrix for clustering
feature_matrix <- as.matrix(curve_features_complete[, 
  .(slope_early, slope_late, curvature, trend_direction)
])
rownames(feature_matrix) <- curve_features_complete$Biomarker

# Standardize features (mean=0, sd=1)
feature_matrix <- scale(feature_matrix)

# HIERARCHICAL CLUSTERING
# Use Ward's method which minimizes within-cluster variance
hc <- hclust(dist(feature_matrix), method = "ward.D2")

# Visualize dendrogram
CairoPDF(file = file.path(output, 
    paste(nameset, biomarker, "hclust.pdf", sep = ".")
))
plot(hc, main = "Hierarchical Clustering of Biomarker Trajectories",
     xlab = "Biomarker", ylab = "Height", cex = 0.7)

# Cut tree into k clusters (try 3-5)
k <- 5
rect.hclust(hc, k = k, border = 2:5)
dev.off()

# Assign clusters
curve_features_complete[, cluster := cutree(hc, k = k)]

# Characterize each cluster
cluster_summary <- curve_features_complete[, .(
  mean_early_slope = mean(slope_early),
  mean_late_slope = mean(slope_late),
  mean_curvature = mean(curvature),
  mean_trend = mean(trend_direction),
  slope_change = mean(slope_late - slope_early),  # How much slope changes
  n = .N
), by = cluster][order(mean_trend)]

print("Cluster characteristics:")
print(cluster_summary)

# Assign descriptive names based on cluster characteristics
cluster_names <- data.table(
  cluster = cluster_summary$cluster,
  trajectory_name = c()  # Will be filled based on patterns
)

# Name clusters based on their characteristics
for(i in 1:nrow(cluster_summary)) {
  clust <- cluster_summary[i]
  
  if(clust$mean_trend > 0.015) {
    name <- "Strongly Increasing"
  } else if(clust$mean_trend < -0.015) {
    name <- "Strongly Decreasing"
  } else if(abs(clust$slope_change) > 0.02) {
    name <- ifelse(clust$slope_change > 0, "Accelerating", "Decelerating")
  } else {
    name <- "Stable/Flat"
  }
  
  cluster_names[i, trajectory_name := name]
}

# Merge names back
curve_features_complete <- merge(
  curve_features_complete, 
  cluster_names, 
  by = "cluster"
)

# Merge trajectory assignment back to long data
a_long <- merge(
  a_long, 
  curve_features_complete[, .(Biomarker, trajectory_name, cluster)], 
  by = "Biomarker", 
  all.x = TRUE
)

# Remove biomarkers without trajectory assignment
a_long <- a_long[!is.na(trajectory_name)]

# add names
a_long <- merge(a_long, names, by.x = 'Biomarker', by.y = 'name')
a_long[, bin := cut(estdtt0, breaks = c(-15, -10, -5, 0, 5, 10, 15))]

# Create faceted plot by trajectory type
p_all <- ggplot(a_long, aes(x = estdtt0, y = Score, fill = Biomarker, color = Biomarker)) +
  geom_smooth(method = "loess", se = TRUE, formula = y ~ x, span = 0.75, alpha = 0.1) +
  facet_wrap(~ trajectory_name, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(x = "Estimated Time to Symptom Onset (years)",
       y = "Module Score",
       title = "Biomarker Trajectories Grouped by Hierarchical Clustering") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7),
        strip.text = element_text(size = 10, face = "bold"))

print(p_all)

# Print detailed cluster membership
print("\nBiomarkers by trajectory type:")
print(curve_features_complete[order(cluster, Biomarker), 
  .(Biomarker, trajectory_name, slope_early, slope_late, curvature, trend_direction)
])

# Create individual plots for each trajectory type
trajectory_types <- unique(a_long$trajectory_name)
plist <- lapply(trajectory_types, function(traj) {
  ggplot(a_long[trajectory_name == traj], 
              aes(x = estdtt0, y = Score, color = term, fill = term)) +
    # geom_point()+
    geom_smooth(method = "loess", se = FALSE, formula = y ~ x, span = 0.75) +
    theme_classic() +
    labs(x = "Estimated Time to Symptom Onset (years)",
         y = "Module Score",
         title = paste0(traj, " Trajectories")) +
    theme(legend.position = "right")
})

wrap_plots(plist)
ggsave(
  plot = wrap_plots(plist),
  filename = file.path(output, 
    paste(nameset, biomarker, "trajectory_plots.pdf", sep = ".")
  ),
  width = 16,
  height = 8,
  device = cairo_pdf
)

# Example plots
ex_terms <- c(
  "Autophagy", "Protein Ubiquitination", "BMP signaling", 
  # "Glutathione metabolism|Oxidative stress", 
  "Synaptic vesicle|SNARE complex",
  "tRNA aminoacylation"
)

# Get the endpoint of each curve for label placement
label_data <- a_long[term %in% ex_terms, 
  .(estdtt0 = max(estdtt0), 
    Score = Score[which.max(estdtt0)]), 
  by = term]

p1 <- ggplot(a_long[term %in% ex_terms & estdtt0 > -10], 
        aes(x = estdtt0, y = Score, color = term, fill = term)) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.15, formula = y ~ x, span = 0.75) +
    geom_text(data = label_data, 
              aes(x = estdtt0, y = Score, label = term),
              hjust = 0, nudge_x = 0.5, size = 3.5) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    theme_classic() +
    labs(x = "Estimated Time to Symptom Onset (years)",
         y = "Module Score",
         title = "Example Trajectories of WSS Module Scores") +
    coord_cartesian(clip = "off") +  # Allow labels to extend beyond plot
    theme(plot.margin = margin(5.5, 80, 5.5, 5.5),
      legend.position = "none",
      text = element_text(family = "Arial", face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold")
    )  # Add right margin for labels

ggsave(
  plot = p1,
  filename = file.path(output, 
    paste(nameset, biomarker, "example_modules.pdf", sep = ".")
  ),
  width = 10,
  height = 6,
  device = cairo_pdf
)


# Testing for bin specific signals of earlier onset ---------------------
# need to adjust for Gender, PCs, weak N though
covs <- c("Gender")

# Map in PCs assumes your primary data table exists 
if (pcs != "") {
  print(sprintf("Reading PCs for inclusion from %s", pcs))
  x <- fread(pcs, colClasses = c(PIDN = "character"))
  a_long <- merge(a_long, x, by = "PIDN")
  covs <- c(covs, names(a_long)[names(a_long) %like% "^PC\\d{1,2}$"]) # flexible for PCs
}

# filter for lowest bin (n = 1)
a_filt <- a_long[bin != "(-15,-10]"][, bin := droplevels(bin)]

# binarize module per estdtt0 bin. Is being in the upper half significantly lower estaget0.
a_filt[, hi.lo := ifelse(Score > median(Score), "hi", "lo"), by = .(Biomarker, bin)]

# trajectory tests for significant biomarkers ---
dt_lrt <- CJ(Biomarker = unique(a_filt$Biomarker))
dt_lrt[, c("deviance", "df", "pval", "nobs") := {
  get_trajectory_lrt(
    a_filt[Biomarker == dt_lrt$Biomarker[.I]],
    outcome = "AGE",
    exposure = "hi.lo",
    time_var = "bin",
    covs = covs[1:3]
  )
}, by = Biomarker]

# Adjust for multiple testing across biomarkers
dt_lrt[, padj := p.adjust(pval, method = "fdr")]

# individual bin tests for sig diff ---
dt1 <- CJ(Biomarker = a_filt$Biomarker, bin = a_filt$bin, unique = TRUE)
dt1[, Sig_Biomarker := Biomarker %in% dt_lrt[padj < 0.1, Biomarker]]

# fast run function
dt1[Sig_Biomarker == TRUE, c("beta", "se", "t stat", "pval", "nobs", "emmeans") := {
  get_linreg(a_filt[Biomarker == dt1$Biomarker[.I] & bin == dt1$bin[.I]], 
    outcome = "AGE", 
    exposure = "hi.lo",
    covs = covs[1:3],
    return_emmeans = TRUE
  )
}, by=.I]
dt1[, padj := p.adjust(pval, method = "bonferroni"), by = Biomarker]

# Create plots of emmeans
sig <- dt1[padj < 0.05]
y_range <- range(a_filt[Biomarker %in% sig$Biomarker, AGE], na.rm = TRUE)

plot_list <- lapply(1:nrow(sig), function(i) {
  bio <- sig$Biomarker[i]
  bin_val <- sig$bin[i]
  term <- unique(a_filt[Biomarker == bio, term])
  
  # Get the data subset
  plot_data <- a_filt[Biomarker == bio & bin == bin_val]
  
  # Get the stored emmeans
  emm_df <- sig$emmeans[[i]]
  
  # Create plot
  bar_swarm(
    dt = plot_data,
    emm = emm_df,
    x = "hi.lo",
    y = "AGE",
    col = "hi.lo",
    titl = paste0(term, "\n", bio, " | ", bin_val, "\npadj = ", 
                  format(sig$padj[i], digits = 3)),
    pal = RColorBrewer::brewer.pal(n = 3, name = "Dark2")[-2]
  ) + 
    coord_cartesian(ylim = y_range)
})

# View all plots
ggsave(
  plot = wrap_plots(plot_list, ncol = nrow(sig), axes = 'collect_y'),
  filename = file.path(output, 
    paste(nameset, biomarker, "hi-lo_bins.pdf", sep = ".")
  ),
  device = cairo_pdf, 
  width = 16
)


# Extract emmeans data for each row
dt1[, `:=`(
  emmean.hi = sapply(emmeans, function(x) if(is.data.frame(x) && nrow(x) >= 2) x$emmean[x$hi.lo == "hi"] else NA_real_),
  emmean.lo = sapply(emmeans, function(x) if(is.data.frame(x) && nrow(x) >= 2) x$emmean[x$hi.lo == "lo"] else NA_real_),
  SE.hi = sapply(emmeans, function(x) if(is.data.frame(x) && nrow(x) >= 2) x$SE[x$hi.lo == "hi"] else NA_real_),
  SE.lo = sapply(emmeans, function(x) if(is.data.frame(x) && nrow(x) >= 2) x$SE[x$hi.lo == "lo"] else NA_real_)
)]

# Remove the emmeans list column
dt1[, emmeans := NULL]

# Write to CSV
fwrite(dt1, file = file.path(output, 
  paste(nameset, biomarker, "hi-lo_pairwise_regression.csv", sep = ".")
))
fwrite(dt_lrt, file = file.path(output, 
  paste(nameset, biomarker, "hi-lo_lrt_trajectory.csv", sep = ".")
))


# Heatmap of progression ------------------------------------------------
# Calculate mean score per biomarker per bin
heatmap_data <- a_long[, .(mean_score = mean(Score, na.rm = TRUE)), 
                       by = .(Biomarker, bin, term, trajectory_name)]

# Remove any NAs and exclude the (-15,-10] bin
heatmap_data <- heatmap_data[!is.na(bin) & bin != "(-15,-10]"]

# Order biomarkers according to biomarker_cols
heatmap_data[, Biomarker := factor(Biomarker, levels = biomarker_cols)]
heatmap_data <- heatmap_data[order(Biomarker)]

# Order bins chronologically (without the (-15,-10] bin)
heatmap_data[, bin := factor(bin, 
  levels = rev(c("(-10,-5]", "(-5,0]", "(0,5]", "(5,10]", "(10,15]"))
)]

# Calculate number of unique individuals per bin (excluding (-15,-10])
bin_counts <- a_long[!is.na(bin) & bin != "(-15,-10]", 
                     .(n_individuals = uniqueN(SampleID)), by = bin]
bin_counts[, bin := factor(bin, 
  levels = rev(c("(-10,-5]", "(-5,0]", "(0,5]", "(5,10]", "(10,15]"))
)]

# Create trajectory annotation data
trajectory_anno <- unique(heatmap_data[, .(Biomarker, trajectory_name)])
trajectory_anno[, Biomarker := factor(Biomarker, levels = biomarker_cols)]

# Get PuOr colors
puor_colors <- brewer.pal(9, "PuOr")

# Main heatmap
p_heatmap_anno <- ggplot(heatmap_data, aes(x = Biomarker, y = bin, fill = mean_score)) +
  geom_tile(color = "white", linewidth = NA) +
  # Add borders for pval significance (gray)
  geom_tile(
      data = dt1[pval < 0.05],
      aes(x = .data[["Biomarker"]], y = .data[["bin"]]),
      fill = NA,
      color = "gray50",
      linewidth = 0.8
  ) +
  # Add borders for padj significance (black)
  geom_tile(
      data = dt1[padj < 0.05],
      aes(x = .data[["Biomarker"]], y = .data[["bin"]]),
      fill = NA,
      color = "black",
      linewidth = 0.8
  ) +
  geom_hline(yintercept = 3.5, linetype = "dashed") +
  scale_fill_gradient2(
    low = puor_colors[9],      # Purple for most negative
    mid = puor_colors[5],      # White at 0
    high = puor_colors[1],     # Orange for most positive
    midpoint = 0,
    name = "Mean\nScore"
  ) +
  scale_y_discrete(
    labels = function(x) {
      counts <- bin_counts$n_individuals[match(x, bin_counts$bin)]
      paste0(x, " (n=", counts, ")")
    }
  ) +
  theme_classic(base_family = "Arial") +
  labs(x = "Module", 
       y = "Time to Threshold (years)",
       title = "Module Score Trajectories Across Disease Progression") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                               size = 8, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  ) +
  scale_x_discrete(labels = setNames(heatmap_data$term, heatmap_data$Biomarker))

# Trajectory annotation bar
p_anno_bar <- ggplot(trajectory_anno, aes(x = Biomarker, y = 1, fill = trajectory_name)) +
  geom_tile(color = "white", linewidth = 0) +
  scale_fill_brewer(palette = "Dark2", name = "Trajectory") +
  scale_x_discrete(labels = setNames(heatmap_data$term, heatmap_data$Biomarker)) +
  theme_void(base_family = "Arial") +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 9, face = "bold", family = "Arial"),
    legend.title = element_text(size = 10, face = "bold", family = "Arial"),
    axis.text.x = element_blank(),
    plot.margin = margin(5.5, 5.5, 0, 5.5)
  )

# Combine plots
p_combined <- p_anno_bar / p_heatmap_anno + plot_layout(heights = c(1, 10))

print(p_combined)

ggsave(
  plot = p_combined,
  filename = file.path(output, 
    paste(nameset, biomarker, "trajectory_heatmap.pdf", sep = ".")
  ),
  width = 14,
  height = 4.5,
  device = cairo_pdf
)

# Follow-up M8 ---------------------------------------------------------------
# Filter to Merged.M8 and remove NA emmeans
m8_data <- dt1[Biomarker == "Merged.M8" & !is.na(emmean.hi)]

# Reshape to long format for plotting
m8_long <- rbind(
  m8_data[, .(bin, group = "hi", emmean = emmean.hi, SE = SE.hi, padj)],
  m8_data[, .(bin, group = "lo", emmean = emmean.lo, SE = SE.lo, padj)]
)

# Create numeric x-axis position for bins
bin_order <- c("(-10,-5]", "(-5,0]", "(0,5]", "(5,10]", "(10,15]")
m8_long[, bin_num := match(bin, bin_order)]

# Create p-value labels
m8_long[, padj_label := ifelse(!is.na(padj), 
                                sprintf("padj = %.3f", padj), 
                                "")]

# Plot
p_m8 <- ggplot(m8_long, aes(x = bin_num, y = emmean, color = group, group = group)) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, linewidth = 0.8) +
  geom_text(aes(label = padj_label), 
            vjust = -1.5, size = 3.5, 
            color = "black", fontface = "bold",
            family = "Arial") +
  scale_x_continuous(breaks = 1:5, labels = bin_order) +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 3, name = "Dark2")[-1],
                     labels = c("High", "Low"),
                     name = "Module Score") +
  theme_classic(base_family = "Arial") +
  labs(x = "Time to Threshold (years)",
       y = "Age at Visit",
       title = "Merged.M8 Module Trajectory by Score Group") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                               size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = c(0.95, 0.05),  # Position in lower right (x=0.85, y=0.15)
    legend.justification = c(1, 0),   # Anchor at bottom-right of legend box
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 11, face = "bold")
  )

print(p_m8)

# Save
ggsave(
  plot = p_m8,
  filename = file.path(output, 
    paste(nameset, biomarker, "M8_trajectory.pdf", sep = ".")
  ),
  width = 4,
  height = 6,
  device = cairo_pdf
)
