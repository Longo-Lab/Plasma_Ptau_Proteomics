#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 02-linreg_glm.R
#'
#' Version: 0.0.8
#'
#' Purpose of script: calculate correlation of plasma ptau with SomaLogic plasma
#' proteomic levels
#'
#' Author: Robert B Butler III, Michael P Brown, Edward N Wilson III
#'
#' Date Created: 2023-08-25
#'
#' Copyright (c) 2023
#' Email: rrbutler@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'   v0.0.8 - fixing logging bugs and adding quasi-complete separation reporting table
#'   v0.0.7 - adding cohen's d and fixing logging bugs
#'   v0.0.6 - adding fallback reporting for quasi-complete separation
#'   v0.0.5 - adding logistic regression option with lme4qtl, relmatGlmer
#'   v0.0.3 - incorporated the lme4qtl model with p-vals from LRT
#'   Linear regression glm for baseline visit only with corrections
#'   for fixed effects: filter covariates based on importance
#'
#'   PlateRunDate & ScannerID couldn't be fit to the model, suggesting some
#'   other issue, so dropping those
#'
#'   This takes a while, so batch unless you like waiting
#'   Usage:
#'     sbatch -J ptau --mem=12G -c 16 -t 06:00:00 -p interactive \
#'       -o %A_02-linreg_lme4.log \
#'       --wrap "ml R/4.2.2; Rscript 02-linreg_glm.R"
#'
#'   REMEMBER, WHEN COMMENTING IN CODE, DO NOT REFERENCE INDIVIDUAL INFORMATION:
#'   SAMPLES/PIDNS, DATES OR ANY INDIVIDUAL LEVEL INFO
#'
#' ---------------------------

#' load up the packages we will need:  (uncomment as required)
suppressPackageStartupMessages({
  library(tictoc)
  tic()
  tic()
  library(data.table)
  library(tools)
  library(stringr)
  library(lme4)
  # library(lmerTest) # loads lme4
  library(lme4qtl)
  library(GENESIS)
  library(pROC)
  library(caret)
  library(future.apply)
  options(future.cmdargs = character(0))
  library(futile.logger)
  library(rstatix)
  library(optparse)
})

# Nameset helper function
nameset_default <- function(input) {
  file_path_sans_ext(basename(input))
}

# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = ".",
    help = "File path to DE input file (*.csv.gz)[default: %default]"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = ".",
    help = "File path to output directory [default: %default]"
  ),
  make_option(
    c("-p", "--prots"),
    type = "character",
    default = "prot_names.txt",
    help = "File path to prots file [default: %default]"
  ),
  make_option(
    c("-r", "--ref"),
    type = "character",
    default = "genomic/mypcrelate.rds",
    help = "File path to directory of reference data [default: %default]"
  ),
  make_option(
    c("-v", "--covs"),
    type = "character",
    default = "",
    help = "Terminal input of covariates in comma-separated format [default: %default]"
  ),
  make_option(
    c("-s", "--scale-cov"),
    type = "character",
    dest = "scale_cov",
    default = "",
    help = "Which covariates to scale for relmatGlmer in comma-separated format [default: %default]"
  ),
  make_option(
    c("-d", "--recodes"),
    type = "character",
    default = "",
    help = "Recode covariates in comma-separated format and pipe between covs [default: %default]"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project [default: Nameset of input]"
  ),
  make_option(
    c("-t", "--threads"),
    default = 0,
    help = "number of threads [default: setDTthreads]"
  ),
  make_option(
    c("-b", "--binomial"),
    action = "store_true",
    default = FALSE,
    help = "Use binomial family for logistic regression [default: %default]"
  ),
  make_option(
    c("-f", "--factors"),
    default = "Q1,Q4",
    help = "when binomial == TRUE, what are the levels (ordered) for comparison [default: %default]"
  ),
  make_option(
    c("-c", "--outcome"),
    type = "character",
    help = "String specifying the column in data containing the outcome variable [required]"
  )
)

args_list <- parse_args(OptionParser(option_list = option_list))

# Nameset default case
if (is.null(args_list$nameset)) {
  args_list$nameset <- nameset_default(args_list$input)
}

# check outpath status
if (!dir.exists(args_list$output)) {
  dir.create(args_list$output, recursive = TRUE)
}

# Set up logger
log_fh <- file.path(
  args_list$output, paste(args_list$nameset, "linreg.log", sep = ".")
)
invisible({
  file.create(log_fh)
  flog.layout(layout.format('[~l] ~m'))
  flog.appender(appender.file(log_fh))
  flog.threshold(INFO)
})

# read in fields
fields <- c(
  "input",
  "output",
  "prots",
  "ref",
  "covs",
  "scale_cov",
  "recodes",
  "nameset",
  "threads",
  "binomial",
  "factors",
  "outcome"
)

# check for optparse values
missing_fields <- setdiff(fields, names(args_list))

if (length(missing_fields) > 0) {
  flog.fatal("Missing required arguments: %s", paste(missing_fields, collapse = ", "))
  stop(parse_args(OptionParser(option_list = option_list), args = "--help"),
    call. = FALSE
  )
}

# Log all provided arguments
flog.info("=== Command-line arguments ===")
for (field in fields) {
  assign(field, args_list[[field]])
  flog.info("  %s: %s", field, get(field))
}
flog.info("==============================")

# Validate and set threads
if (threads == 0) {
  threads <- as.integer(parallelly::availableCores() - 1)
} else if (!is.numeric(threads)) {
  flog.fatal("threads must be an integer value, got: %s", class(threads))
  stop(parse_args(OptionParser(option_list = option_list), args = "--help"),
    call. = FALSE
  )
} else {
  threads <- as.integer(threads)
}

flog.info("Setting up parallel processing with %d workers", threads)
plan(multicore, workers = threads)

# Split list variables
covs <- strsplit(covs, ",")[[1]]
scale_cov <- strsplit(scale_cov, ",")[[1]]
prots <- readLines(prots)
recodes <- strsplit(unlist(strsplit(recodes, "\\|")), ",")
factors <- strsplit(factors, ",")[[1]]


# Functions ----------------------------------------------------------------
#' Recodes dummy categorical variables and can subset by vector of categories
#'
#' @param dat Data.table with the original data
#' @param dcol String representing the column to be split
#' @param keep Vector of strings that are the categories within dcol to extract. default is all c()
#'
#' @returns modified data.table with new dummy variable columns named "dcol.category" after original column
recode_dummies <- function(dat, dcol, keep = c()) {
  dt <- dat

  # check input
  if (length(keep) == 0) {
    inds <- unique(dt[[dcol]])
  } else if (length(keep) > 0 & all(keep %in% unique(dt[[dcol]]))) {
    inds <- unique(dt[[dcol]])[unique(dt[[dcol]]) %in% keep]
  } else {
    flog.error("Invalid 'keep' argument for column '%s'. Provided: %s, Available: %s", 
               dcol, paste(keep, collapse = ", "), paste(unique(dt[[dcol]]), collapse = ", "))
    stop("Keep must be a vector of names in the original column")
  }

  # split out dummies
  dt[, c(paste(dcol, inds, sep = ".")) := lapply(inds, function(x) ifelse(get(dcol) == x, 1, 0))]
  # reorder behind original variable
  v <- grep(paste0("^", dcol, "$"), names(dt))[1]
  inds <- paste(dcol, inds, sep = ".")
  setcolorder(dt, c(names(dt)[1:v], inds, names(dt)[(v + 1):(ncol(dt) - length(inds))]))
  return(dt)
}


#' Check for quasi-complete separation in logistic regression
#'
#' Evaluates whether an exposure variable exhibits quasi-complete or complete
#' separation with respect to a binary outcome. This is done by creating a
#' binary split of the exposure at its median and checking if any cell in the
#' resulting 2x2 contingency table has 2 or fewer observations.
#'
#' @param dat Data.table or data.frame containing the data
#' @param outcome String specifying the column name of the binary outcome variable
#' @param exposure String specifying the column name of the exposure variable
#'
#' @return Logical. TRUE if quasi-complete separation is detected (min cell count <= 2),
#'   FALSE otherwise
#'
#' @details
#' Quasi-complete separation occurs when the outcome categories can be nearly
#' perfectly predicted by the exposure variable. This causes logistic regression
#' models to fail to converge as coefficient estimates approach infinity.
#' 
#' The function uses a median split of the exposure variable as a simple heuristic.
#' A minimum cell count of 2 or fewer in the 2x2 table is used as the threshold
#' for flagging separation issues.
#'
#' @examples
#' \dontrun{
#' # Check if protein shows separation with outcome
#' is_separated <- check_separation(dat, outcome = "disease_status", 
#'                                   exposure = "protein_ABC")
#' }
#'
#' @export
check_separation <- function(dat, outcome, exposure) {
  exposure_binary <- dat[[exposure]] > median(dat[[exposure]], na.rm = TRUE)
  tab <- table(dat[[outcome]], exposure_binary)
  min_cell <- min(tab)
  
  return(min_cell <= 2)
}


#' Generalized linear model based on a flexible formula
#'
#' @param dat Data.table with the original data
#' @param m0 lmerMod baseline model with without the exposure element (linear only)
#' @param outcome String specifying the column in dat containing the outcome variable
#' @param exposure String specifying the column in dat containing the exposure variable
#' @param covs Vector of column names containing covariates to be included in the formula
#' @param mat dsCMatrix a Generalized relatedness matrix
#' @param logistic booliean indicating logistic regression instead of gaussian (default)
#'
#' @returns A list: 'Estimate' (beta), 'Std. Error' (se), 't value' (test statistic), 'Pr(>Chisq)' (p-value), 'nobs' (number of observations)
get_linreg <- function(dat, m0 = NULL, outcome, exposure, covs = c(), mat, logistic = FALSE) {
  
  log_msgs <- character(0)  # Collect log messages
  
  if (logistic) {
    # Check for quasi-complete separation
    if (check_separation(dat, outcome, exposure)) {
      log_msgs <- c(log_msgs, sprintf("Quasi-complete separation detected for %s.", exposure))
      return(list(
        Estimate = NA,
        `Std. Error` = NA,
        `z value` = NA,
        `Pr(>|z|)` = 0,
        nobs = nrow(dat),
        flag = "separated",
        logs = log_msgs
      ))
    }
    
    # Try with random effect
    m <- tryCatch({
      suppressWarnings({
        relmatGlmer(
          as.formula(paste(outcome, "~", paste(c(exposure, covs, "(1|PIDN)"), collapse = " + "))),
          data = dat,
          relmat = list(PIDN = mat),
          family = binomial()
        )
      })
    }, error = function(e) {
      log_msgs <<- c(log_msgs, sprintf("relmatGlmer failed for %s: %s", exposure, e$message))
      NULL
    })
    
    # Check for convergence issues
    use_glm <- FALSE
    if (is.null(m)) {
      use_glm <- TRUE
    } else if (isSingular(m)) {
      log_msgs <- c(log_msgs, sprintf("Singular fit detected for %s, using GLM", exposure))
      use_glm <- TRUE
    } else {
      # Check for convergence warnings
      if (!is.null(m@optinfo$conv$lme4$messages)) {
        log_msgs <- c(log_msgs, sprintf("Convergence issues for %s: %s. Using GLM", 
                       exposure, 
                       paste(m@optinfo$conv$lme4$messages, collapse = "; ")))
        use_glm <- TRUE
      }
    }
    
    # Fall back to GLM if needed
    if (use_glm) {
      m <- tryCatch({
        glm(
          as.formula(paste(outcome, "~", paste(c(exposure, covs), collapse = " + "))),
          data = dat,
          family = binomial()
        )
      }, error = function(e) {
        log_msgs <<- c(log_msgs, sprintf("GLM also failed for %s: %s", exposure, e$message))
        NULL
      })
      
      fl <- "glm"
    } else {
      fl <- "relmatGlmer"
    }
    
    # If both models failed
    if (is.null(m)) {
      return(list(
        Estimate = NA,
        `Std. Error` = NA,
        `z value` = NA,
        `Pr(>|z|)` = 1,
        nobs = nrow(dat),
        flag = "failed",
        logs = log_msgs
      ))
    }
    
    result <- as.list(c(summary(m)$coefficients[2, ], nobs = nobs(m), flag = fl))
    result$logs <- log_msgs
    return(result)
    
  } else {
    # Linear regression case
    m <- tryCatch({
      withCallingHandlers({
        relmatLmer(
          as.formula(paste(outcome, "~", paste(c(exposure, covs, "(1|PIDN)"), collapse = " + "))),
          data = dat,
          relmat = list(PIDN = mat)
        )
      }, warning = function(w) {
        log_msgs <<- c(log_msgs, w$message)
        invokeRestart("muffleWarning")
      }, message = function(m) {
        log_msgs <<- c(log_msgs, m$message)
        invokeRestart("muffleMessage")
      })
    }, error = function(e) {
      log_msgs <<- c(log_msgs, sprintf("Error for %s: %s", exposure, e$message))
      NULL
    })
    
    if (is.null(m)) {
      return(list(
        Estimate = NA,
        `Std. Error` = NA,
        `t value` = NA,
        `Pr(>Chisq)` = 1,
        nobs = nrow(dat),
        flag = "failed",
        logs = log_msgs
      ))
    }
    
    fl <- "relmatLmer"
    
    # compare to baseline model and fetch result
    t <- suppressWarnings(suppressMessages(anova(m0, m)))
    a <- summary(m)
    
    result <- as.list(c(a$coefficients[2, ], t$`Pr(>Chisq)`[2], nobs = nobs(m), flag = fl))
    result$logs <- log_msgs
    return(result)
  }
}

#' Comprehensive biomarker performance analysis
#'
#' Performs a complete classification and diagnostic analysis for a continuous
#' biomarker against a binary outcome. Calculates ROC-AUC, optimal thresholds,
#' classification metrics, effect sizes, and statistical tests.
#'
#' @param dat Data.table or data.frame containing the data
#' @param outcome String specifying the column name of the binary outcome variable.
#'   Should have exactly two levels (e.g., "Q1" and "Q4", "Case" and "Control")
#' @param exposure String specifying the column name of the continuous exposure/
#'   biomarker variable
#'
#' @return A list containing biomarker performance metrics (see details)
#'
#' @details
#' Returns a list with ROC metrics, classification performance, effect sizes,
#' statistical tests, and group characteristics.
#'
#' @importFrom pROC roc auc ci.auc coords
#' @importFrom caret confusionMatrix
#' @importFrom stats wilcox.test var median sd
#' @importFrom rstatix cohens_d
#'
#' @export
analyze_biomarker <- function(dat, outcome, exposure) {
  # Load required packages
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required but not installed.")
  }
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' is required but not installed.")
  }
  if (!requireNamespace("rstatix", quietly = TRUE)) {
    stop("Package 'rstatix' is required but not installed.")
  }
  
  # Get outcome levels and ensure they're character
  outcome_levels <- sort(unique(as.character(dat[[outcome]])))
  if (length(outcome_levels) != 2) {
    stop("Outcome variable must have exactly 2 levels")
  }
  
  # 1. AUC-ROC
  roc_obj <- pROC::roc(
    response = dat[[outcome]], 
    predictor = dat[[exposure]], 
    levels = outcome_levels,
    quiet = TRUE
  )
  auc_val <- pROC::auc(roc_obj)
  ci_auc <- pROC::ci.auc(roc_obj)
  
  # 2. Optimal threshold (Youden's index)
  coords_obj <- pROC::coords(roc_obj, "best", best.method = "youden", 
                             transpose = FALSE)
  threshold <- coords_obj$threshold
  sensitivity <- coords_obj$sensitivity
  specificity <- coords_obj$specificity
  
  # 3. Classification at median threshold
  median_threshold <- median(dat[[exposure]], na.rm = TRUE)
  
  # Create complete cases dataset
  complete_idx <- !is.na(dat[[exposure]]) & !is.na(dat[[outcome]])
  dat_complete <- dat[complete_idx, ]
  
  # Create predictions and true values as factors with same levels
  pred_class <- factor(
    ifelse(dat_complete[[exposure]] > median_threshold, outcome_levels[2], outcome_levels[1]),
    levels = outcome_levels
  )
  true_class <- factor(as.character(dat_complete[[outcome]]), levels = outcome_levels)
  
  # Calculate confusion matrix
  cm <- caret::confusionMatrix(
    data = pred_class, 
    reference = true_class, 
    positive = as.character(outcome_levels[2])
  )
  
  # 4. Effect size: Cohen's d
  q1_vals <- dat[[exposure]][dat[[outcome]] == outcome_levels[1]]
  q4_vals <- dat[[exposure]][dat[[outcome]] == outcome_levels[2]]
  
  # Remove NAs
  q1_vals <- q1_vals[!is.na(q1_vals)]
  q4_vals <- q4_vals[!is.na(q4_vals)]
  
  cohens_d <- cohens_d(
    data = dat, 
    formula = as.formula(paste(exposure, "~", outcome)),
    ref.group = factors[2]
  )$effsize
  
  # 5. Wilcoxon test
  wilcox_test <- wilcox.test(dat[[exposure]] ~ dat[[outcome]])
  
  # Return as data.table with one row
  data.table(
    biomarker = exposure,
    # ROC metrics
    AUC = as.numeric(auc_val),
    AUC_CI_lower = as.numeric(ci_auc[1]),
    AUC_CI_upper = as.numeric(ci_auc[3]),
    # Optimal threshold
    optimal_threshold = threshold,
    sensitivity_optimal = sensitivity,
    specificity_optimal = specificity,
    # Median threshold performance
    accuracy = as.numeric(cm$overall["Accuracy"]),
    accuracy_CI_lower = as.numeric(cm$overall["AccuracyLower"]),
    accuracy_CI_upper = as.numeric(cm$overall["AccuracyUpper"]),
    PPV = as.numeric(cm$byClass["Pos Pred Value"]),
    NPV = as.numeric(cm$byClass["Neg Pred Value"]),
    sensitivity_median = as.numeric(cm$byClass["Sensitivity"]),
    specificity_median = as.numeric(cm$byClass["Specificity"]),
    # Effect size
    cohens_d = cohens_d,
    # Statistical test
    wilcox_p = wilcox_test$p.value,
    # Group statistics
    Q1_mean = mean(q1_vals),
    Q1_sd = sd(q1_vals),
    Q1_n = length(q1_vals),
    Q4_mean = mean(q4_vals),
    Q4_sd = sd(q4_vals),
    Q4_n = length(q4_vals)
  )
}



# Main -------------------------
flog.info("=== Setting up for analysis ===")
dat <- fread(input)
names(dat) <- make.names(names(dat))
names(prots) <- make.names(prots)
# similar to other name shielding, but reversed for ease of use
orig_outcome <- outcome
outcome <- make.names(outcome)

# add in single category columns with feature importance
# PlateId was important, but P0028449, P0028565, P0028437, P0028430, P0028452, P0028448 were top30 in either glmcross-validated or plr estimates
# dat = recode_dummies(dat, "PlateId", keep=c("P0028451", "P0028463", "P0028462", "P0028438", "P0028430", "P0028437"))
if(length(recodes) > 0) {
  new_covs <- c()
  for(recode_item in recodes) {
    if(recode_item[1] %in% covs) {
      dat <- recode_dummies(dat, 
        dcol = recode_item[1],
        keep = recode_item[2:length(recode_item)]
      )
      new_covs <- c(new_covs, paste(recode_item[1], recode_item[2:length(recode_item)], sep = "."))
      covs <- covs[covs != recode_item[1]]
    } else {
      flog.info(sprintf("%s not found in input, dropping & skipping recode...", recode_item[1]))
      covs <- covs[covs != recode_item[1]]
    }
  }
  covs <- c(covs, new_covs)
}
# list covariates
flog.info(paste("Running covs:", paste(covs, collapse = ", ")))

# read in and generate the relatedness matrix
flog.info("Computing relatedness matrix...")
pcr <- readRDS(file.path(ref))

# Capture output (both stdout and stderr)
all_output <- capture.output({
  grm <- pcrelateToMatrix(pcr, 
    sample.include = as.character(unique(dat$PIDN)), 
    thresh = 2^(-11/2), 
    scaleKin = 2
  )
}, type = "message")

# Log it
if (length(all_output) > 0) {
  flog.info("pcrelateToMatrix output:\n%s", paste(all_output, collapse = "\n"))
}
flog.info("Relatedness matrix dimensions: %d x %d", nrow(grm), ncol(grm))

# make output table
dt1 <- CJ(dx = outcome, ix = names(prots))

# generate the baseline model w/o proteins
if (binomial) {
  # If outcome is character, convert to factor:
  if(is.character(dat[[outcome]]) & length(unique(dat[[outcome]])) == 2) {
    unique_vals <- unique(dat[[outcome]])
    if (
      length(unique_vals) != 2 
      || !all(unique_vals %in% factors) 
      || !all(factors %in% unique_vals)
    ) {
      flog.error("Invalid factors specified for '%s'. Provided: %s, Available: %s", 
        outcome, 
        paste(factors, collapse = ", "), 
        paste(unique_vals, collapse = ", ")
      )
      stop("exactly two factors must be specified for binomial models")
    }
    dat[[outcome]] <- factor(dat[[outcome]], levels = factors)
  }
  
  # scale covariates for logistic regression
  dat[, (scale_cov) := lapply(.SD, function(x) scale(x)[,1]), .SDcols = scale_cov]
  cols <- c("beta", "se", "z stat", "pval", "nobs", "flag")
} else {
  # scale covariates 
  if (length(scale_cov > 0)) {
    dat[, (scale_cov) := lapply(.SD, function(x) scale(x)[,1]), .SDcols = scale_cov]
  }
  # baseline model
  cols <- c("beta", "se", "z stat", "pval", "nobs", "flag")
  m0 <- relmatLmer(
    as.formula(paste(outcome, "~", paste(c(covs, "(1|PIDN)"), collapse = " + "))),
    data = dat,
    relmat = list(PIDN = grm)
  )

  cols <- c("beta", "se", "t stat", "pval", "nobs", "flag")
}


elapsed <- toc(log = TRUE, quiet = TRUE)
flog.info("=== Analysis started: %s ===", elapsed$callback_msg)
flog.info(
  "Outcome: %s, Proteins: %d, Cores: %d", outcome, nrow(dt1), threads
)
# fast run function
results <- future_lapply(1:nrow(dt1), function(i) {
  get_linreg(
    dat = dat,
    m0 = if(binomial) NULL else m0,
    outcome = dt1$dx[i],
    exposure = dt1$ix[i],
    covs = covs,
    mat = grm,
    logistic = binomial
  )
}, future.seed = TRUE, future.stdout = FALSE, future.conditions = character(0))

# Extract all collected log messages and write them to the log file
all_logs <- lapply(seq_along(results), function(i) {
  if (!is.null(results[[i]]$logs) && length(results[[i]]$logs) > 0) {
    # Add context about which iteration
    sprintf("[Row %d: %s vs %s] %s", 
            i, 
            dt1$dx[i], 
            dt1$ix[i], 
            results[[i]]$logs)
  }
})

# Flatten and write to log
all_logs <- unlist(all_logs)
if (length(all_logs) > 0) {
  for (msg in all_logs) {
    flog.info(msg)
  }
}

# Remove logs from results to clean up the output
results <- lapply(results, function(r) {
  r$logs <- NULL
  r
})

# Convert results to data.table and bind
results_dt <- rbindlist(results)
dt1[, (cols) := results_dt]


if (binomial) {
  # adding cohens_d for comparison
  dt1[flag != "failed", cohens_d := cohens_d(
    data = dat, 
    formula = as.formula(paste(ix, "~", dx)),
    ref.group = factors[2]
  )$effsize, by = .I]

  # Generate a summary report for quasi-complete separation subset
  flog.info("=== Analyzing quasi-complete separation genes ===")
  b_res <- rbindlist(lapply(dt1[flag == "separated", ix], function(k) {
    analyze_biomarker(dat, outcome, exposure = k)
  }))

  # Remove Xs needed for formula
  b_res$biomarker <- plyr::mapvalues(
    x = b_res$biomarker, 
    from = names(prots), 
    to = prots,
    warn_missing = FALSE
  )
  fwrite(b_res, 
    file = file.path(output, paste(nameset, "separation_genes", "csv.gz", sep = "."))
  )
}

# Summary
# remove Xs needed for formula
dt1$ix <- plyr::mapvalues(x = dt1$ix, from = names(prots), to = prots)
dt1$dx <- plyr::mapvalues(x = dt1$dx, from = outcome, to = orig_outcome)
# dt1[startsWith(ix, "X."), ix := gsub("^X.", ".", ix)]
fwrite(dt1, 
  file = file.path(output, paste(nameset, "lme4qtl", "csv.gz", sep = "."))
)

flog.info("Summary: Analyzed %d biomarkers", nrow(dt1))
invisible({
  dt1[, .N, by = flag][order(-N)] |> 
    apply(1, function(x) flog.info("  %s: %s", x["flag"], x["N"]))
})

# cleanup
end <- toc(log = TRUE, quiet = TRUE)
flog.info("=== Total run time: %s ===", end$callback_msg)


# # model testing ---------------------
# # spot check 10 proteins
# set(123)
# prots = sample(names(dat)[grep(first_prot, names(dat))[1]:last(grep(last_prot, names(dat)))], 10)

# model = glm(paste("PlasmaPTau181", "~", paste(c(prots[2], cols), collapse=" + ")),
# na.action=na.exclude, data=dat)
# summary(model)
# plot(fitted(model), resid(model, type = "pearson"))
# abline(0,0, col="red")
# shapiro.test(resid(model, type="pearson"))
