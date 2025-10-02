#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 08-Cell-type_scoring.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Score celltype for each sample using weighted GSEA method
#'  
#'
#' Author: Robert B Butler III
#'
#' Date Created: 2025-09-30
#'
#' Copyright (c) 2025
#' Email: rrbutler@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'   
#'   
#'   
#'   
#'
#' ---------------------------

#' load up the packages we will need:  (uncomment as required)
library(data.table)
library(optparse)
library(ggplot2)
library(patchwork)
library(future.apply)
library(tictoc)

# Nameset helper function
nameset_default <- function(input){
  file_path_sans_ext(basename(input))
}

# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Biomarker results to use for probe filtering, keeps most significant"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project [default: Nameset of input]"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "celltype_scores",
    help = "File path to output directory [default: %default]"
  ),
  make_option(
    c("-r", "--ref"),
    type = "character",
    default = "cell_type_prot/Sharma_etal.EWCE.human.lvl1.spec_matr.txt",
    help = "Path to specificity matrix [default: %default]"
   ),
  make_option(
    c("-t", "--threads"),
    default = 0,
    help = "num threads [default: parallelly::availableCores()]"
   ),
  make_option(
    c("-p", "--plot"), 
    action = "store_true", 
    default = FALSE,
    help = "Plot a random sampling of 10 individuals in output folder"
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
  "output",
  "ref",
  "threads",
  "nameset",
  "plot"
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

# check the current active plan
if (threads == 0) {
  threads <- as.integer(parallelly::availableCores()) - 2
} else if (is.numeric(threads)) {
  threads <- as.integer(threads)
} else {
  cat("threads must be an integer value")
  stop(parse_args(OptionParser(option_list = option_list), args = "--help"),
    call. = FALSE
  )
}
print(sprintf("Setting up %d forks for samples", threads))
setDTthreads(1)
plan('multicore', workers = threads)

# input variables
input <- strsplit(input, ",")[[1]]

# check outpath status
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}


# Functions -------------------------
#' Calculate Cell Type Enrichment Using Specificity-Weighted GSEA
#'
#' This function implements a weighted Gene Set Enrichment Analysis (GSEA) approach
#' for proteomic data using cell type specificity scores as weights. It calculates
#' enrichment scores for each cell type based on protein expression rankings and
#' their cell type specificity.
#'
#' @param expr_matrix Numeric matrix of protein expression values. Rows are proteins,
#'        columns are samples.
#' @param specificity_matrix Numeric matrix of cell type specificity scores (0-1 scale).
#'        Rows are proteins, columns are cell types.
#' @param sample_col Integer or character. Which column of expr_matrix to use for ranking.
#'        Default is 1.
#' @param specificity_threshold Numeric or character. Threshold for considering proteins as "specific"
#'        to a cell type. Can be a numeric value (0-1) or "median" to use the median of non-zero
#'        specificity values for each cell type. Default is 0 (no threshold).
#' @param nperm Integer. Number of permutations for significance testing. Default is median.
#'        Set to 0 to skip permutation testing.
#' @param p_weight Numeric. Exponent for the ranking metric in the enrichment score 
#'        calculation. Default is 1.
#' @param min_size Integer. Minimum number of proteins with non-zero specificity to 
#'        consider a cell type for analysis. Default is 5.
#' @param verbose Logical. Whether to print progress messages. Default is TRUE.
#'
#' @return A list containing:
#'   \item{enrichment_scores}{Numeric vector of enrichment scores for each cell type}
#'   \item{normalized_scores}{Numeric vector of normalized enrichment scores}
#'   \item{p_values}{P-values from permutation testing, if nperm > 0}
#'   \item{fdr}{False discovery rates, adjusted p-values}
#'   \item{leading_edge}{List of leading edge proteins for each significant cell type}
#'   \item{running_scores}{Matrix of running enrichment score values for plotting}
#'   \item{thresholds_used}{Numeric vector of specificity thresholds used for each cell type}
#'
#' @examples
#' # Create example data
#' proteins <- paste0("protein", 1:100)
#' cell_types <- paste0("celltype", 1:5)
#' samples <- paste0("sample", 1:3)
#' 
#' # Example expression matrix (proteins x samples)
#' expr_data <- matrix(rnorm(100*3), nrow=100, ncol=3)
#' rownames(expr_data) <- proteins
#' colnames(expr_data) <- samples
#' 
#' # Example specificity matrix (proteins x cell types)
#' spec_data <- matrix(runif(100*5), nrow=100, ncol=5)
#' rownames(spec_data) <- proteins
#' colnames(spec_data) <- cell_types
#' 
#' # Run the analysis
#' results <- specificityWeightedGSEA(expr_data, spec_data, sample_col=1, nperm=100)
#' 
#' # Plot top enriched cell type
#' plotEnrichment(results, cell_type=results$enrichment_scores[1])
#'
#' @references
#' Adapted from the Gene Set Enrichment Analysis methodology:
#' Subramanian, A., et al. (2005). Gene set enrichment analysis: A knowledge-based 
#' approach for interpreting genome-wide expression profiles. PNAS, 102(43), 15545-15550.
#'
#' @export

specificityWeightedGSEA <- function(expr_matrix, 
                                    specificity_matrix, 
                                    sample_col = 1,
                                    specificity_threshold = "median", 
                                    nperm = 1000, 
                                    p_weight = 1, 
                                    min_size = 5, 
                                    verbose = TRUE
                                    ) {
  # Input validation
  if (!is.matrix(expr_matrix) || !is.matrix(specificity_matrix)) {
    stop("Expression and specificity inputs must be matrices")
  }
  
  if (is.character(sample_col)) {
    if (!sample_col %in% colnames(expr_matrix)) {
      stop("Sample column name not found in expression matrix")
    }
    sample_col <- which(colnames(expr_matrix) == sample_col)
  }
  
  # Check for matching protein names
  if (is.null(rownames(expr_matrix)) || is.null(rownames(specificity_matrix))) {
    stop("Both matrices must have rownames identifying proteins")
  }
  
  # Find common proteins
  common_proteins <- intersect(rownames(expr_matrix), rownames(specificity_matrix))
  if (length(common_proteins) == 0) {
    stop("No common proteins found between expression and specificity matrices")
  }
  if (verbose) {
    message(paste("Found", length(common_proteins), "common proteins for analysis"))
  }
  
  # Subset matrices to common proteins
  expr_matrix <- expr_matrix[common_proteins, , drop = FALSE]
  specificity_matrix <- specificity_matrix[common_proteins, , drop = FALSE]
  
  # Sort proteins by rank
  ranked_values <- expr_matrix[, sample_col]
  protein_ranks <- sort(ranked_values, decreasing = TRUE, index.return = TRUE)
  ranked_proteins <- rownames(expr_matrix)[protein_ranks$ix]
  
  # Calculate enrichment scores for each cell type
  cell_types <- colnames(specificity_matrix)
  n_proteins <- length(ranked_proteins)
  
  # Initialize results containers
  enrichment_scores <- numeric(length(cell_types))
  names(enrichment_scores) <- cell_types
  running_scores <- matrix(0, nrow = n_proteins, ncol = length(cell_types))
  rownames(running_scores) <- ranked_proteins
  colnames(running_scores) <- cell_types
  leading_edge <- list()
  thresholds_used <- numeric(length(cell_types))
  names(thresholds_used) <- cell_types
  
  # Calculate enrichment score for each cell type
  for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]
    if (verbose && (i %% 10 == 0 || i == 1)) {
      message(paste("Processing cell type:", i, "of", length(cell_types)))
    }
    
    # Get specificity scores for this cell type
    spec_scores <- specificity_matrix[, cell_type]
    
    # Determine threshold for this cell type
    if (is.character(specificity_threshold) && specificity_threshold == "median") {
      # Use median of non-zero specificity values
      non_zero_specs <- spec_scores[spec_scores > 0]
      if (length(non_zero_specs) > 0) {
        actual_threshold <- median(non_zero_specs)
      } else {
        actual_threshold <- 0
      }
      if (verbose) {
        message(paste("Using median threshold for", cell_type, ":", actual_threshold))
      }
    } else {
      # Use provided numeric threshold
      actual_threshold <- specificity_threshold
    }
    
    # Store the threshold used
    thresholds_used[i] <- actual_threshold
    
    # Skip cell types with too few specific proteins
    if (sum(spec_scores > 0) < min_size) {
      if (verbose) {
        message(paste("Skipping", cell_type, "- too few specific proteins"))
      }
      enrichment_scores[i] <- NA
      next
    }
    
    # Calculate hit_indices (where each protein is in the ranked list)
    hit_indices <- match(names(spec_scores), ranked_proteins)
    
    # Calculate running sum
    running_sum <- 0
    max_deviation <- 0
    max_position <- 0
    
    # Calculate sum of weights for normalization (only for proteins above threshold)
    N_R <- sum((spec_scores[spec_scores > actual_threshold])^p_weight)
    
    # Count proteins below threshold for penalty calculation
    N_non_specific <- sum(spec_scores <= actual_threshold)
    
    # Walk down the ranked list
    for (j in 1:n_proteins) {
      protein <- ranked_proteins[j]
      
      # If protein has specificity for this cell type, add weighted score
      if (spec_scores[protein] > actual_threshold) {
        # Scale the positive steps
        running_sum <- running_sum + (spec_scores[protein]^p_weight / N_R)
      } else {
        # All other proteins count as "misses"
        running_sum <- running_sum - (1 / N_non_specific)
      }
      
      # Record running sum
      running_scores[j, i] <- running_sum
      
      # Check if this is the maximum deviation
      if (abs(running_sum) > abs(max_deviation)) {
        max_deviation <- running_sum
        max_position <- j
      }
    }
    
    # Store enrichment score
    enrichment_scores[i] <- max_deviation
    
    # Store leading edge proteins if enrichment score is positive
    if (max_deviation > 0) {
      leading_edge[[cell_type]] <- ranked_proteins[1:max_position]
    } else {
      leading_edge[[cell_type]] <- ranked_proteins[max_position:n_proteins]
    }
  }
  
  # Perform permutation testing if requested
  p_values <- rep(NA, length(cell_types))
  names(p_values) <- cell_types
  normalized_scores <- enrichment_scores
  
  if (nperm > 0) {
    if (verbose) {
      message("Performing permutation testing...")
    }
    
    # Calculate null distribution by permuting protein rankings
    null_scores <- matrix(0, nrow = nperm, ncol = length(cell_types))
    
    for (p in 1:nperm) {
      if (verbose && (p %% 100 == 0)) {
        message(paste("Permutation:", p, "of", nperm))
      }
      
      # Permute protein ranks
      perm_ranks <- sample(n_proteins)
      perm_proteins <- ranked_proteins[perm_ranks]
      
      # Calculate enrichment score for each cell type using permuted ranks
      for (i in seq_along(cell_types)) {
        if (is.na(enrichment_scores[i])) next
        
        cell_type <- cell_types[i]
        spec_scores <- specificity_matrix[, cell_type]
        actual_threshold <- thresholds_used[i]
        
        # Calculate running sum for permuted ranks
        running_sum <- 0
        max_deviation <- 0
        
        # Calculate sum of weights for normalization (using same threshold)
        N_R <- sum((spec_scores[spec_scores > actual_threshold])^p_weight)
        N_non_specific <- sum(spec_scores <= actual_threshold)
        
        for (j in 1:n_proteins) {
          protein <- perm_proteins[j]
          
          if (spec_scores[protein] > actual_threshold) {
            running_sum <- running_sum + (spec_scores[protein]^p_weight / N_R)
          } else {
            running_sum <- running_sum - (1 / N_non_specific)
          }
          
          if (abs(running_sum) > abs(max_deviation)) {
            max_deviation <- running_sum
          }
        }
        
        null_scores[p, i] <- max_deviation
      }
    }
    
    # Calculate p-values and normalize enrichment scores
    for (i in seq_along(cell_types)) {
      if (is.na(enrichment_scores[i])) next
      
      # Calculate p-value with pseudocount
      if (enrichment_scores[i] >= 0) {
        p_values[i] <- (sum(null_scores[, i] >= enrichment_scores[i]) + 1) / (nperm + 1)
      } else {
        p_values[i] <- (sum(null_scores[, i] <= enrichment_scores[i]) + 1) / (nperm + 1)
      }
      
      # Normalize enrichment score
      null_mean <- mean(null_scores[, i])
      null_sd <- sd(null_scores[, i])
      if (null_sd > 0) {
        normalized_scores[i] <- (enrichment_scores[i] - null_mean) / null_sd
      }
    }
    
    # Adjust p-values for multiple testing
    fdr <- p.adjust(p_values, method = "BH")
  } else {
    fdr <- rep(NA, length(cell_types))
    names(fdr) <- cell_types
  }
  
  # Return results
  results <- list(
    enrichment_scores = enrichment_scores,
    normalized_scores = normalized_scores,
    p_values = p_values,
    fdr = fdr,
    leading_edge = leading_edge,
    running_scores = running_scores,
    cell_types = cell_types,
    ranked_proteins = ranked_proteins,
    thresholds_used = thresholds_used
  )
  
  return(results)
}


#' Plot Cell Type Enrichment Results
#'
#' Creates a GSEA-style enrichment plot for a specified cell type from
#' specificityWeightedGSEA results.
#'
#' @param gsea_results List. Results from specificityWeightedGSEA function.
#' @param cell_type Character. Name of cell type to plot.
#' @param show_pval Logical. Whether to display p-value on plot. Default is TRUE.
#' @param title Character. Optional title for the plot.
#'
#' @return A ggplot2 object with the enrichment plot.
#'
#' @import ggplot2
#' @export
plotEnrichment <- function(gsea_results, cell_type, show_pval = TRUE, 
                          title = NULL) {
  # Check if required packages are installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is needed for this function")
  }
  
  # Extract from data.table if not list
  if (is.data.table(gsea_results)) {
    gsea_results <- unlist(gsea_results, recursive = FALSE)
  }
  
  # Input validation
  if (!cell_type %in% gsea_results$cell_types) {
    stop("Cell type not found in results")
  }
  
  # Get cell type index
  ct_idx <- which(gsea_results$cell_types == cell_type)
  
  # Extract running scores
  run_scores <- gsea_results$running_scores[, ct_idx]
  
  # Create data frame for plotting
  plot_data <- data.table(
    position = 1:length(run_scores),
    running_score = run_scores
  )
  
  # Create annotation text
  if (show_pval) {
    pval <- gsea_results$p_values[ct_idx]
    pval_text <- paste("p-value =", format(pval, digits = 3))
    es <- gsea_results$enrichment_scores[ct_idx]
    nes <- gsea_results$normalized_scores[ct_idx]
    es_text <- paste("ES =", format(es, digits = 3), 
                    "NES =", format(nes, digits = 3))
  }
  
  # Set plot title
  if (is.null(title)) {
    plot_title <- paste("Enrichment plot:", cell_type)
  } else {
    plot_title <- title
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = position, y = running_score)) +
    ggplot2::geom_line(color = "#40004b", size = 1) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(
      xintercept = plot_data[which.max(running_score), position], 
      linetype = "dashed",
      color = "#1b7837"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = plot_title,
      x = "Protein Rank",
      y = "Enrichment Score"
    )
  
  # Add p-value annotation
  if (show_pval) {
    p <- p + ggplot2::annotate(
      "text", 
      x = max(plot_data$position) * 0.8, 
      y = max(abs(plot_data$running_score)) * 0.8,
      label = paste(pval_text, "\n", es_text)
    )
  }
  
  return(p)
}


#' Print multiple plots across pages in a grid layout
#'
#' This function takes a list of ggplot2 objects and arranges them in a 
#' specified grid layout per page, then prints them to a PDF file. By default, 
#' plots are arranged in a 2x2 grid (4 plots per page).
#'
#' @param plot_list A list of ggplot2 objects to be arranged and printed
#' @param plots_per_page Integer specifying number of plots per page (default: 4)
#' @param file Character string for the output PDF filename (default: "plots_by_page.pdf")
#' @param width PDF width in inches (default: 10)
#' @param height PDF height in inches (default: 8)
#' @param title Optional overall title for each page (default: NULL)
#' @param tag_levels Whether to add tags (labels) to plots. Can be 'A' for uppercase letters,
#'        'a' for lowercase letters, '1' for numbers, or NULL for no tags (default: NULL)
#'
#' @return Invisibly returns TRUE if successful
#'
#' @details
#' Currently only supports layouts with 4 plots per page in a 2x2 grid.
#' Empty plots are added automatically if the number of plots doesn't fill 
#' the last page completely.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' 
#' # Create a list of 10 plots
#' plot_list <- list()
#' for (i in 1:10) {
#'   plot_list[[i]] <- ggplot(mtcars, aes(hp, mpg)) + 
#'     geom_point() + 
#'     ggtitle(paste("Plot", i))
#' }
#' 
#' # Print plots 4 per page
#' print_plots_per_page(plot_list)
#' 
#' # With custom file name and page labels
#' print_plots_per_page(plot_list, 
#'                      file = "my_analysis_results.pdf", 
#'                      tag_levels = 'A')
#' }
#'
#' @import ggplot2
#' @import patchwork
#' @export
print_plots_per_page <- function(plot_list, 
                                 plots_per_page = 4, 
                                 file = "plots_by_page.pdf",
                                 width = 11, 
                                 height = 8.5,
                                 title = NULL,
                                 tag_levels = NULL) {
  
  # Check inputs
  if (!inherits(plot_list, "list")) {
    stop("plot_list must be a list of ggplot objects")
  }
  
  if (plots_per_page != 4) {
    stop("Currently only supports 4 plots per page")
  }
  
  if (length(plot_list) == 0) {
    warning("Empty plot list provided, no output generated")
    return(invisible(FALSE))
  }
  
  # Calculate pages needed
  n_plots <- length(plot_list)
  n_pages <- ceiling(n_plots / plots_per_page)
  
  # Open PDF device
  pdf(file, width = width, height = height)
  
  # For each page
  for (i in 1:n_pages) {
    # Get plot indices for this page
    start_idx <- (i-1) * plots_per_page + 1
    end_idx <- min(i * plots_per_page, n_plots)
    
    # Select plots for current page
    current_plots <- plot_list[start_idx:end_idx]
    
    # Fill with empty plots if needed to complete the grid
    if (length(current_plots) < plots_per_page) {
      empty_plots <- plots_per_page - length(current_plots)
      for (j in 1:empty_plots) {
        current_plots <- c(current_plots, list(ggplot() + theme_void()))
      }
    }
    
    # Create a 2x2 grid for this page
    layout <- (current_plots[[1]] + current_plots[[2]]) / 
              (current_plots[[3]] + current_plots[[4]])
    
    # Add annotations if requested
    if (!is.null(title) || !is.null(tag_levels)) {
      layout <- layout + patchwork::plot_annotation(
        title = title,
        tag_levels = tag_levels
      )
    }
    
    # Print the current page
    print(layout)
  }
  
  # Close the device
  dev.off()
  
  # Return success invisibly
  return(invisible(TRUE))
}


# Main -----------------------------------------------------------------------

# Stage data for analysis --------------------------
# DE file locations
fh <- file.path("baseline_PCs_relatedness", input,
  paste("dataProt_SS-205063.PC.related", input, "csv", sep = ".")
)
result <- lapply(fh, fread)
setattr(result, 'names', input)

# merge
result <- rbindlist(result)
input <- gsub(".APOE", "", input)

# filtered set of genes with at least one significant result, top sig duplicate probe
all_genes <- unique(result[order(padj), .(Symbol, ix, padj)], by = 'Symbol')
# bio_genes <- unique(all_genes[padj < 0.05, .(Symbol, ix)], by = 'Symbol')

# get specificity index
spec <- fread(ref)
cols <- names(spec)[-1]
# spec <- as.matrix(spec[!GENE %in% bio_genes$Symbol, (cols) := 0], rownames = 'GENE')
spec <- as.matrix(spec, rownames = 'GENE')

# get expression data and filter, rename, flip
prots <- fread("dataProt_SS-205063.merged.QC.csv")
cols <- c("SampleID", all_genes$ix)
prots <- prots[, ..cols]
setnames(prots, old = all_genes$ix, new = all_genes$Symbol)
prots <- transpose(prots, keep.names = 'Symbol', make.names = 'SampleID')
# prots <- as.matrix(prots[Symbol %in% bio_genes$Symbol], rownames = 'Symbol')
prots <- as.matrix(prots, rownames = 'Symbol')


# Build result table(s) & run function -------------
# Create grouping factor for columns
groups <- ceiling(seq_len(ncol(prots)) / ceiling(ncol(prots) / threads))
dt_list <- lapply(seq(threads), function(i) {
  data.table('SampleID' = colnames(prots)[groups == i])
})
dt_list <- future_lapply(seq(dt_list), function(i) {
  tic()
  x <- dt_list[[i]][, {
    # Run the analysis
    results <- specificityWeightedGSEA(
      expr_matrix = prots,
      specificity_matrix = spec,
      sample_col = as.character(SampleID),
      nperm = 1000, 
      verbose = FALSE
    )
    
    # Extract each component of results as a separate column
    # Note: for list/complex columns, we'll convert to JSON strings or serialize as needed
    list(
      # Scalar or simple vector results
      cell_types = list(results$cell_types),
      enrichment_scores = list(results$enrichment_scores),
      normalized_scores = list(results$normalized_scores),
      p_values = list(results$p_values),
      fdr_values = list(results$fdr),
      thresholds_used = list(results$thresholds_used),
      
      # For the more complex elements, we can either:
      # 1. Extract specific information of interest
      top_cell_type = names(results$enrichment_scores)[which.max(results$enrichment_scores)],
      top_score = max(results$enrichment_scores, na.rm = TRUE),
      significant_cell_types = list(names(which(results$fdr < 0.05))),
      
      # 2. Store as serialized/JSON strings
      leading_edge = list(jsonlite::toJSON(results$leading_edge)),
      running_scores = list(results$running_scores)
    )
  }, by = 'SampleID']
  print(sprintf("Group %s: %s", i, toc()$callback_msg))
  return(x)
}, future.seed = TRUE)


# join and save files ------------------------------
dt <- rbindlist(dt_list)
saveRDS(dt, file = file.path(output, paste0(nameset, ".scores.rds")))

if (plot) {
  # Set seed for reproducible sampling
  set.seed(123)
  sampled_dt <- sample(dt[["SampleID"]], size = 10)
  lapply(sampled_dt, function(i) {
    plist <- lapply(unlist(dt[SampleID == i, cell_types]), function(j) {
      plotEnrichment(dt[SampleID == i], j, title = paste(j, "in", i))
    })
    print_plots_per_page(plist, 
      file = file.path(output, paste(nameset, i, "plots.pdf", sep = "."))
    )
  })
}

