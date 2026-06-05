#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 08-String_module_construction.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Extract STRING network modules and prepar for Scoring with weighted GSEA
#'  
#'
#' Author: Robert B Butler III
#'
#' Date Created: 2025-11-05
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
library(ggplot2)
library(ggridges)
library(patchwork)
library(gtools)
library(igraph)
library(stringr)


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
    default = "module_scores",
    help = "File path to output directory [default: %default]"
  ),
  make_option(
    c("-r", "--ref"),
    type = "character",
    default = "manuscript/Figure_2/dataProt_SS-205063.PC.related.all_bios.overlap_enrich_clusters.csv",
    help = "Path to overlap enrichment table [default: %default]"
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
  "prots",
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

# input variables
input <- strsplit(input, ",")[[1]]

# check outpath status
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}


# Functions -------------------------
#' Calculate Node Centrality Metrics by Cluster
#'
#' This function builds a graph from edge and vertex data, subsets it by cluster,
#' and calculates various centrality metrics (degree, strength, and alpha centrality
#' with different alpha values) for nodes within each cluster.
#'
#' @param edges_dt A data.frame or data.table containing edge information for graph construction.
#'   Must be compatible with \code{\link[igraph]{graph_from_data_frame}}. Should contain
#'   a column specified by \code{score_col} that will be used as edge weights.
#' @param exo_dt A data.frame or data.table containing vertex information with vertex attributes.
#'   Must include a column specified by \code{exo_col} (values used for alpha centrality calculations)
#'   and a cluster column specified by the \code{clus} parameter.
#' @param clus Character string specifying the name of the column in \code{exo_dt} that
#'   contains cluster assignments. Default is "cluster".
#' @param score_col Character string specifying the name of the column in \code{edges_dt}
#'   that contains edge weights/scores. Default is "stringdb.score".
#' @param exo_col Character string specifying the name of the column in \code{exo_dt}
#'   that contains exogenous values for alpha centrality calculation. Default is "abs.beta".
#'
#' @return A data.table with the following columns:
#'   \describe{
#'     \item{clus}{Cluster identifier (column name matches the \code{clus} parameter)}
#'     \item{node}{Node/vertex name}
#'     \item{degree}{Node degree (number of connections)}
#'     \item{strength}{Node strength (sum of edge weights)}
#'     \item{alpha_0.1}{Alpha centrality with alpha = 0.1}
#'     \item{alpha_0.5}{Alpha centrality with alpha = 0.5}
#'     \item{alpha_1}{Alpha centrality with alpha = 1.0}
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Handles missing values in the exogenous column by replacing NAs with
#'         floor(min(exo_col) * 10) / 10
#'   \item Constructs an igraph object using edges_dt and exo_dt
#'   \item Sets edge weights from the specified score column
#'   \item For each unique cluster, creates an induced subgraph
#'   \item Calculates centrality metrics for nodes within each cluster subgraph
#'   \item Combines results across all clusters into a single data.table
#' }
#'
#' Alpha centrality is calculated using the exogenous column values as exogenous factors,
#' with three different alpha parameters (0.1, 0.2, and 0.4) to capture different
#' aspects of node importance.
#'
#' @note
#' \itemize{
#'   \item Requires the igraph and data.table packages
#'   \item The specified score column in edges_dt is used as edge weights
#'   \item Missing values in the exogenous column are handled by imputation to prevent calculation errors
#'   \item The function assumes that vertex names in edges_dt match those in exo_dt
#' }
#'
#' @examples
#' \dontrun{
#' # Create example data
#' edges_dt <- data.table(
#'   from = c("A", "B", "C"),
#'   to = c("B", "C", "A"),
#'   stringdb.score = c(0.8, 0.6, 0.9),
#'   correlation = c(0.7, 0.5, 0.8)
#' )
#' 
#' exo_dt <- data.table(
#'   name = c("A", "B", "C"),
#'   abs.beta = c(1.2, 0.8, 1.5),
#'   log_pvalue = c(2.3, 1.8, 2.7),
#'   cluster = c(1, 1, 2)
#' )
#' 
#' # Calculate node metrics by cluster (default columns)
#' result1 <- node_metrics(edges_dt, exo_dt, clus = "cluster")
#' 
#' # Using different score and exogenous columns
#' result2 <- node_metrics(edges_dt, exo_dt, 
#'                        clus = "cluster",
#'                        score_col = "correlation", 
#'                        exo_col = "log_pvalue")
#' 
#' # Using a different cluster column
#' exo_dt$group <- c("X", "X", "Y")
#' result3 <- node_metrics(edges_dt, exo_dt, 
#'                        clus = "group",
#'                        score_col = "stringdb.score",
#'                        exo_col = "abs.beta")
#' }
#'
#' @seealso
#' \code{\link[igraph]{graph_from_data_frame}}, \code{\link[igraph]{induced_subgraph}},
#' \code{\link[igraph]{degree}}, \code{\link[igraph]{strength}}, \code{\link[igraph]{alpha_centrality}}
#'
#' @importFrom igraph graph_from_data_frame induced_subgraph degree strength alpha_centrality vertex_attr V E
#' @importFrom data.table data.table rbindlist setattr
#'
#' @export
node_metrics <- function(edges_dt, exo_dt, clus = "cluster", 
                        score_col = "stringdb.score", exo_col = "abs.beta") {
  
  # Input validation
  if (!score_col %in% names(edges_dt)) {
    stop(paste("Column", score_col, "not found in edges_dt"))
  }
  if (!exo_col %in% names(exo_dt)) {
    stop(paste("Column", exo_col, "not found in exo_dt"))
  }
  if (!clus %in% names(exo_dt)) {
    stop(paste("Column", clus, "not found in exo_dt"))
  }
  
  # prevent NAs in exogenous column
  exo_dt[is.na(abs.beta)] <- floor(min(exo_dt$abs.beta)*10)/10
  
  # build graph
  g <- graph_from_data_frame(
    d = edges_dt,
    vertices = exo_dt
  )
  
  # Setting weight column for current version of igraph. 
  # In the future, can directly use score_col in alpha_centrality(g, weights = score_col) 
  # https://github.com/igraph/rigraph/issues/915
  E(g)$weight <- edge_attr(g, score_col)
  
  # subset by clus column
  x <- lapply(unique(exo_dt[[clus]]), function(i) {
    g2 <- induced_subgraph(g, V(g)[vertex_attr(g, clus) == i])
    
    # calculate scores
    data.table(
      node = V(g2)$name, 
      degree = degree(g2), 
      strength = strength(g2),
      alpha_0.1 = alpha_centrality(g2, 
        alpha = 0.1, 
        exo = vertex_attr(g2, exo_col)
      ),
      alpha_0.2 = alpha_centrality(g2, 
        alpha = 0.2, 
        exo = vertex_attr(g2, exo_col)
      ),
      alpha_0.4 = alpha_centrality(g2, 
        alpha = 0.4, 
        exo = vertex_attr(g2, exo_col)
      )
    )
  })
  setattr(x, 'names', unique(exo_dt[[clus]]))
  x <- rbindlist(x, idcol = clus)
  
  return(x)
}

#' Joy plot of specificity matrix
#'
#' @param x, data.table with specificity matrix and gene name column
#' @param gen_col, string of gene name column
#' @param outname, string of output prefix to use
#'
#' @return None. Saves a joyplot in pdf format
spec_joy <- function(x, gen_col = "GENE", outname = "ewce_spec"){
  # prepping data (melt and label strip)
  ggmelt <- melt(x, id.vars = gen_col)
  ggmelt[ggmelt$value == 0, "value"] <- 5e-08 # log scale cannot be 0

  # plotting
  p1 <- ggplot(
      ggmelt,
      aes(x = value, y = variable, fill = factor(after_stat(quantile)))
  ) +
    stat_density_ridges(
      geom = "density_ridges_gradient",
      calc_ecdf = TRUE,
      quantiles = 10
    ) +
    scale_fill_brewer(type = "div", palette = "RdYlBu", direction = -1) +
    scale_x_log10(limits = c(5e-03, 1)) +
    labs(x = "Specificity", y = "Cell type", fill = "Decile") +
    theme_classic()

  ggsave(
    filename = paste(outname, "joy_plot.pdf", sep = "."),
    plot = p1,
    device = cairo_pdf
  )
}

#' Scale numeric vector to range [0,1]
#'
#' This function performs min-max normalization, scaling all values in a numeric vector
#' to the range [0,1]. The transformation is done using the formula (x-min(x))/(max(x)-min(x)).
#'
#' @param x A numeric vector to be scaled to the range [0,1]
#'
#' @return A numeric vector of the same length as \code{x} with values scaled to the range [0,1]
#'   where the minimum value in the original vector is mapped to 0 and the maximum value is mapped to 1.
#'   If all values in \code{x} are equal (i.e., min(x) = max(x)), the function returns a vector of the same
#'   length containing all 0s to avoid division by zero.
#'
#' @examples
#' range01(c(1, 2, 3, 4, 5))  # Returns (0, 0.25, 0.5, 0.75, 1)
#' range01(c(0, 50, 100))     # Returns (0, 0.5, 1)
#' range01(c(-10, 0, 10))     # Returns (0, 0.5, 1)
#' 
#' # Handling edge case
#' range01(rep(5, 3))         # Returns (0, 0, 0) when all values are identical
#'
#' @seealso \code{\link{scale}} for z-score standardization, 
#'          \code{\link{normalize}} in some packages for similar functionality
#'
#' @export
range01 <- function(x){
  # Check if all values are the same to avoid division by zero
  if(max(x) == min(x)) {
    return(rep(0, length(x)))
  }
  (x-min(x))/(max(x)-min(x))
}


#' Merge Overlapping Clusters Based on Intersection Table
#'
#' This function identifies connected components of overlapping clusters using
#' graph analysis and merges them into consolidated clusters. Clusters that
#' overlap according to the go_table are merged into single clusters with
#' Biomarker "Merged" and new cluster names ordered by size.
#'
#' @param scores A data.table containing cluster data with columns:
#'   Biomarker, symbol, cluster, and other associated data (GENE, ix, degree, 
#'   strength, alpha_scaled, etc.)
#' @param go_table A data.table containing overlap information with columns:
#'   dx.1, clust.1, dx.2, clust.2. Additional columns (like intersection) 
#'   are ignored. Each row represents an overlap between two clusters.
#'
#' @return A data.table with the same structure as the input scores table,
#'   plus an additional 'orig_groups' column. Overlapping clusters are merged
#'   into "Merged" biomarker clusters numbered M1, M2, M3... in descending 
#'   order by size. Non-overlapping clusters remain unchanged.
#'
#' @details 
#' The function performs the following steps:
#' \itemize{
#'   \item Creates an undirected graph from the overlap relationships
#'   \item Identifies connected components using graph analysis
#'   \item Groups clusters that belong to the same connected component
#'   \item Merges all scores data for clusters in each component
#'   \item Assigns new cluster names (M1, M2, ...) based on cluster size
#'   \item Creates orig_groups column with original Biomarker-cluster combinations
#'   \item Preserves all non-overlapping clusters with their original names
#'   \item Returns combined dataset with proper key ordering
#' }
#'
#' The orig_groups column format is "Biomarker1-cluster1;Biomarker2-cluster2"
#' for merged clusters, and "Biomarker-cluster" for non-overlapping clusters.
#'
#' @note Requires the igraph package for connected components analysis.
#'   Empty go_table will result in all clusters remaining unchanged.
#'
#' @examples
#' # Basic usage with overlap table
#' result <- merge_overlapping_clusters(scores, go_table)
#' 
#' # Check merged clusters
#' result[Biomarker == "Merged", .N, by = cluster][order(cluster)]
#' 
#' # View original groups for merged clusters  
#' result[Biomarker == "Merged", .(cluster, orig_groups)] |> unique()
#' 
#' # Count original vs merged clusters
#' result[, .N, by = .(Biomarker == "Merged")]
#'
#' @seealso 
#' \code{\link{merge_by_term_reorder}} for additional term-based merging
#' \code{\link[igraph]{components}} for connected components algorithm
#'
#' @import data.table
#' @import igraph
#' @export
merge_overlapping_clusters <- function(scores, go_table) {
  
  # Input validation
  if (!is.data.table(scores)) {
    stop("scores must be a data.table")
  }
  
  if (!is.data.table(go_table)) {
    stop("go_table must be a data.table")
  }
  
  # Check required columns in scores
  required_scores_cols <- c("Biomarker", "symbol", "cluster")
  missing_scores_cols <- setdiff(required_scores_cols, names(scores))
  if (length(missing_scores_cols) > 0) {
    stop("Missing required columns in scores: ", paste(missing_scores_cols, collapse = ", "))
  }
  
  # Check required columns in go_table
  required_go_cols <- c("dx.1", "clust.1", "dx.2", "clust.2")
  missing_go_cols <- setdiff(required_go_cols, names(go_table))
  if (length(missing_go_cols) > 0) {
    stop("Missing required columns in go_table: ", paste(missing_go_cols, collapse = ", "))
  }
  
  # Handle empty go_table case
  if (nrow(go_table) == 0) {
    stop("Empty go_table provided...")
  }
  
  message("Processing ", nrow(go_table), " overlap relationships...")
  
  # Create graph from overlaps
  edges <- go_table[, .(
    from = paste(dx.1, clust.1, sep = "-"), 
    to = paste(dx.2, clust.2, sep = "-")
  )]
  
  # Remove self-loops if any
  edges <- edges[from != to]
  
  if (nrow(edges) == 0) {
    message("No valid edges found after removing self-loops. Adding orig_groups column.")
    result <- copy(scores)
    result[, orig_groups := paste(Biomarker, cluster, sep = "-")]
    setkey(result, Biomarker, symbol, cluster)
    return(result)
  }
  
  # Create undirected graph and find connected components
  g <- graph_from_data_frame(edges, directed = FALSE)
  components_result <- components(g)
  
  message("Found ", components_result$no, " connected components")
  
  # Create component membership table
  comp_membership <- data.table(
    biomarker_cluster = names(components_result$membership),
    component = components_result$membership
  )
  
  # Split biomarker_cluster back into separate columns
  comp_membership[, c("Biomarker", "cluster") := tstrsplit(biomarker_cluster, "-", fixed = TRUE)]
  comp_membership[, biomarker_cluster := NULL]
  
  # Add component information to scores table
  scores_with_comp <- merge(scores, comp_membership, 
                           by = c("Biomarker", "cluster"), 
                           all.x = TRUE)
  
  # Separate overlapping and non-overlapping clusters
  overlapping <- scores_with_comp[!is.na(component)]
  non_overlapping <- scores_with_comp[is.na(component)]
  
  message("Clusters to merge: ", nrow(overlapping), 
          ", Clusters to keep unchanged: ", nrow(non_overlapping))
  
  # Process overlapping clusters
  if (nrow(overlapping) > 0) {
    # Create orig_groups string from unique biomarker-cluster combinations
    merged_clusters <- overlapping[, `:=`(
      orig_groups = paste(
        unique(paste(Biomarker, cluster, sep = "-")), collapse = ";"
      ),
      orig_terms = paste(unique(term), collapse = ";")
    ),
    by = component]
    
    # Calculate cluster sizes for ranking
    cluster_sizes <- merged_clusters[, .(size = .N), by = component]
    cluster_sizes[, new_cluster := paste0("M", frank(-size, ties.method = "first"))]
    
    # Merge cluster size info and assign final names
    merged_clusters <- merge(merged_clusters, 
                            cluster_sizes[, .(component, new_cluster)], 
                            by = "component")
    
    # Set final cluster attributes
    merged_clusters[, `:=`(
      Biomarker = "Merged",
      cluster = new_cluster,
      term = orig_terms,
      component = NULL,
      new_cluster = NULL
    )]
  } else {
    # Create empty merged_clusters with correct structure
    merged_clusters <- data.table()
  }
  
  # Process non-overlapping clusters (add orig_groups & orig_terms columns)
  if (nrow(non_overlapping) > 0) {
    non_overlapping[, `:=`(
      orig_groups = paste(Biomarker, cluster, sep = "-"),
      orig_terms = term,
      component = NULL
    )]
  }
  
  # Combine results
  if (nrow(merged_clusters) > 0 && nrow(non_overlapping) > 0) {
    # Get common columns
    common_cols <- intersect(names(merged_clusters), names(non_overlapping))
    result <- rbind(
      merged_clusters[, ..common_cols],
      non_overlapping[, ..common_cols]
    )
  } else if (nrow(merged_clusters) > 0) {
    result <- merged_clusters
  } else if (nrow(non_overlapping) > 0) {
    result <- non_overlapping
  } else {
    stop("Unexpected error: both merged and non-overlapping clusters are empty")
  }
  
  # Set key for efficient operations
  setkey(result, Biomarker, symbol, cluster)
  
  message("Merging complete. Final result has ", nrow(result), " rows with ",
          length(unique(result[Biomarker == "Merged"]$cluster)), " merged clusters")
  
  return(result)
}

#' Merge Clusters by Term Pattern and Reorder by Size
#'
#' This function identifies clusters that match a given regex pattern in their 
#' term column, merges them into a single new cluster with Biomarker "Merged", 
#' and reorders all merged clusters by size (largest first).
#'
#' @param result_table A data.table containing cluster data with columns:
#'   Biomarker, symbol, cluster, term, orig_groups, and other associated data
#' @param term_regex A character string containing a regular expression pattern
#'   to match against the term column (case-insensitive matching)
#' @param new_term_name A character string specifying the new term name for 
#'   the merged cluster
#'
#' @return A data.table with the same structure as input, but with matching 
#'   clusters merged into a new "Merged" biomarker cluster. All "Merged" 
#'   clusters are renumbered (M1, M2, M3...) based on size in descending order.
#'   The orig_groups column is updated to contain all original biomarker-cluster
#'   combinations that were merged.
#'
#' @details 
#' The function performs the following steps:
#' \itemize{
#'   \item Identifies clusters where the term column matches the regex pattern
#'   \item Merges all matching clusters into a single cluster
#'   \item Combines orig_groups from all merged clusters
#'   \item Assigns "Merged" as the biomarker and the new_term_name as the term
#'   \item Recalculates cluster numbers for all "Merged" clusters based on size
#'   \item Returns the updated table with proper key ordering
#' }
#'
#' @examples
#' # Merge inflammation-related clusters
#' result <- merge_by_term_reorder(result, 
#'                                term_regex = "inflamm|cytokine|interleukin", 
#'                                new_term_name = "Inflammation")
#'
#' # Merge metabolic clusters
#' result <- merge_by_term_reorder(result,
#'                                term_regex = "metabol|glucose|insulin|lipid",
#'                                new_term_name = "Metabolism")
#'
#' # Merge transport-related clusters
#' result <- merge_by_term_reorder(result,
#'                                term_regex = "transport|channel|carrier",
#'                                new_term_name = "Transport")
#'
#' @seealso 
#' \code{\link{merge_overlapping_clusters}} for merging based on overlap tables
#'
#' @import data.table
#' @export
merge_by_term_reorder <- function(result_table, term_regex, new_term_name) {
  
  # Input validation
  if (!is.data.table(result_table)) {
    stop("result_table must be a data.table")
  }
  
  required_cols <- c("Biomarker", "term", "orig_groups", "orig_terms")
  missing_cols <- setdiff(required_cols, names(result_table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (!is.character(term_regex) || length(term_regex) != 1) {
    stop("term_regex must be a single character string")
  }
  
  if (!is.character(new_term_name) || length(new_term_name) != 1) {
    stop("new_term_name must be a single character string")
  }
  
  # Identify matching clusters
  matching_clusters <- result_table[grepl(term_regex, term)]
  non_matching_clusters <- result_table[!grepl(term_regex, term)]
  
  if (nrow(matching_clusters) == 0) {
    message("No clusters found matching the regex pattern: ", term_regex)
    return(result_table)
  }
  
  message("Found ", nrow(matching_clusters), " rows matching pattern '", 
          term_regex, "', merging into '", new_term_name, "'")
  
  # Create merged cluster
  merged_by_term <- matching_clusters[, {
    # Handle orig_groups combination based on whether clusters are already merged
    if (all(Biomarker == "Merged")) {
      # All are already merged clusters, combine their orig_groups
      combined_orig_groups = paste(unique(orig_groups), collapse = ";")
      combined_orig_terms = paste(unique(orig_terms), collapse = ";")
    } else {
      # Mix of merged and original clusters
      orig_parts <- ifelse(
        Biomarker == "Merged", 
        unique(orig_groups),
        unique(paste(Biomarker, cluster, sep = "-"))
      )
      combined_orig_groups = paste(unique(orig_parts), collapse = ";")
      # For orig_terms, use existing terms or the original term names
      term_parts <- ifelse(
        Biomarker == "Merged",
        unique(orig_terms),
        unique(term)
      )
      combined_orig_terms = paste(unique(term_parts), collapse = ";")
    }
    
    # Return all rows with updated grouping information
    .SD[, .(symbol, GENE, ix, beta, degree, strength, alpha_scaled, 
            combined_orig_groups, combined_orig_terms)]
  }]
  
  # Set merged cluster attributes
  merged_by_term[, `:=`(
    Biomarker = "Merged",
    term = new_term_name,
    orig_groups = combined_orig_groups,
    orig_terms = combined_orig_terms,
    combined_orig_groups = NULL,
    combined_orig_terms = NULL
  )]
  
  # Combine all clusters
  temp_result <- rbind(non_matching_clusters, merged_by_term, fill = TRUE)
  
  # Reorder all Merged clusters by size (largest first)
  merged_clusters <- temp_result[Biomarker == "Merged"]
  other_clusters <- temp_result[Biomarker != "Merged"]
  
  if (nrow(merged_clusters) > 0) {
    # Calculate cluster sizes and assign new numbers
    cluster_sizes <- merged_clusters[, .N, by = term]
    cluster_sizes[, new_cluster := paste0("M", frank(-N, ties.method = "first"))]
    
    # Update cluster names
    merged_clusters <- merge(merged_clusters, 
                            cluster_sizes[, .(term, new_cluster)], 
                            by = "term")
    merged_clusters[, `:=`(cluster = new_cluster, new_cluster = NULL)]
  }
  
  # Create final result
  result <- rbind(other_clusters, merged_clusters, fill = TRUE)
  
  # Set key for efficient operations
  setkey(result, Biomarker, symbol, cluster)
  
  return(result)
}

#' Complete workflow for merging clusters
#'
#' @param scores Original scores data.table
#' @param go_table Overlap table for initial merging
#' @param term_merges Named list of regex patterns and new term names
#'
#' @return Fully processed data.table with all merges applied
#'
#' @examples
#' # Define term-based merges
#' term_merges <- list(
#'   "Immune Response" = "immune|inflamm|cytokine|interleukin",
#'   "Metabolism" = "metabol|glucose|insulin|lipid|energy",
#'   "Transport" = "transport|channel|carrier|pump",
#'   "Neuronal" = "neuron|synap|axon|dendrit|neural"
#' )
#' 
#' result <- process_all_merges(scores, go_table, term_merges)
process_all_merges <- function(scores, go_table, term_merges) {
  
  # Start with overlap-based merging
  result <- merge_overlapping_clusters(scores, go_table)
  
  # Apply each term-based merge
  for (new_term in names(term_merges)) {
    regex_pattern <- term_merges[[new_term]]
    result <- merge_by_term_reorder(result, regex_pattern, new_term)
  }
  
  return(result)
}


# Main -----------------------------------------------------------------------
# Stage data for analysis --------------------------
# set constant values
tbls <- c("nodes", "edges", "members")
path <- file.path("baseline_PCs_relatedness", input, "STRING")

# read in for each biomarker and join
lapply(tbls, function(i){
  fh <- file.path(path, paste0("Cluster_", i, ".csv"))
  x <- lapply(fh, fread)
  setattr(x, 'names', input)
  x <- rbindlist(x, idcol = 'Biomarker', fill = TRUE, use.names = TRUE)
  assign(i, x, envir = .GlobalEnv)
})

# filter and simplify tables
edges <- edges[selected == TRUE]
nodes <- nodes[selected == TRUE]
edges[, c("shared name", "shared interaction", "selected") := NULL]
nodes <- nodes[, .(Biomarker, beta, `display name`, GENE, ix, name)]


# Join and check -----------------------------------
scores <- merge(
  x = members, 
  y = nodes, 
  by.x = c('Biomarker', 'symbol'),
  by.y = c('Biomarker', 'display name'),
  all.x = TRUE
)
setnames(scores, "name", "string.name")
scores[, abs.beta := abs(beta)]
cols = c("symbol", names(scores)[names(scores) != "symbol"])
scores <- scores[, ..cols]

# Testing --------------------------------
# # check duplicates missing N
# scores[ 
  # duplicated(scores[,.(Biomarker, symbol)]) 
  # | duplicated(scores[,.(Biomarker, symbol)], fromLast=TRUE)
# ]
# # check missing N
# scores[, .N, by = .(Biomarker, cluster)
       # ][, .SD[mixedorder(cluster)], by = Biomarker
         # ][, print(.SD, nrows = Inf)]

# # check missing edges N
# edges[, .N, by = .(Biomarker)
       # ][, print(.SD, nrows = Inf)]

# This assumes non-overlapping cluster membership, 
# otherwise edges need a more complex cluster matching
# -----------------

# give edge table node1 & node2 as symbols
edges[, c("name1", "name2") := tstrsplit(name, " (pp) ", fixed = TRUE)]
mp <- unique(scores[, .(symbol, string.name)]) # mapping table
edges[mp, node1 := symbol, on = c(name1 = "string.name")]
edges[mp, node2 := symbol, on = c(name2 = "string.name")]

# simplify edge table
setnames(edges, "stringdb::score", "stringdb.score")
edges <- edges[, .(node1, node2, Biomarker, stringdb.score)]
fwrite(edges, file = file.path(output, paste0("STRING_raw_edge_table.csv")))


# Node strength/degree calculation -------------------------------------
results <- lapply(input, function(i) {
  node_metrics(edges[Biomarker == i], scores[Biomarker == i])
})
setattr(results, 'names', input)
results <- rbindlist(results, idcol = "Biomarker")

# merge into scoring table
scores <- merge(
  x = scores, 
  y = results, 
  by.x = c("Biomarker", "symbol", "cluster"), 
  by.y = c("Biomarker", "node", "cluster"),
  all.x = TRUE
)
scores[, alpha_scaled := range01(alpha_0.2), by = Biomarker]
mat <- dcast(scores, symbol ~ Biomarker, value.var = 'alpha_scaled', fill = 0)
spec_joy(mat, 
  gen_col = "symbol", 
  outname = file.path(output,
    paste0("STRING_alpha_scaled_raw")
  )
)
fwrite(scores, file = file.path(output, paste0("STRING_raw_scaled_scores.csv")))

# Collapsing shared modules by overlap enrichment ------------------
# read in overlaps from manuscript
go_table <- fread(ref)
go_table <- go_table[padj < 0.05 & intersection > 2]

# simplify scores and add annotations
fh <- file.path(path, paste0("Annotated_clusters.csv"))
x <- lapply(fh, fread)
setattr(x, 'names', input)
x <- rbindlist(x, idcol = 'Biomarker', fill = TRUE, use.names = TRUE)
scores[, c("abs.beta", "alpha_0.1", "alpha_0.2", "alpha_0.4", "string.name") := NULL]
scores <- scores[x[, .(Biomarker, cluster, term)], on = c("Biomarker", "cluster")]

# Basic overlap merging only
result <- merge_overlapping_clusters(scores, go_table)

# Check for overlapping terms for additional term merging
term_patterns <- list(
  "Actin dynamics" = "Actin",
  "BMP signaling" = "BMP signaling",
  "Carbon metabolism" = "Carbon",
  "Chaperone|Protein folding" = "Chaperone|Chaperonin TCP",
  "Cadherin|Cell adhesion" = "Cell adhes",
  "Focal adhesion|Integrins" = "Focal|Integr",
  "Complement cascade" = "Complement",
  "Cell cycle" = "Cell cycle",
  "Glutathione metabolism|Oxidative stress" = "Glutathione",
  "Glycosylation" = "Glycosylation",
  "ER-Golgi transport" = "Golgi",
  "Hippo signaling" = "Hippo signaling",
  "Adaptive immunity|Immunoglobulins|Cytokines" = "Immunoglobulin",
  "Matrix metalloprotease" = "Matrix metalloprotease",
  "Neuroactive ligand-receptor interaction" = "Neuroactive",
  "NFkB|APC/C degradation" = "NFkB\\|APC/C degradation",
  "Neddylation" = "Neddylation",
  "Nucleolus|Translation" = "Nucleolus\\|Translation",
  "RHO GTPase|Neurotrophin signaling" = "75NTR|Neurotrop",
  "RNA-binding" = "RNA-binding",
  "Serine protease inhibitor" = "Serine protease inhibitor",
  "Synaptic vesicle|SNARE complex" = "Synaptic vesicle\\|SNARE complex",
  "Synapse assembly" = "Synapse assembly",
  "Translation" = "Translation"
)

result <- process_all_merges(scores, go_table, term_patterns)

# Drop Ambiguous & non-Brain terms ------------
result <- result[term != "Ambiguous"]

# Merged set naming
result[, term := str_replace_all(term, fixed(c(
  "TGF-beta signaling;PI3K-Akt|MAPK signaling;JAK-STAT signaling;Ambiguous;PI3K-Akt|JAK-STAT|MAPK signaling;Cytokine|JAK-STAT signaling;Hormone response"="PI3K-Akt|JAK-STAT|MAPK signaling",
  "Post-chaperonin tubulin folding;Nucleolus"="Post-chaperonin tubulin folding"
)))]

# dropping non-brain associated terms
drops <- c(
  "Muscle contraction",
  "Blood microparticle",
  "Leukocyte transendothelial migration",
  "Digestion of dietary lipid",
  "Cell wall disruption",
  "Antrax receptor",
  "Antigen processing and presentation;Graft-versus-host disease"
)
result <- result[!term %in% drops]

# Calculating mean values for (node degree, strength, alpha) -------
cols <- c("Biomarker", "symbol", "cluster", "GENE", "ix", "term", "orig_groups", "orig_terms")
mean_scores <- result[, .(
  mean_beta = mean(beta),
  mean_degree = mean(degree), 
  mean_strength = mean(strength),
  mean_alpha_scaled = mean(alpha_scaled) 
), by = cols]
# inside clusters, take the highest mean_node_degree
setorder(mean_scores, -mean_alpha_scaled)
mean_scores[, .SD[1], by = .(Biomarker, symbol, cluster)]

# STRING db Symbol error, labeled CYRIB as FAM49B
mean_scores[symbol == "FAM49B", symbol := "CYRIB"]
# This one is on the SomaLogic folks (should be COXFA4L3)
mean_scores[symbol == "C15orf48", symbol := "NMES1"]

fwrite(mean_scores, 
  file = file.path(output, paste0("STRING_merged_scaled_scores.csv"))
)
# term mappings for WSS
eig <- unique(mean_scores[, .(Biomarker, cluster, term)])
eig[, name := paste(Biomarker, cluster, sep = ".")]
fwrite(eig[, .(name, term, group = Biomarker)], 
  file = file.path(output, paste0("STRING.INT.WSS.terms.csv"))
)

# Save and generate specificity matrices/joy plots  ----------------
mat <- dcast(mean_scores, 
  symbol ~ Biomarker, 
  value.var = 'mean_alpha_scaled', 
  fill = 0,
  fun.aggregate = first
)

spec_joy(mat, 
  gen_col = "symbol", 
  outname = file.path(output,
    paste0("STRING_alpha_scaled_merged")
  )
)


# module split terms -------------------
# Limit to terms where N >= 10 for GSEA
top <- mean_scores[,.N, by = term][N >= 10, term]
mean_scores <- mean_scores[term %in% top]

mat <- dcast(mean_scores, 
  symbol ~ term, 
  value.var = 'mean_alpha_scaled', 
  fill = 0,
  fun.aggregate = first
)
spec_joy(mat, 
  gen_col = "symbol", 
  outname = file.path(output,
    paste0("STRING_alpha_scaled_merged_terms_n10")
  )
)

# included SOMA proteins that not in any cluster as background
names(mat)[1] <- "GENE"
prot <- fread(prots)
# split string and explode piped terms into separate rows
mat2 <- merge(prot[, .(GENE)], mat, all = TRUE)
mat2[is.na(mat2)] <- 0

fwrite(mat2, file = file.path(output, "STRING_alpha_scaled.spec_matr.txt"))


# Expore WSS for clean distributions...
