#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 08a-String_module_testing.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Test STRING network modules and prepar for Scoring
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
library(data.table)
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(igraph)
library(RColorBrewer)
library(Cairo)

# read in data (hard coded for now)
test <- fread("module_scores/STRING_merged_scaled_scores.csv")

# Corrected function with proper factor ordering
create_ranking_flow_data <- function(dt, alpha_cols = c("alpha_0.1", "alpha_0.2", "alpha_0.4")) {
  
  # Calculate ranks for each alpha (higher centrality = lower rank number)
  dt_ranks <- copy(dt)
  for (col in alpha_cols) {
    rank_col <- paste0(col, "_rank")
    dt_ranks[, (rank_col) := rank(-get(col), ties.method = "min")]
  }
  
  # Create ranking bins with proper factor levels
  n_nodes <- nrow(dt_ranks)
  
  # Define the correct order for ranking tiers
  tier_levels <- c("Top 10%", "Top 25%", "Top 50%", "Top 75%", "Bottom 25%")
  
  # Create bins for each alpha
  for (alpha_col in c("alpha_0.1", "alpha_0.2", "alpha_0.4")) {
    rank_col <- paste0(alpha_col, "_rank")
    bin_col <- paste0(alpha_col, "_bin")
    
    dt_ranks[, (bin_col) := cut(get(rank_col), 
      breaks = c(0, n_nodes*0.1, n_nodes*0.25, n_nodes*0.5, n_nodes*0.75, n_nodes),
      labels = tier_levels,
      include.lowest = TRUE)
    ]
    
    # Explicitly set factor levels in correct order
    dt_ranks[, (bin_col) := factor(get(bin_col), levels = tier_levels)]
  }
  
  return(dt_ranks)
}

# Create the flow data
setorder(test, -alpha_0.1)
test_with_ranks <- create_ranking_flow_data(test)

# Prepare data for alluvial plot
prepare_alluvial_data <- function(dt_ranks) {
  
  # Reshape for alluvial plot
  flow_data <- dt_ranks[, .(
    id,
    alpha_0.1 = alpha_0.1_bin,
    alpha_0.2 = alpha_0.2_bin, 
    alpha_0.4 = alpha_0.4_bin
  )]
  
  # Convert to long format
  flow_long <- melt(flow_data, 
    id.vars = "id", 
    variable.name = "alpha_setting", 
    value.name = "rank_bin"
  )
  
  return(flow_long)
}
alluvial_data <- prepare_alluvial_data(test_with_ranks)

# Create the alluvial plot
plot_ranking_flow <- function(alluvial_data) {
  # Create origin rank_bin (rank_bin at alpha_0.1 for each id)
  origin_mapping <- alluvial_data[alpha_setting == "alpha_0.1", .(id, origin_rank = rank_bin)]
  
  # Merge back to get origin for all rows
  alluvial_data <- merge(alluvial_data, origin_mapping, by = "id")  

  ggplot(alluvial_data, 
      aes(x = alpha_setting, stratum = rank_bin, alluvium = id)) +
    geom_flow(aes(fill = origin_rank),
      stat = "alluvium", 
      lode.guidance = "frontback",
      # color = "darkgray", 
      alpha = 0.7
    ) +
    geom_stratum(aes(fill = rank_bin), alpha = 0.9, color = "white", size = 0.5) +
    geom_label(
        aes(label = rank_bin, fill = rank_bin), 
      stat = "stratum", 
      size = 3, 
      color = "white", 
      fontface = "bold"
    ) +
    scale_fill_viridis_d(
      name = "Ranking Tier", 
      option = "plasma", 
      direction = -1
    ) +
    labs(
      title = "How Node Rankings Change Across Alpha Values",
      subtitle = "Flow shows movement between ranking tiers as alpha increases",
      x = "Alpha Setting",
      y = "Number of Nodes"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom"
    ) +
    scale_x_discrete(labels = c("α = 0.1", "α = 0.2", "α = 0.4"))
}


# Create the plot
ranking_flow_plot <- plot_ranking_flow(alluvial_data)
print(ranking_flow_plot)






# Analyze specific movement patterns -------------------------------
analyze_ranking_movements <- function(dt_ranks) {
  
  # Calculate rank changes
  dt_ranks[, `:=`(
    change_0.1_to_0.2 = alpha_0.2_rank - alpha_0.1_rank,
    change_0.2_to_0.4 = alpha_0.4_rank - alpha_0.2_rank,
    total_change = alpha_0.4_rank - alpha_0.1_rank
  )]
  
  # Categorize movements
  dt_ranks[, movement_type := fcase(
    total_change < -50, "Major Improver",
    total_change < -10, "Moderate Improver", 
    abs(total_change) <= 10, "Stable",
    total_change > 10 & total_change <= 50, "Moderate Decliner",
    total_change > 50, "Major Decliner"
  )]
  
  # Summary statistics
  movement_summary <- dt_ranks[, .(
    count = .N,
    avg_change = mean(total_change),
    median_change = median(total_change)
  ), by = movement_type]
  
  return(list(
    detailed = dt_ranks,
    summary = movement_summary
  ))
}
movement_analysis <- analyze_ranking_movements(test_with_ranks)
print(movement_analysis$summary)

# Plot movement distribution
plot_movement_distribution <- function(movement_data) {
  
  ggplot(movement_data$detailed, aes(x = total_change, fill = movement_type)) +
    geom_histogram(bins = 50, alpha = 0.7, color = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
    labs(
      title = "Distribution of Ranking Changes (α=0.1 to α=0.4)",
      subtitle = "Negative values = improved ranking, Positive values = declined ranking",
      x = "Change in Rank Position", 
      y = "Number of Nodes",
      fill = "Movement Type"
    ) +
    scale_fill_viridis_d(option = "plasma") +
    theme_minimal() +
    facet_wrap(~movement_type, scales = "free_y")
}
movement_plot <- plot_movement_distribution(movement_analysis)
print(movement_plot)



# Top movers analysis -------------------------------------------
identify_top_movers <- function(dt_ranks, n_top = 10) {
  
  # Top improvers (biggest negative change = moved up most)
  top_improvers <- dt_ranks[order(total_change)][1:n_top, 
    .(id, alpha_0.2_rank, alpha_0.4_rank, total_change, improvement = -total_change)
  ]
  
  # Top decliners (biggest positive change = moved down most)  
  top_decliners <- dt_ranks[order(-total_change)][1:n_top,
    .(id, alpha_0.2_rank, alpha_0.4_rank, total_change, decline = total_change)
  ]
  
  return(list(improvers = top_improvers, decliners = top_decliners))
}
top_movers <- identify_top_movers(test_with_ranks, n_top = 15)
top_movers

# Create focused alluvial plot for top movers
plot_top_movers_flow <- function(dt_ranks, top_movers, n_show = 20) {
  
  # Get IDs of biggest movers
  mover_ids <- c(
    top_movers$improvers[1:(n_show/2)]$id,
    top_movers$decliners[1:(n_show/2)]$id
  )
  
  # Filter data to only show these nodes
  focused_data <- dt_ranks[id %in% mover_ids]
  
  # Prepare for alluvial
  flow_focused <- focused_data[, .(
    id,
    alpha_0.1 = paste("Rank", alpha_0.1_rank),
    alpha_0.2 = paste("Rank", alpha_0.2_rank),
    alpha_0.4 = paste("Rank", alpha_0.4_rank)
  )]
  
  flow_long_focused <- melt(flow_focused, 
    id.vars = "id", 
    variable.name = "alpha_setting", 
    value.name = "rank_position"
  )
  
  # Add movement type for coloring
  flow_long_focused[, movement_type := ifelse(
    id %in% top_movers$improvers$id, "Top Improvers", "Top Decliners"
  )]
  
  ggplot(flow_long_focused, 
      aes(x = alpha_setting, stratum = rank_position, alluvium = id, fill = movement_type)
    ) +
    geom_flow(stat = "alluvium", alpha = 0.8) +
    geom_stratum(alpha = 0.5, color = "white") +
    scale_fill_manual(values = c("Top Improvers" = "#2E8B57", 
                                "Top Decliners" = "#DC143C")) +
    labs(
      title = "Ranking Flows for Biggest Movers",
      subtitle = paste("Top", n_show, "nodes with largest ranking changes"),
      x = "Alpha Setting",
      fill = "Movement Type"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}
top_movers_plot <- plot_top_movers_flow(test_with_ranks, top_movers, n_show = 20)
print(top_movers_plot)


# Example merged plot --------------------------------------------------------
# Not used for module scores, but might be a good visual representation. 
# Would need Mean PPI, derived from...
edges <- fread("module_scores/STRING_raw_edge_table.csv")

# lookup key for gene membership in original clusters
scores <- fread("module_scores/STRING_raw_scaled_scores.csv")
scores <- scores[, bc := paste(Biomarker, cluster, sep = "-")][, .(symbol, bc)]

# build networks for mean_PPI as weights ----------------
# Source - https://stackoverflow.com/q/78808739
# Posted by dragon951, modified by community. See post 'Timeline' for change history
# Retrieved 2026-05-12, License - CC BY-SA 4.0

# explode columns with multiple entries
mycol <- c("orig_groups", "orig_terms")
test_exp <- test[, lapply(.SD, \(x) unlist(tstrsplit(x, ";"))), 
  by = c(names(test)[!names(test) %in% mycol]), 
  .SDcols = mycol[1]
]

# Create lookup tables (orig_groups must contain symbol in scores)
lookup1 <- unique(
  test_exp[scores, 
    .(symbol, orig_groups, merge_bio = Biomarker, merge_clus = cluster),
    on = .(symbol = symbol, orig_groups = bc)
  ]
)
setnames(lookup1, c("node1", "bc1", "merge_bio", "merge_clus"))

lookup2 <- unique(
  test_exp[scores, .(symbol, orig_groups), on = .(symbol = symbol, orig_groups = bc)]
)
setnames(lookup2, c("node2", "bc2"))

# Join edges with both lookups
edge_exp <- edges[lookup1, on = "node1", allow.cartesian = TRUE, nomatch = 0]
edge_exp <- edge_exp[lookup2, on = "node2", allow.cartesian = TRUE, nomatch = 0]

# keep only edges where both nodes are in the same orig_group-cluster and map to the original edge table Biomarker
edge_filt <- na.omit(
  edge_exp[bc1 == bc2 & mapply(grepl, Biomarker, bc1, fixed = TRUE)]
)

# Calculate mean PPI for each edge within each biomarker-cluster group
edge_merged <- edge_filt[, .(
  mean_PPI = mean(stringdb.score),
  n_sources = .N,
  source_biomarkers = paste(unique(bc1), collapse = ";")
), by = .(node1, node2, merge_bio, merge_clus)]


# build node graphs -------------------------
# for genes in multiple modules, keep by top mean_alpha_specificity
result <- test[order(-mean_alpha_scaled)][, .SD[1], by = symbol]
result[, module := paste(Biomarker, cluster, sep = ".")]

# filter singles/doubles
keep <- result[, .N, by = module][N > 2, module]
result <- result[module %in% keep, 
  .(symbol, Biomarker, cluster, ix, GENE, mean_beta, mean_alpha_scaled, module, term)
]

# Filter edges to only include those where BOTH nodes are in the correct module
valid_set <- result[, paste(symbol, Biomarker, cluster, sep = ":::")]
edge_final <- edge_merged[
  paste(node1, merge_bio, merge_clus, sep = ":::") %in% valid_set &
  paste(node2, merge_bio, merge_clus, sep = ":::") %in% valid_set
]

# first edges
g_all <- graph_from_data_frame(
  d = edge_final[, .(
    node1, 
    node2, 
    weight = mean_PPI, 
    Biomarker = merge_bio, 
    cluster = merge_clus,
    source_biomarkers
  )],
  directed = FALSE
)

# Add node attributes directly from result
V(g_all)$Biomarker <- result$Biomarker[match(V(g_all)$name, result$symbol)]
V(g_all)$cluster <- result$cluster[match(V(g_all)$name, result$symbol)]
V(g_all)$ix <- result$ix[match(V(g_all)$name, result$symbol)]
V(g_all)$GENE <- result$GENE[match(V(g_all)$name, result$symbol)]
V(g_all)$mean_beta <- result$mean_beta[match(V(g_all)$name, result$symbol)]
V(g_all)$mean_alpha_scaled <- result$mean_alpha_scaled[match(V(g_all)$name, result$symbol)]
V(g_all)$module <- result$module[match(V(g_all)$name, result$symbol)]
V(g_all)$term <- result$term[match(V(g_all)$name, result$symbol)]



# Plot the graphs ---------------------------
# Create color palettes
# For nodes: gradient based on abs(mean_beta)
node_colors <- colorRampPalette(c("white", "#8d191d"))(100)

# For edges: gradient based on mean_PPI
edge_colors <- colorRampPalette(c("lightgray", "#7c7559"))(100)

# Function to map values to colors
map_to_color <- function(values, colors) {
  # Normalize values to 0-1 range
  normalized <- (values - min(values, na.rm = TRUE)) / 
                (max(values, na.rm = TRUE) - min(values, na.rm = TRUE))
  # Map to color indices
  color_indices <- cut(normalized, breaks = length(colors), labels = FALSE, include.lowest = TRUE)
  colors[color_indices]
}

# Assign colors to nodes and edges
V(g_all)$color <- map_to_color(abs(V(g_all)$mean_beta), node_colors)
E(g_all)$color <- map_to_color(E(g_all)$weight, edge_colors)

# Set edge width based on weight
E(g_all)$width <- scales::rescale(E(g_all)$weight, to = c(1, 5))

# # Set node size based on mean_alpha_scaled
# V(g_all)$size <- scales::rescale(V(g_all)$mean_alpha_scaled, to = c(3, 10))


# Function to create grid layout with clusters
layout_clusters_grid <- function(g, cluster_attr = "module", 
                                 n_cols = 10, 
                                 cluster_spacing = 3) {
  
  # Get cluster assignments
  clusters <- vertex_attr(g, cluster_attr)
  
  # Count nodes per cluster and sort by size (largest first)
  cluster_sizes <- table(clusters)
  cluster_order <- names(sort(cluster_sizes, decreasing = TRUE))
  
  # Initialize layout matrix
  layout_matrix <- matrix(0, nrow = vcount(g), ncol = 2)
  
  # Calculate grid positions for each cluster
  for (i in seq_along(cluster_order)) {
    cluster_name <- cluster_order[i]
    
    # Get nodes in this cluster
    cluster_nodes <- which(clusters == cluster_name)
    n_nodes <- length(cluster_nodes)
    
    # Calculate grid position (row, col) for this cluster
    grid_row <- floor((i - 1) / n_cols)
    grid_col <- (i - 1) %% n_cols
    
    # Create subgraph for this cluster
    g_sub <- induced_subgraph(g, cluster_nodes)
    
    # Layout this cluster using FR or another algorithm
    if (n_nodes > 1) {
      sub_layout <- layout_with_fr(g_sub)
      # Normalize to unit square
      sub_layout <- apply(sub_layout, 2, function(x) {
        (x - min(x)) / (max(x) - min(x) + 1e-10)
      })
    } else {
      sub_layout <- matrix(c(0.5, 0.5), nrow = 1)
    }
    
    # Offset by grid position
    offset_x <- grid_col * cluster_spacing
    offset_y <- -grid_row * cluster_spacing  # Negative for top-to-bottom
    
    layout_matrix[cluster_nodes, 1] <- sub_layout[, 1] + offset_x
    layout_matrix[cluster_nodes, 2] <- sub_layout[, 2] + offset_y
  }
  
  return(layout_matrix)
}

# Create the layout
custom_layout <- layout_clusters_grid(g_all, 
  cluster_attr = "module", 
  cluster_spacing = 2.5
)



# Save full network as PDF
CairoPDF("module_scores/STRING_merged_network_all.pdf", width = 16, height = 20)
par(mar = c(0, 0, 2, 0))

plot(g_all,
   vertex.color = V(g_all)$color,
   vertex.size = 2,
   # vertex.label = NA,
   vertex.label.cex = 0.5,
   vertex.label.color = "black",
   vertex.frame.color = "#8d191d",
   edge.color = E(g_all)$color,
   edge.width = E(g_all)$width,
   layout = custom_layout,
   main = "Merged Protein Networks\nNode color: |mean_beta|, Edge color: mean_PPI"
)

# Add legend for nodes
legend("topright", 
  legend = c(
    paste("High |β|:", round(max(abs(V(g_all)$mean_beta)), 2)),
    paste("Low |β|:", round(min(abs(V(g_all)$mean_beta)), 2))
  ),
  col = c(node_colors[100], node_colors[1]),
  pch = 19, cex = 0.8, pt.cex = 2,
  title = "Node Color")

# Add legend for edges
legend("bottomright",
  legend = c(
    paste("High PPI:", round(max(E(g_all)$weight), 2)),
    paste("Low PPI:", round(min(E(g_all)$weight), 2))
  ),
  col = c(edge_colors[100], edge_colors[1]),
  lwd = 3, cex = 0.8,
  title = "Edge Color"
)

dev.off()


# Extract subgraph for a specific module
plot_single_module <- function(g, biomarker, cluster, 
                               node_col, edge_col,
                               output_file = NULL) {
  
  # Get nodes for this specific biomarker-cluster
  module_nodes <- which(V(g)$Biomarker == biomarker & V(g)$cluster == cluster)
  
  if (length(module_nodes) == 0) {
    cat("No nodes found for", biomarker, cluster, "\n")
    return(NULL)
  }
  
  # Create subgraph
  g_sub <- induced_subgraph(g, vids = module_nodes)
  
  # Skip if too few nodes or edges
  if (vcount(g_sub) < 3) {
    cat("Too few nodes (", vcount(g_sub), ") for", biomarker, cluster, "\n")
    return(NULL)
  }
  
  # Recolor for this subgraph (rescale to this subgraph's range)
  V(g_sub)$color <- map_to_color(abs(V(g_sub)$mean_beta), node_col)
  E(g_sub)$color <- map_to_color(E(g_sub)$weight, edge_col)
  E(g_sub)$width <- scales::rescale(E(g_sub)$weight, to = c(1, 5))
  
  # Build title
  plot_title <- paste0(biomarker, " - ", cluster, "\n",
                      unique(V(g_sub)$term)[1], "\n",
                      vcount(g_sub), " nodes, ", ecount(g_sub), " edges")
  
  # Create plot
  CairoPDF(output_file, width = 7, height = 7)
  
  par(mar = c(1, 1, 3, 1))
  plot(g_sub,
       vertex.color = V(g_sub)$color,
       vertex.size = 5,
       vertex.label = V(g_sub)$symbol,
       vertex.frame.color = "#8d191d",
       edge.color = E(g_sub)$color,
       edge.width = E(g_sub)$width,
       layout = layout_nicely(g_sub),
       main = plot_title
  )
  
  # Add legends
  legend("topright", 
         legend = c(sprintf("High |β|: %.2f", max(abs(V(g_sub)$mean_beta))),
                    sprintf("Low |β|: %.2f", min(abs(V(g_sub)$mean_beta)))),
         col = c(node_col[100], node_col[1]),
         pch = 19, cex = 0.8, pt.cex = 2,
         title = "Node Color", bty = "n")
  
  legend("bottomright",
         legend = c(sprintf("High PPI: %.2f", max(E(g_sub)$weight)),
                    sprintf("Low PPI: %.2f", min(E(g_sub)$weight))),
         col = c(edge_col[100], edge_col[1]),
         lwd = 3, cex = 0.8,
         title = "Edge Color", bty = "n")
  
  dev.off()
  cat("Saved:", output_file, "\n")
}

# Plot Merged M2
plot_single_module(g_all, 
  biomarker = "Merged", 
  cluster = "M1", 
  node_col = node_colors, 
  edge_col = edge_colors,
  output_file = "module_scores/STRING_merged_network_Merged.M1.pdf"
)

plot_single_module(g_all, 
  biomarker = "Merged", 
  cluster = "M2", 
  node_col = node_colors, 
  edge_col = edge_colors,
  output_file = "module_scores/STRING_merged_network_Merged.M2.pdf"
)

plot_single_module(g_all, 
  biomarker = "Merged", 
  cluster = "M8", 
  node_col = node_colors, 
  edge_col = edge_colors,
  output_file = "module_scores/STRING_merged_network_Merged.M8.pdf"
)
plot_single_module(g_all, 
  biomarker = "Merged", 
  cluster = "M26", 
  node_col = node_colors, 
  edge_col = edge_colors,
  output_file = "module_scores/STRING_merged_network_Merged.M26.pdf"
)

# # Simple one-liner to plot Merged.M2
# CairoPDF("cell_type_prot/STRING_merged_network_Merged.M2.pdf", width = 7, height = 7)
# plot(induced_subgraph(g_all, V(g_all)[Biomarker == "Merged" & cluster == "M2"]),
  # vertex.size = 5, 
  # vertex.label = NA, 
  # vertex.color = V(g_all)$color,
  # edge.width = E(g_all)$width,
  # main = "Merged.M2"
# )
# dev.off()

