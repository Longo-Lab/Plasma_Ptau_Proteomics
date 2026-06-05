#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 01-Top_enrichments.R
#'
#' Version: 0.0.3
#'
#' Purpose of script: Compare enrichment hits for biomarker associations in Fig 2
#'
#' Author: Robert B Butler III, Michael Peter Brown
#'
#' Date Created: 2024-09-24
#'
#' Copyright (c) 2024
#' Email: rrbutler@stanford.edu, mpb@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'
#'
#' ---------------------------
#' load up the packages we will need:  (uncomment as required)
library(optparse)
library(gtools)
# library(ggplot2)
library(RColorBrewer)
library(stringr)
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
library(data.table)
library(patchwork)
library(Cairo)
library(GeneOverlap)


# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "List of input datasets comma separated"
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
  )
)

args_list <- parse_args(OptionParser(option_list = option_list))

# Nameset default case
if(is.null(args_list$nameset)){
  args_list$nameset <- file_path_sans_ext(basename(args_list$input))
}

# read in fields
fields <- c(
  "input",
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
#' Plots a circos barplot with links and labels for a given datasets
#'
#' @param dt data.table
#' @param dx_col string, column name for the compared datasets
#' @param x_col int, rank order column name of terms within each dataset (for position in barplot)
#' @param term_col string, column name for item labels 
#' @param sig_col string, column name with significance values
#' @param col_col string, column name for values to color by
#' @param wrap integer, width of characters to wrap text labels
#' @param mtype string, match type for connections. "name", or "enrichment"
#' @param mtable data.table, overlap enrichment table, requires match type "enrichment"
#'
#' @return a base graphics object 
plot_circos_enr <- function(dt, 
                            dx_col = "dx", 
                            x_col = "x", 
                            term_col = "term_name",
                            sig_col = "logp",
                            col_col = "query",
                            wrap = 24,
                            mtype = "name",
                            mtable = FALSE
                            ) {
  # color building
  if (col_col == "up-down") {
    col_col <- 'query'
    bar_cols <- brewer.pal(n = 5, name = "RdBu")[c(1, 5)]
    names(bar_cols) <- c("up", "down")
  } else {
    bar_cols <- brewer.pal(length(unique(dt[[col_col]])), "Paired")
    names(bar_cols) <- unique(dt[[col_col]])
  }
  
  circos.clear()
  circos.par(
    "track.height" = 0.05, 
    gap.degree = 5,
    start.degree = 90,
    circle.margin = c(0.2, 0.2, 0.2, 0.2)
  )
  circos.initialize(
    sectors = levels(dt[[dx_col]]), 
    xlim = cbind(0, dt[, max(get(x_col)), by = dx_col]$V1)
  )

  # term name labels
  circos.labels(
    sectors = dt[[dx_col]], 
    x = dt[[x_col]] - 0.5, 
    labels = str_wrap(dt[[term_col]], width = wrap), 
    side = "outside"
  )

  # bar plots of significance
  circos.track(
    sectors = dt[[dx_col]], 
    ylim = c(0, ceiling(max(dt[[sig_col]]))),
    panel.fun = function(x, y) {
      s.idx <- dt[[dx_col]] == CELL_META$sector.index
      pos <- dt[[x_col]][s.idx] - 0.5
      circos.barplot(
        value = dt[[sig_col]][s.idx], 
        pos = pos, 
        bar_width = 0.6, 
        col = bar_cols[match(dt[[col_col]][s.idx], names(bar_cols))]
      )
  })
  circos.yaxis(sector.index = "GFAP")


  # biomarker ring
  cols <- brewer.pal(n = length(levels(dt[[dx_col]])), name = "Dark2")
  names(cols) <- levels(dt[[dx_col]])
  circos.track(
    ylim = c(0, 1), 
    panel.fun = function(x, y) {
      chr = CELL_META$sector.index
      xlim = CELL_META$xlim
      ylim = CELL_META$ylim
      circos.rect(
        xleft = xlim[1], 
        ybottom = 0, 
        xright = xlim[2], 
        ytop = 1, 
        col = cols[match(chr, names(cols))]
      )
      circos.text(
        x = mean(xlim), 
        y = mean(ylim), 
        labels = chr, 
        col = "white", 
        font = 2, 
        cex = 0.7, 
        facing = "bending.inside", 
        niceFacing = TRUE
      )
    }, 
    track.height = 0.08, 
    bg.border = NA
  )

  # adding links based on the type 
  if (mtype == "enrichment" & is.data.table(mtable)) {
    dlist <- get_links(mtable, dx_col, x_col, term_col, mtype = mtype)
  } else if (mtype == "name") {
    dlist <- get_links(dt, dx_col, x_col, term_col, mtype = mtype)
  } else {
    stop(
      sprintf("incorrect options for mtype: %s and mtable: %s", mtype, mtable),
      call. = FALSE
    )
  }
  lapply(seq(nrow(dlist)), function(i) {
    circos.link(
      sector.index1 = dlist[[i, 1]], 
      point1 = dlist[[i, 2]] - 0.5, 
      sector.index2 = dlist[[i, 3]], 
      point2 = dlist[[i, 4]] - 0.5, 
      col = dlist[[i, 5]], 
      lwd = 5
    )
  })
  
  # adding a legend
  lgd <- Legend(
    at = names(bar_cols), 
    title = col_col, 
    legend_gp = gpar(fill = bar_cols)
  )
  lgd_list_vertical = packLegend(lgd)
  draw(
    lgd_list_vertical, 
    x = unit(4, "mm"), 
    y = unit(4, "mm"), 
    just = c("left", "bottom")
  )
}

#' Generates paired links table for use in circlize
#'
#' @param dt data.table
#' @param dx_col string, column name for the compared datasets
#' @param x_col int, rank order column name of terms within each dataset (for position in barplot)
#' @param term_col string, column name for item labels 
#' @param mtype string, match type for connections. "name", or "enrichment"
#'
#' @return a data.table with the cirlize link format.
get_links <- function(dt, dx_col, x_col, term_col, mtype = "name") {
  if (mtype == "name") {
    # identify duplicates and colors
    dups <- dt[duplicated(get(term_col)) 
              | duplicated(get(term_col), fromLast = TRUE)][order(get(term_col))]
    cols <- colorRampPalette(brewer.pal(9, "PuBu")[-1])(length(unique(dups[[term_col]])))
    names(cols) <- unique(dups[[term_col]])
    dups[, tcols := plyr::mapvalues(
      x = get(term_col), 
      from = names(cols), 
      to = cols
    )]
    # grab pairwise matches
    dlist <- unlist(
      lapply(unique(dups[[term_col]]), function(j) {
        lapply(seq(nrow(dups[get(term_col) == j]) - 1), function(i) {
          cbind(
            dups[get(term_col) == j, c(dx_col, x_col), with = FALSE][1], 
            dups[get(term_col) == j, c(dx_col, x_col, "tcols"), with = FALSE][i+1]
          )
        })
      }), 
      recursive = FALSE
    )
    return(rbindlist(dlist))
  } else if (mtype == "enrichment") {
    # if you are passing a link table, it just goes through with thresh filter
    return(dt)
  }
}


# Main ------------------------------------------
# Gprofiler pathways -------------------------------------------
# DE file locations
fh <- file.path("../..", "baseline_PCs_relatedness", input,
  paste("dataProt_SS-205063.PC.related", input, "gprofiler.csv", sep = ".")
)
result <- lapply(fh, fread)
setattr(result, 'names', input)

# merge
result <- rbindlist(result, idcol = "dx")

# filter
result <- result[query %in% c("up", "down")]
result <- result[source %in% c("GO:BP", "REAC", "WP", "KEGG")]
setorder(result, dx, term_name, -query)
result[, logp := -log10(p_value)]
result[, x := frank(.(term_name), ties.method = "dense"), by = c("dx")]

# ready for plotting
a <- result[, .(dx, x, logp, term_name, query)]
a$dx <- droplevels(factor(a$dx, levels = input))
levels(a$dx) <- gsub("Plasma[_]*", "", levels(a$dx)) # drop Plasma label from plots
setorder(a, dx, x)

# plot func
CairoPDF(
  file = file.path(output, paste(nameset, "overlap_gprofiler.pdf", sep=".")),
  width = 10,
  height = 10
)
plot_circos_enr(a, 
  dx_col = "dx", 
  x_col = "x", 
  term_col = "term_name",
  sig_col = "logp",
  col_col = "up-down"
)
dev.off()

# All STRING enrich pathways -------------------------------------------
# DE file locations
fh <- file.path("../..", "baseline_PCs_relatedness", input, "STRING",
  "All_Overall_Enrich.csv"
)
result <- lapply(fh, fread)
setattr(result, 'names', input)

# merge
result <- rbindlist(result, idcol = "dx")

# filter
cats <- c(
  "GO Biological Process", 
  "GO Cellular Component", 
  "GO Molecular Function", 
  "Reactome Pathways"
)
result <- result[category %in% cats]
setorder(result, dx, description)
result[, logp := -log10(`FDR value`)]
result[, x := frank(.(description), ties.method = "dense"), by = c("dx")]

# ready for plotting
a <- result[, .(dx, x, logp, description, category)]
a$dx <- droplevels(factor(a$dx, levels = input))
levels(a$dx) <- gsub("Plasma[_]*", "", levels(a$dx)) # drop Plasma label from plots
setorder(a, dx, x)

# plot func
CairoPDF(
  file = file.path(output, paste(nameset, "overlap_STRING_overall.pdf", sep=".")),
  width = 20,
  height = 20
)
plot_circos_enr(a, 
  dx_col = "dx", 
  x_col = "x", 
  term_col = "description",
  sig_col = "logp",
  col_col = "category"
)
dev.off()


# STRING clusters annotated -------------------------------------------
# trying with just cluster name matching (imperfect)
# annotation file locations
fh <- file.path("../..", "baseline_PCs_relatedness", input, "STRING",
  "Annotated_clusters.csv"
)
result <- lapply(fh, fread)
setattr(result, 'names', input)

# merge
result <- rbindlist(result, idcol = "dx", fill = TRUE)

# annotation file locations
fh <- file.path("../..", "baseline_PCs_relatedness", input, "STRING",
  "Cluster_members.csv"
)
g_lists <- lapply(fh, fread)
setattr(g_lists, 'names', input)

# collapse & merge
g_lists <- lapply(g_lists, function(i) i[, .(symbol = list(symbol)), by = cluster])
g_lists <- rbindlist(g_lists, idcol = "dx")
result <- merge(result, g_lists, by = c("dx", "cluster"))

# fill sizes to reflect gene lists
result[, size := lengths(symbol)]
result[, logppi := -log10(PPI)]

# # filter
# cats <- c(
  # "GO Biological Process", 
  # "GO Cellular Component", 
  # "GO Molecular Function", 
  # "Reactome Pathways"
# )
# result <- result[category %in% cats]
a <- result[term != "Ambiguous"]
a$cluster <- factor(a$cluster, levels = mixedsort(unique(a$cluster)))
setorder(a, dx, cluster)
a[, x := frank(.(cluster), ties.method = "dense"), by = c("dx")]

# ready for plotting
a <- a[, .(dx, x, logppi, term, cluster)]
a$dx <- droplevels(factor(a$dx, levels = input))
levels(a$dx) <- gsub("Plasma[_]*", "", levels(a$dx)) # drop Plasma label from plots
setorder(a, dx, x)

# plot func
CairoPDF(
  file = file.path(output, paste(nameset, "overlap_STRING_clusters.pdf", sep=".")),
  width = 30,
  height = 30
)
plot_circos_enr(a, 
  dx_col = "dx", 
  x_col = "x", 
  term_col = "term",
  sig_col = "logppi",
  col_col = "dx",
  wrap = 40
)
dev.off()


# STRING clusters annotated -------------------------------------------
# trying with overlap enrichment using same result data as above
gs <- 6377 # size of all proteins

# Initialize results list for enrichment scores
go_table <- data.table(
  dx.1 = character(0),
  clust.1 = character(0),
  dx.2 = character(0),
  clust.2 = character(0),
  intersection = integer(0),
  pval = numeric(0),
  padj = numeric(0),
  logFDR = numeric(0)
)

# Nested loops for pairwise comparisons without repeating pairs
# # test use
# k <- input[1]
# j <- input[2]
for (k in 1:(length(input) - 1)) {
  for (j in (k + 1):length(input)) {
    # setup each dx with module gene lists
    x = result[dx == input[k]]
    y = result[dx == input[j]]
    go <- newGOM(
      sapply(
        X = seq_len(nrow(x)), 
        FUN = function(i) setNames(x$symbol[i], x$cluster[i])
      ),
      sapply(
        X = seq_len(nrow(y)), 
        FUN = function(i) setNames(y$symbol[i], y$cluster[i])
      ),
      genome.size = gs
    )

    # extract values and merge
    go_list <- lapply(c("intersection", "pval"), function(i) {
      melt(
        data = data.table(getMatrix(go, i), keep.rownames = TRUE), 
        id.vars = 'rn',
        value.name = i,
        variable.factor = FALSE
      )
    })
    go_list <- merge(go_list[[1]], go_list[[2]], by = c("rn", "variable"))
    go_list[, padj := p.adjust(pval, method = "fdr")]
    go_list[, logFDR := -log10(padj)]

    # add dx labels and rename reorder
    go_list[, dx.1 := input[k]]
    go_list[, dx.2 := input[j]]
    setnames(go_list, old = c("rn", "variable"), new = c("clust.1", "clust.2"))
    setcolorder(go_list, 
      neworder = c("dx.1", "clust.1", "dx.2", "clust.2", "intersection", "pval", "padj", "logFDR")
    )
    go_table <- rbindlist(list(go_table, go_list))
  }
}

# reorder at the end
go_table$clust.1 <- factor(go_table$clust.1, 
  levels = mixedsort(unique(go_table$clust.1))
)
go_table$clust.2 <- factor(go_table$clust.2, 
  levels = mixedsort(unique(go_table$clust.2))
)
setorder(go_table, dx.1, clust.1, dx.2, clust.2)

fwrite(go_table, 
  file = file.path(output, paste(nameset, "overlap_enrich_clusters.csv", sep="."))
)

# now need to convert into circos.link format and modify the function to pass it straight to it
# if no link table, separate function calculates it from matching terms

# Map values to colors
cols <- colorRampPalette(c("#ffffff", "#08306b"))
go_table[, tcols := cols(nrow(go_table))[as.numeric(cut(logFDR, breaks = nrow(go_table)))] ]
setorder(go_table, logFDR)
b <- go_table[padj < 0.05 & intersection > 2, .(dx.1, clust.1, dx.2, clust.2, tcols)]
b[, clust.1 := as.integer(gsub("M", "", clust.1))]
b[, clust.2 := as.integer(gsub("M", "", clust.2))]
b[, dx.1 := gsub("Plasma[_]*", "", dx.1)]
b[, dx.2 := gsub("Plasma[_]*", "", dx.2)]

# ready result table for plotting
a <- result
a[, x := as.integer(gsub("M", "", cluster))]
a$dx <- droplevels(factor(a$dx, levels = input))
levels(a$dx) <- gsub("Plasma[_]*", "", levels(a$dx)) # drop Plasma label from plots
a[, term_name := paste(cluster, term)]
setorder(a, dx, x)

# plot func
CairoPDF(
  file = file.path(output, paste(nameset, "overlap_enrich_clusters.pdf", sep=".")),
  width = 30,
  height = 30
)
plot_circos_enr(a, 
  dx_col = "dx", 
  x_col = "x", 
  term_col = "term_name",
  sig_col = "size",
  col_col = "dx",
  wrap = 40,
  mtype = "enrichment",
  mtable = b
)
dev.off()


# also need legend
# change bar plots to heatmaps of logppi and size 