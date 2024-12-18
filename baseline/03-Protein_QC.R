#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 03-Protein_QC.R
##
## Version: 0.0.4
##
## Purpose of script: Prepare LnPlasmaPTau data for gene set analyses
##
## Author: Robert B Butler III, Michael P Brown
##
## Date Created: 2022-12-21
##
## Copyright (c) 2022
## Email: rrbutler@stanford.edu, mpb@stanford.edu
##
## ---------------------------
##
## Notes:
##   Post linear regression, filters based on protein qualities, dropping
##   non-human and non-protein markers. Then matches to Ensembl gene IDs.
##
##   For spots detecting gene complexes it splits the genes to duplicate rows,
##   then for each gene, it selects the row with the highest effect size by
##   absolute value of beta. Interesting genes should be checked to see if they
##   were detected across multiple spots.
##
##   Finally, prints some plots, and calculates some lists of genes for
##   downstream analysis, including a gprofiler enrichment check using whole
##   SomaLogic set as custom background
##
##   REMEMBER, WHEN COMMENTING IN CODE, DO NOT REFERENCE INDIVIDUAL INFORMATION:
##   SAMPLES/PIDNS, DATES OR ANY INDIVIDUAL LEVEL INFO
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(data.table)
library(biomaRt)
library(ggplot2)
library(ggpubr)
library(gprofiler2)
library(patchwork)
library(qqman)
library(optparse)
library(tools)


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
    help = "File path to DE input file [default: %default]"
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
    default = "..",
    help = "File path to directory of reference data [default: %default]"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project [default: Nameset of input]"
  )
)

args_list <- parse_args(OptionParser(option_list = option_list))

# Nameset default case
if (is.null(args_list$nameset)) {
  args_list$nameset <- nameset_default(args_list$input)
}

# read in fields
fields <- c(
  "input",
  "output",
  "ref",
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


ensver <- "112"


# Functions -------------------------------------
#' for a given column of NCBI IDs, replace them with ensembl gene ids
#' includes a filter for autosomal genes that are not pseudo or TEC
#'
#' @param dt data table to replace names
#' @param colname string, name of the column containing NCBI gene IDs
#' @param symname string, name of the column containing gene Symbols 
#' @param keep.symbols boolean, retain the symbols column?
#'
#' @return dt with a GENE column
convertGeneSymbol <- function(dt, colname, symname, keep.symbols = FALSE) {
  require("biomaRt")
  print(paste(length(unique(dt[[colname]])), "genes to look up"))
  ensembl <- useEnsembl("ensembl", 
    dataset = "hsapiens_gene_ensembl", 
    version = ensver
  )
  # match with ensembl_gene_id across both
  x <- data.table(getBM(
    attributes = c(
      "ensembl_gene_id", 
      "external_gene_name",
      "entrezgene_id",
      "gene_biotype", 
      "chromosome_name"
    ),
    filters = "entrezgene_id",
    values = unique(dt[[colname]]),
    mart = ensembl
  ))
  print(paste(nrow(x), "total matches returned"))
  # ensembl transcript types to keep, excluding pseudogenes and TECs
  keep_categ <- grep("pseudo", unique(x$gene_biotype), value = TRUE, invert = TRUE)
  keep_categ <- grep("TEC", keep_categ, value = TRUE, invert = TRUE)
  x <- x[gene_biotype %in% keep_categ]
  print(paste(nrow(x), "genes correct biotype"))
  # drop duplicates on patches and scaffolds if exist on main chromosomes
  x[, chr_length := nchar(chromosome_name)]
  x <- unique(
    x[order(entrezgene_id, external_gene_name, chr_length)], 
    by = c("entrezgene_id", "external_gene_name")
  )
  setnames(x, 
    old = c("ensembl_gene_id", "entrezgene_id", "external_gene_name"), 
    new = c("GENE", "NCBI_ID", "Symbol")
  )
  # replace gene names
  dt <- merge(
    x[, .(GENE, NCBI_ID, Symbol)], 
    dt, 
    all.y = TRUE, 
    by.x = c("NCBI_ID", "Symbol"),
    by.y = c(colname, symname)
  )

  # protein complexes will not match by symbol, need to update them separately
  # convert Symbol to list and match by NCBI and any symbol in list
  cols <- names(dt)[names(dt) != "Symbol"]
  dt[, Symbol := as.list(strsplit(Symbol, "\\|"))]
  dt[UniProt %like% "\\|", GENE := mapply(function(z, y) {
    return(x[NCBI_ID == z & Symbol %in% y, GENE])
  }, NCBI_ID, Symbol), by = .I]
  # collapse Symbol back to a character
  dt[, Symbol := paste(unlist(Symbol), collapse = "|"), by = .I]
  dt$Symbol <- unlist(dt$Symbol) # for the life of me this unlist will not work inside the line above

  # write to file
  print(paste(length(unique(dt[["GENE"]])), "unique genes with ensembl ids"))
  if (keep.symbols == FALSE) dt[, NCBI_ID := NULL]
  return(unique(dt))
}

#' for a string of ens ids: explode, look up symbols, return symbols as string
#'
#' @param string comma separated ens_ids
#' @param dds DESeq Dataset with annotations in mcol
#' @param col string, colname to grab
#'
#' @return string of concatenated, comma separated symbols
get_mcols <- function(str, dds, col) {
  require("stringr")
  y <- dds[[col]][match(unlist(str_split(str, ",")), dds[["GENE"]])]
  return(paste0(y, collapse = ","))
}


# Main ------------------------------------------
# Load assocaition table
a <- fread(input)
a <- a[ix != ""]
# a$variable = gsub("^X\\.", "\\.", a$variable)

# map to reference key
b <- fread(ref)
b <- b[Organism == "Human" & Type == "Protein"]
b[, c("Organism", "Units", "Type", "EntrezGeneSymbol_v2") := NULL]
c <- merge(a, b, by.x = "ix", by.y = "Key_2")

# padjust and clean up (total set of human protein probes)
c[, padj := p.adjust(pval, method = "fdr")]
c[, pbonf := p.adjust(pval, method = "bonferroni")]
c <- c[EntrezGeneSymbol != ""]

# fix piped rows of multiple NCBI IDs (explode to multiple rows)
c$EntrezGeneID <- gsub("\\s+", "|", c$EntrezGeneID, perl = TRUE)
# c$EntrezGeneID = sub("\\|.*$", "", c$EntrezGeneID) # take only the first NCBI_ID in the row
cols <- names(c)[names(c) != "EntrezGeneID"]
c <- c[, .(EntrezGeneID = tstrsplit(EntrezGeneID, "\\|")), by = cols]
c$EntrezGeneID <- as.integer(c$EntrezGeneID)

# drop lower duplicates, use absolute value of beta
# c[duplicated(EntrezGeneSymbol) | duplicated(EntrezGeneSymbol, fromLast=T)][order(EntrezGeneSymbol, -beta)]
c <- c[order(-abs(beta))][, head(.SD, 1), by = "EntrezGeneID"]

# lookup NCBI IDs for matching
c <- convertGeneSymbol(c, 
  colname = "EntrezGeneID", 
  symname = "EntrezGeneSymbol", 
  keep.symbols = TRUE
)

fwrite(c, file = file.path(output, paste(nameset, "csv", sep = ".")))


# plot histograms ----------------------------------------------
# all betas
p1 <- ggplot(c, aes(x = beta)) +
  geom_density(fill = "#999999") +
  ggtitle("All proteins") +
  theme_classic() +
  xlim(-0.5, 0.5) +
  theme(plot.title = element_text(hjust = 0.5))

# calculate genomic inflation
lambda.c <- median(qchisq(c$pval, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
# plot QQ
p3 <- ggplot(c, aes(sample = -log10(pval))) +
  stat_qq(distribution = qexp, dparams = list(rate = log(10))) +
  geom_abline(col = "red", intercept = 0, slope = 1) +
  ggtitle(sprintf("QQ (Lambda: %.3f)", lambda.c)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

# check for significance before writing genesets, pathway enrich, etc
if (nrow(c[padj < 0.05]) < 2) {
  # write plots
  ggsave(
    plot = p1 + p3,
    filename = paste(nameset, "beta_distributions.pdf", sep = "."),
    path = output,
    width = 10,
    height = 6,
    device = cairo_pdf
  )

  # exit early
  stop("fewer than 2 genes are padj < 0.05, exiting...", call. = FALSE)
}

# Bonferroni significant betas
c_cut <- c[padj < 0.05]
# make density distribution for quantiles
dens <- density(c_cut$beta)
qs <- quantile(c_cut$beta, c(0.2, 0.8))
dt <- data.frame(x = dens$x, y = dens$y)
dt$quant <- factor(
  findInterval(dt$x, c(min(qs[[1]], -.Machine$double.xmin), 0, qs[[2]]))
)
# plot
p2 <- ggplot(dt, aes(x, y)) +
  geom_line() +
  labs(x = "beta", y = "density", title = "Proteins with Padj < 0.05") +
  geom_ribbon(aes(ymin = 0, ymax = y, fill = quant)) +
  scale_fill_brewer(
    type = "div", 
    palette = "RdYlBu", 
    direction = -1, 
    name = "",
    labels = c("Bottom 10%", "Less than 0", "Greater than 0", "Top 10%")
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.position = c(0.9, 0.5)
  )

# write plots
ggsave(
  plot = p2 +
    inset_element(p1, left = 0.05, top = 1, bottom = 0.6, right = 0.3) +
    inset_element(p3, left = 0.7, top = 1, bottom = 0.6, right = 0.95),
  filename = paste(nameset, "beta_distributions.pdf", sep = "."),
  path = output, 
  width = 10, 
  height = 6, 
  device = cairo_pdf
)

# write gene lists
sets <- list(
  all = c[!is.na(GENE) & padj < 0.05, .(GENE, beta)][order(-abs(beta))],
  up = c[!is.na(GENE) & padj < 0.05 & beta > 0, .(GENE, beta)][order(-beta)],
  down = c[!is.na(GENE) & padj < 0.05 & beta < 0, .(GENE, beta)][order(beta)],
  top20 = c[!is.na(GENE) & padj < 0.05 & beta > qs[[2]], .(GENE, beta)
            ][order(-beta)],
  bottom20 = c[!is.na(GENE) & padj < 0.05 & beta < qs[[1]], .(GENE, beta)
               ][order(beta)],
  top.bonf = c[!is.na(GENE) & pbonf < 0.05 & beta > 0, .(GENE, beta)
               ][order(-beta)],
  bottom.bonf = c[!is.na(GENE) & pbonf < 0.05 & beta < 0, .(GENE, beta)
                  ][order(beta)]
)
lapply(seq(sets), function(i) {
  fwrite(sets[[i]],
    file = file.path(output, paste(nameset, names(sets)[i], "txt", sep = ".")),
    sep = "\t", 
    col.names = FALSE
  )
})

# gprofiler enrichment
g <- gost(
  query = sets, 
  organism = "hsapiens", 
  ordered_query = FALSE, 
  evcodes = TRUE,
  custom_bg = unique(c[!is.na(GENE), GENE])
)
# write to file
if (is.list(g$result)) {
  # get symbols for terms
  g$result$symbols <- lapply(g$result$intersection, get_mcols, c, "Symbol")
  # write to file
  fwrite(g$result, file = file.path(output, paste0(nameset, ".gprofiler.csv")))
}
