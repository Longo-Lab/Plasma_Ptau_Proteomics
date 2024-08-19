#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 12-Combined-Enrich.R
#'
#' Version: 0.0.4
#'
#' Purpose of script: Show Fisher's Exact Enrichment of proteomic modules
#' based on user input
#'
#' Author: Robert B Butler III, Michael Peter Brown
#'
#' Date Created: 2024-7-8
#'
#' Copyright (c) 2022
#' Email: rrbutler@stanford.edu, mpb@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'   v0.0.4 - added new option and extra error handling 
#'   v0.0.4 - added output path to each instance of fwrite/pdf
#'   v0.0.3 - added subset to reference genes that are only in SomaLogic spots
#'   v0.0.3 - run styler
#'   v0.0.3 - added command line support for biodomain folder
#'   v0.0.3 - added default options
#'   v0.0.2 - added command line support
#'   v0.0.2 - sourced external functions for heatmaps
#'   v0.0.2 - consolidated enrichment to single script
#'
#'
#' ---------------------------
#' load up the packages we will need:  (uncomment as required)
library(optparse)
library(ggplot2)
library(ggpubr)
library(data.table)
library(stringr)
library(plyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(GeneOverlap)
library(tools)

# Nameset helper function
nameset_default <- function(input){
  file_path_sans_ext(basename(input))
}

# read in args
option_list <- list(
  make_option(
    c("-m", "--modules"),
    type = "character",
    help = "Name of module script set desired to run"
  ),
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
    default = "../biodomains",
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
if(is.null(args_list$nameset)){
  args_list$nameset <- nameset_default(args_list$input)
}

# read in fields
fields <- c(
  "modules",
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

# Source() file path
source_file <- file.path(ref, paste0("biodomain_", modules, ".R"))

# Source() error/warning handling
tryCatch(
  expr = {
    source(source_file)
  },
  error = function(e){
    message('Error with source file:')
    print(e)
  },
  warning = function(w){
    message('Warning messages from source file:')
    print(w)
  }
)

# Functions ---------------------------
gommat_2_dt <- function(go, dt, dtype = c(
                          "pval", "odds.ratio",
                          "intersection", "union", "Jaccard"
                        )) {
  # add gom object to a data.table
  # gom.obj: GeneOverlap object
  # dt: data.table to add to
  # dtype: which GO measure to extract
  d <- melt(data.table(getMatrix(go, name = dtype), keep.rownames = "sample"),
    variable.name = "group", value.name = dtype, id.vars = "sample"
  )
  dt <- merge(dt, d, by = c("sample", "group"))
  return(dt)
}

get_mcols <- function(str, result, col) {
  # for a string of ens ids: explode, look up symbols, return symbols as string
  # string: comma separated ens_ids
  # result: data table with ens_ids in Gene_id
  # col: colname to grab
  require("stringr")
  y <- result[[col]][match(unlist(str_split(str, ",")), result$GENE)]
  return(paste0(y, collapse = ","))
}

gomtable <- function(list1, list2, background, outname, width = 7,
                     height = 7, ...) {
  # Takes two gene overlap lists and writes a table of stats, makes list object
  # list1: named list object
  # list2: named list object
  # outname: name nameset for file
  # width: width of pdf files to be plotted in inches*2
  # height: height of pdf files
  # returns a list object of table, r_names, c_names
  # enrichment basic
  gom.obj <- newGOM(list1, list2, background)
  r_names <- names(list1)
  c_names <- names(list2)
  # # geneoverlap heatmap
  # pdf(file=paste(nameset, nameset, outname, "pval.pdf", sep="."),
  # width=width*2, height=height)
  # drawHeatmap(gom.obj, grid.col="Blues", cutoff=1)
  # dev.off()
  # write DE vals to file
  b <- melt(data.table(getMatrix(gom.obj, name = "pval"), keep.rownames = "sample"),
    variable.name = "group", value.name = "P", id.vars = "sample"
  )
  # write dtypes to table
  b <- gommat_2_dt(gom.obj, b, dtype = "odds.ratio")
  b <- gommat_2_dt(gom.obj, b, dtype = "Jaccard")
  b <- gommat_2_dt(gom.obj, b, dtype = "intersection")
  b <- gommat_2_dt(gom.obj, b, dtype = "union")
  # add fdr
  b[, FDR := p.adjust(P, method = "fdr")]
  b[, LOGP := -log10(P)]
  b[, LOGFDR := -log10(FDR)]
  # get intersection genes
  nest <- getNestedList(gom.obj, name = "intersection")
  nest2 <- lapply(seq(nest), function(i) {
    lapply(nest[[i]], function(j) paste0(j, collapse = ","))
  })
  nest3 <- lapply(seq(nest2), function(i) stack(nest2[[i]]))
  setattr(nest3, "names", names(nest))
  dt <- rbindlist(nest3, idcol = T)
  setnames(dt, c("group", "ens_id", "sample"))
  dt$symbol <- lapply(dt$ens_id, get_mcols, result, "EntrezGeneSymbol")
  b <- merge(b, dt, by = c("sample", "group"), all.x = T, sort = F)
  # write and return
  fwrite(b,
    file = file.path(output, paste(nameset, modules, outname, "txt", sep = ".")), sep = "\t",
    quote = F
  )
  return(list("b" = b, "r_names" = r_names, "c_names" = c_names))
}

enrich <- function(list1, list2, gs, outname, width = 7, height = 7, fdr = TRUE, ...) {
  # function runs enrichment and saves to file for odds ratio and padj/p
  # takes two set lists, background gene count and output name and pdf width
  # DE enrich
  de_enrich <- gomtable(list1, list2, gs, outname, width = width, height = height)

  # fetch Jaccard index
  f <- as.matrix(dcast(de_enrich$b, sample ~ group, value.var = "Jaccard"),
    rownames = "sample"
  )
  f <- f[de_enrich$r_names, de_enrich$c_names, drop = F]
  # log10p heatmap
  if (isTRUE(fdr)) {
    e <- as.matrix(dcast(de_enrich$b, sample ~ group, value.var = "LOGFDR"),
      rownames = "sample"
    )
    e <- e[de_enrich$r_names, de_enrich$c_names, drop = F]
    # filename and heatmap
    pdf( 
      file = file.path(output, paste(nameset, modules, outname, "logfdr.pdf", sep = ".")), width = width,
      height = height
    )
    heat(e, e, f, "-log FDR", ...)
    dev.off()
  } else {
    e <- as.matrix(dcast(de_enrich$b, sample ~ group, value.var = "LOGP"),
      rownames = "sample"
    )
    e <- e[de_enrich$r_names, de_enrich$c_names, drop = F]
    # filename and heatmap
    pdf(
      file = file.path(output, paste(nameset, modules, outname, "logp.pdf", sep = ".")), width = width,
      height = height
    )
    heat(e, e, f, "-log P", ...)
    dev.off()
  }
  # OR heatmap
  c <- as.matrix(dcast(de_enrich$b, sample ~ group, value.var = "odds.ratio"),
    rownames = "sample"
  )
  c <- c[de_enrich$r_names, de_enrich$c_names, drop = F]
  # filename and heatmap
  pdf(
    file = file.path(output, paste(nameset, modules, outname, "logOR.pdf", sep = ".")), width = width,
    height = height
  )
  heat(c, e, f, "Odds Ratio", ...)
  dev.off()
}


# Main ------------------------------------------
# DE file locations
result <- fread(input)
f_names <- sub(paste0(".", modules), "", modules)
replacement <- f_names

# Fisher exact testing for enrichment of DE sets -----------------
# generate list of gene sets for AD genes
# DE genes (all)
de_list <- lapply(f_names, function(i) {
  result[padj < 0.05, GENE] # geneoverlap errors with only one list, so run it with the up/down
})
# names(de_list) = c(replacement, paste0(replacement, "_"))

# DE genes (split up/down)
up_list <- lapply(f_names, function(i) {
  result[padj < 0.05 & beta > 0, GENE]
})
down_list <- lapply(f_names, function(i) {
  result[padj < 0.05 & beta < 0, GENE]
})

# DE genes (split up/down)
qs <- quantile(result[padj < 0.05, beta], c(0.1, 0.9))
top10_list <- lapply(f_names, function(i) {
  result[padj < 0.05 & beta > qs[[2]], GENE]
})
bottom10_list <- lapply(f_names, function(i) {
  result[padj < 0.05 & beta < qs[[1]], GENE]
})

# join and name
updown_list <- c(top10_list, up_list, de_list, down_list, bottom10_list)
f_names_updown <- c("Top 10%", "Up", "All", "Down", "Bottom 10%")
names(updown_list) <- paste(replacement, f_names_updown)


# Enrichments -----------------------------------
ad_list <- get_modules_list(ref)

# filter to only genes with protein spots
ad_list <- lapply(ad_list, function(i) { i[i %in% result[, GENE]] })

# background DE genes list
gs <- length(unique(c(result[, GENE]))) # only the Somalogic set are detectable

# enrich for up/down genes
enrich(
  list1 = updown_list, 
  list2 = ad_list, 
  gs = gs, 
  outname = "Modules_Up-Down", 
  width = w, 
  height = h1, 
  path = ref
)
