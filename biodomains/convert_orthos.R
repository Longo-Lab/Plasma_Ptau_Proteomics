#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: convert_orthos.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Map orthologous EnsemblIDs, Symbols across species
#'
#' Author: Robert R Butler III
#'
#' Date Created: 2024-06-13
#'
#' Copyright (c) 2024
#' Email: rrbutler@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'   Breaking out ortho functions from 01-all_gene_lists.R
#'
#'   Now handles keeping the biomart lookup in memory
#'     if looking back and forth across species, set store_dt = FALSE to get the
#'     full lookup in each direction. (ie if multiple mouse -> human lookups, 
#'     storing the table is useful, but if doing both mouse -> human and 
#'     human -> mouse lookups, you will end up getting mouse -> human -> mouse 
#'     on the second, which will lose some terms unnecessarily
#'
#'   Updated to handle any conversion in a single function: convert_orthos
#'    prevents the need to store both tables.
#'   Will handle one2many orthologs by keeping the top othro based on ensembl criterion
#'    sans "tree compliant" Min GOC > Min WGA > Min %identity
#'    http://mart.ensembl.org/info/genome/compara/Ortholog_qc_manual.html
#'
#'
#' ---------------------------

library(data.table)
library(biomaRt)
library(readxl)
library(tidyverse)

# wk_dir <- "/labs/flongo/reference/biodomains"
# ensembl_ver <- "111"


# Functions -------------------------------------------------------------

#' Filter genes table by gene_biotype and chromosome_name columns
#'
#' @param dt1 Data table fetched from Ensembl containing gene_biotype and chromosome_name columns
#' @param species string to define number of chromosomes
#'
#' @return Filtered dt1
filter_genes <- function(dt1, species = c("mmusculus", "hsapiens")) {
  rlang::arg_match(species)
  # Ensembl transcript types to keep, excluding pseudogenes and TECs
  keep_categ <- grep("pseudo", unique(dt1$gene_biotype), value = TRUE, invert = TRUE)
  keep_categ <- grep("TEC", keep_categ, value = TRUE, invert = TRUE)
  dt1 <- dt1[gene_biotype %in% keep_categ]
  # print(paste(nrow(dt1), "genes correct biotype"))
  
  if (species == "mmusculus") {
    # Use only canonical chromosomes
    use_chrs <- c(1:19, "X", "Y", "MT")
    dt1 <- dt1[chromosome_name %in% use_chrs]
    # print(paste(nrow(dt1), "genes in chr 1:19, X, Y, MT"))
  } else if (species == "hsapiens") {
    # Use only canonical chromosomes
    use_chrs <- c(1:22, "X", "Y", "MT")
    dt1 <- dt1[chromosome_name %in% use_chrs]
    # print(paste(nrow(dt1), "genes in chr 1:22, X, Y, MT"))
  } else {
    stop(
      "incorrect species in gene filter, hsapiens or mmusculus...", 
      call. = FALSE
    )
  }
  
  dt1
}

#' Handle input column name for outfile
#'   remove or keep, if it has the same name as output column, rename to origID
#'
#' @param dt1 Data table
#' @param keep_original Boolean, retain the original ids in origID column?
#' @param in_type string, type of input column (current name of input column) 
#' @param in_col string, original name of input column
#' @param out_col string, name of output column
#'
#' @return data table with renamed or dropped input col
rename_in_col <- function(dt1, keep_original = FALSE, in_type, in_col, out_col){
  if (keep_original && in_col == out_col) {
    setnames(dt1, in_type, "origID")
  } else if (keep_original && in_col != out_col) {
    setnames(dt1, in_type, in_col)
  } else {
    dt1[, c(in_type) := NULL]
  }
  return(dt1)
}


#' Map species gene IDs/symbols to either human or mouse, optionally keeping original IDs
#'
#' @param dt1 Data table to add human Ensembl ID column based on gene symbol
#' @param in_col Name of the column containing id to convert
#' @param in_species string, hsapiens(human) or mmusculus(mouse)
#' @param in_type string, input id is either "ensID" or "symbol"
#' @param out_col Name of the column containing output IDs
#' @param out_species string, hsapiens(human) or mmusculus(mouse)
#' @param out_type string, output id is either "ensID" or "symbol"
#' @param ortho_type string, how to filter orthologs, either "top" or "one2one"
#' @param filter Boolean, filter by output biotype and chromosome?
#' @param keep_original Boolean, retain the original ids in origID column?
#' @param ver string, either a specific ensembl version, or NULL for latest
#' @param store_dt Boolean, keep mapping table as "gene_info" variable for script
#'
#' @return dt1 with added out_type column
convert_orthos <- function(dt1, 
                           in_col = "GENE",
                           in_species = c("mmusculus", "hsapiens"),
                           in_type = c("ensID", "symbol"), 
                           out_col = "GENE",
                           out_species = c("mmusculus", "hsapiens"),
                           out_type = c("ensID", "symbol"), 
                           ortho_type = c("top", "one2one"),
                           filter = TRUE, 
                           keep_original = FALSE,
                           ver = NULL,
                           store_dt = FALSE) {
  require("biomaRt")
  in_species <- rlang::arg_match(in_species)
  in_type <- rlang::arg_match(in_type)
  out_species <- rlang::arg_match(out_species)
  out_type <- rlang::arg_match(out_type)
  ortho_type <- rlang::arg_match(ortho_type)
  # lapply(c(in_species, in_type, out_species, out_type, ortho_type), rlang::arg_match)
  if (in_species == out_species && in_type == out_type) {
    print("Same species and same ID type, returning original data")
    return(dt1)
  }
  
  # connect to mart
  mart <- useEnsembl(
    biomart = "genes", 
    dataset = paste0(in_species, "_gene_ensembl"),
    version = ver
  )

  # May be many-to-many join (each human ensembl ID can map to multiple gene symbol and belong to multiple modules)
  # Handling lookups as either same species (simple), existing (in memory), or else full lookup
  if (in_species == out_species) {
    print(sprintf("Same species, looking up %s given %s", out_type, in_type))
    
    # build homolog lookup table
    x <- data.table(
      getBM(mart = mart, 
        attributes = c(
          "ensembl_gene_id", 
          "external_gene_name", 
          "chromosome_name", 
          "gene_biotype"
        )
      )
    )
    
    # rename and drop missing
    setnames(x, 
      old = c("ensembl_gene_id", "external_gene_name"), 
      new = c("ensID", "symbol")
    )
    x <- x[ensID != ""]
    x <- x[symbol != ""]

    # filter based on gene biotype and chromosome
    if (filter) {
      x <- filter_genes(x, species = in_species)
    }
    
    # merge gene names
    old_cols <- colnames(dt1)[colnames(dt1) != in_col]
    x_cols <- c("ensID", "symbol")
    dt1 <- merge(x[, ..x_cols], dt1, all.y = TRUE, by.x = in_type, by.y = in_col)
    setnames(dt1, ifelse(out_type == "symbol", "symbol", "ensID"), out_col)
    dt1 <- dt1[!is.na(get(out_col))]
    print(
      sprintf(
        "%d unique %ss found for %d %ss",
        length(unique(dt1[[out_col]])), 
        out_type, 
        length(unique(dt1[[in_type]])), 
        in_type
      )
    )
    
    # set column names for detection of duplicates and output table
    out_cols <- union(c(out_col, in_type), old_cols)
    check_cols <- union(in_type, old_cols)
    setorderv(dt1, cols = out_col)
    if (anyDuplicated(dt1, by = check_cols) > 0) {
      cat(
        "Warning: duplicate IDs found in the lookup matches...\n",
        "These could be duplicate rows in your table, or issues with\n",
        "database annotation keeping only the lowest ID (alphanumerically)\n",
        "if keep_original = TRUE, dups will be kept with the original column\n"
      )
    }
    if (!keep_original) {
      dt1 <- dt1[, head(.SD, 1), by = check_cols]
    }
    # prepare column names for outfile
    dt1 <- rename_in_col(dt1[, ..out_cols], 
      keep_original = keep_original, 
      in_type = in_type, 
      in_col = in_col, 
      out_col = out_col
    )
    return(unique(dt1))
  } else if (exists("gene_info") && is.data.table(gene_info)) {
    x <- gene_info
  } else {
    print(
      sprintf(
        "Looking up %s %s given %s %s", 
        out_species, 
        out_type, 
        in_species, 
        in_type
      )
    )
    # build homolog lookup table
    x <- data.table(
      getBM(mart = mart, 
        attributes = c(
          "ensembl_gene_id", 
          "external_gene_name", 
          "chromosome_name", 
          paste0(out_species, "_homolog_ensembl_gene"), 
          paste0(out_species, "_homolog_associated_gene_name"), 
          paste0(out_species, "_homolog_chromosome"),
          paste0(out_species, "_homolog_goc_score"), 
          paste0(out_species, "_homolog_orthology_type"), 
          paste0(out_species, "_homolog_perc_id"), 
          paste0(out_species, "_homolog_perc_id_r1"), 
          paste0(out_species, "_homolog_wga_coverage")
        )
      )
    )
    
    # rename and drop missing
    setnames(x, 
      old = c(
        "ensembl_gene_id",
        "external_gene_name",
        paste0(out_species, "_homolog_ensembl_gene"), 
        paste0(out_species, "_homolog_associated_gene_name"), 
        paste0(out_species, "_homolog_chromosome"),
        paste0(out_species, "_homolog_goc_score"),
        paste0(out_species, "_homolog_orthology_type"), 
        paste0(out_species, "_homolog_perc_id"),
        paste0(out_species, "_homolog_perc_id_r1"), 
        paste0(out_species, "_homolog_wga_coverage")
      ), 
      new = c(
        "ensID",
        "symbol",
        "h_ensID",
        "h_symbol",
        "h_chromosome",
        "goc",
        "orthology_type",
        "perc_id",
        "perc_id_r1", 
        "wga"
      )
    )
    x <- x[h_ensID != ""]
    x <- x[h_symbol != ""]

    # define min %identity
    x[, min_perc_id := min(perc_id, perc_id_r1, na.rm = TRUE), by = .I]

    # one2many orthologs are kept based on top hit by Min GOC > Min WGA > Min %identity or else one2one only
    setorderv(x, 
      cols = c(in_type, "goc", "wga", "min_perc_id"), 
      order = c(1, -1, -1, -1), 
      na.last = TRUE
    )
    if (ortho_type == "one2one") {
      x <- x[orthology_type == "ortholog_one2one"]
    } else {
      x <- x[, head(.SD, 1), by = in_type]
    }

    # filter based on gene biotype and chromosome
    if (filter) {
      x[, chromosome_name := h_chromosome]
      y <- data.table(
        getBM(mart = mart, 
          filter = "ensembl_gene_id",
          values = x$ensID,
          attributes = c("ensembl_gene_id", "gene_biotype")
        )
      )
      setnames(y, "ensembl_gene_id", "ensID")
      x <- filter_genes(merge(x, y), species = out_species)
    }

    # store file if script is still running
    if (store_dt) {
      assign("gene_info", x, envir = .GlobalEnv)
    }
  }

  # merge gene names
  old_cols <- colnames(dt1)[colnames(dt1) != in_col]
  x_cols <- c(in_type, "h_ensID", "h_symbol", "goc", "wga", "min_perc_id")
  dt1 <- merge(x[, ..x_cols], dt1, all.y = TRUE, by.x = in_type, by.y = in_col)
  setnames(dt1, ifelse(out_type == "symbol", "h_symbol", "h_ensID"), out_col)
  dt1 <- dt1[!is.na(get(out_col))]
  print(
    sprintf(
      "%d unique %s %ss found for %d %s %ss using %s matching",
      length(unique(dt1[[out_col]])), 
      out_species,
      out_type, 
      length(unique(dt1[[in_type]])), 
      in_species,
      in_type,
      ortho_type
    )
  )
  
  # set column names for detection of duplicates and output table
  out_cols <- union(c(out_col, in_type), old_cols)
  check_cols <- union(out_col, old_cols)
  
  # resolve many2one issues in the same manner as one2many above 
  # (check other columns of dt1 to make sure duplicates aren't purposeful)
  setorderv(dt1, 
    cols = c(out_col, "goc", "wga", "min_perc_id"), 
    order = c(1, -1, -1, -1), 
    na.last = TRUE
  )
  if (ortho_type == "one2one" & anyDuplicated(dt1, by = check_cols) > 0) {
    cat(
      "Warning: some many2one duplicates found in the one2one matches...\n",
      "These could be duplicate names in your table,\n",
      "or issues with orthology annotation\n"
    )
    cat(
      "  trimming to top hit for each ortholog and : ", 
      unique(dt1[duplicated(dt1, by = check_cols), get(out_col)]),
      "\n"
    )
  }
  dt1 <- dt1[, head(.SD, 1), by = check_cols]
  
  # prepare column names for outfile
  dt1 <- rename_in_col(dt1[, ..out_cols], 
    keep_original = keep_original, 
    in_type = in_type, 
    in_col = in_col, 
    out_col = out_col
  )
  return(unique(dt1))
}


