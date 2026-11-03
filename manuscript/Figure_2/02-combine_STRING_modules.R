#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 02-combine_STRING_modules.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Combine STRING module membership (Cluster_members.csv)
#'  and module annotation (Annotated_clusters.csv) across biomarkers into a
#'  single Supplementary Table for Fig 2
#'
#' Author: Robert B Butler III
#'
#' Date Created: 2026-07-29
#'
#' Copyright (c) 2026
#' Email: rrbutler@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'   v0.0.1 - initial build
#'   Reads Cluster_members.csv + Annotated_clusters.csv from each biomarker's
#'   baseline_PCs_relatedness/<Biomarker>/STRING folder (.APOE-suffixed
#'   folders excluded), merges by `cluster`, and stacks across biomarkers
#'   with rbindlist(idcol = "Biomarker").
#'
#' ---------------------------
#' load up the packages we will need:  (uncomment as required)
library(data.table)

base_dir <- "../../baseline_PCs_relatedness"

# Biomarker subfolders that contain a STRING module analysis.
# (Folders with the ".APOE" suffix are excluded.)
biomarkers <- c(
  "PlasmaGFAP",
  "PlasmaNfL",
  "PlasmapTau217_AB42Ratio",
  "PlasmaABRatio",
  "PlasmaPTau181",
  "PlasmaAB140P",
  "PlasmaAB142P",
  "PlasmapTau217"
)

names(biomarkers) <- biomarkers

# For each biomarker: read cluster membership + cluster annotation,
# and merge them on `cluster`.
mod.list <- lapply(biomarkers, function(bm) {

  string_dir <- file.path(base_dir, bm, "STRING")

  members <- fread(file.path(string_dir, "Cluster_members.csv"))
  clusters <- fread(file.path(string_dir, "Annotated_clusters.csv"))

  merge(members, clusters, by = "cluster", all.x = TRUE)
})

## Stack all biomarkers into one table, tagging rows with their
## biomarker of origin.
supp_table <- rbindlist(mod.list, idcol = "Biomarker", fill = TRUE)
# fill in missing cluster sizes
supp_table[, size := length(size), by = .(Biomarker, cluster)]

## Tidy column order/names for the supplementary table.
setnames(supp_table,
         old = c("cluster", "symbol", "term", "PPI"),
         new = c("Cluster", "Symbol", "Term", "PPI"),
         skip_absent = TRUE)

setcolorder(supp_table,
            intersect(c("Biomarker", "Cluster", "Symbol", "Term",
                        "PPI_pvalue", "size"),
                      names(supp_table)))

out_dir <- file.path(".")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

fwrite(supp_table,
       file.path(out_dir, "SupplementaryTable_STRING_modules.csv"))
