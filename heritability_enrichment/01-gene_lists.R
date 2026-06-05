#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 01-gene_lists.R
##
## Version: 0.0.2
##
## Purpose of script: Generate DE genesets for enrichment analysis 
##
## Author: Robert B Butler III
##
## Date Created: 2024-01-17
##
## Copyright (c) 2024
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##   Currently hardcodes the input files (good target for options funtion)
##
##     Looks at DE gene sets, drug comparison (in stim), stim comparison, wgcna modules
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(data.table)
library(optparse)
# library(biomaRt)

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
    c("-r", "--ref"),
    default = "ref/ENSG_coord_vGRCh37_PC.GRCh37.txt",
    help = "Full path to reference genome coors file[default %default]"
  ),
  make_option(
    c("-c", "--covs"),
    default = "hla,apoe",
    help = "covariates to filter for analysis [default %default]"
  ),
  make_option(
    c("-g", "--gwas"),
    default = "ref",
    help = "Full path to gwas data directory [default %default]"
  )
)

args_list <- parse_args(OptionParser(option_list = option_list))

# read in fields
fields <- c(
  "input",
  "output",
  "ref",
  "covs",
  "gwas"
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
covs <- strsplit(covs, ",")[[1]]


# Functions ------------------------------------------------------------------


# Main -----------------------------------------------------------------------
# DE file locations
fh <- file.path("..", "baseline_PCs_relatedness", input,
  paste("dataProt_SS-205063.PC.related", input, "csv", sep = ".")
)
result <- lapply(fh, fread)
setattr(result, 'names', input)

# merge
result <- rbindlist(result)
setnames(result, "t stat", "t_stat")

# Write DE sets ---------------------------------------------------------

# DE genes
# generate list of gene sets for AD genes
lapply(input, function(i) {
  fwrite(result[dx == i & padj < 0.05, "GENE"], 
    quote = FALSE, 
    col.names = FALSE, 
    file = paste0("padj_", i, ".txt")
  )
})


# top 10% genes
# generate list of gene sets for AD genes
lapply(input, function(i) {
  fwrite(result[dx == i
                ][order(-t_stat)
                  ][t_stat > quantile(t_stat, probs = 0.8)
                    ][, "GENE"], 
    quote = FALSE, 
    col.names = FALSE, 
    file = paste0("top20_", i, ".txt")
  )
})

# # DE genes (split up/down)
# # up
# lapply(replace, function(i) {
  # fwrite(result[Gene_set == i & padj < 0.05 & log2FoldChange > 0, "GENE"],
         # quote=F, col.names=F, file=paste0(i, ".up.txt"))
# })

# # down
# lapply(replace, function(i) {
  # fwrite(result[Gene_set == i & padj < 0.05 & log2FoldChange < 0, "GENE"],
         # quote=F, col.names=F, file=paste0(i, ".down.txt"))
# })

# linear correlations -----------------------------------------------------
# Pvals only for analysis convert to -log(P) per magma manual
result[, LOGP := -log10(pval) ]

# Pvals with L2FC sign
result[, signLOGP := sign(beta)*LOGP]

# Pvals with beta & LOGP
result[, betaLOGP := beta*LOGP]

# drop those missing a ENSG
result <- result[GENE != ""]

lapply(c("betaLOGP", "signLOGP", "t_stat"), function(i) {
  fwrite(dcast(result[], GENE ~ dx, value.var = i), 
    quote = FALSE, 
    sep="\t",
    na = "NA",
    file=paste0(i, ".txt")
  )
})


# Filter ENSG coords file
coords <- fread(ref)[V1 %in% unique(result$GENE)]
fwrite(coords,
  file = gsub(".txt", ".filtered.txt", ref),
  quote = FALSE,
  sep = "\t", 
  col.names = FALSE
)

lapply(covs, function(i){
  x <- readLines(file.path(gwas, paste0("cov_", i, ".geneset")))
  x <- x[x %in% unique(result$GENE)]
  writeLines(x, con = file.path(gwas, paste0("cov_", i, ".filtered.geneset")))
})



# # cleanup
# system("squeue -j $SLURM_JOBID --Format=TimeUsed")



















