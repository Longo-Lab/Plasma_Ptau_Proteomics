#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 01-all_gene_lists.R
#'
#' Version: 0.0.6
#'
#' Purpose of script: Save modules list from multiple reference sources
#'
#' Author: Robert R Butler III, Crystal Han
#'
#' Date Created: 2024-01-22
#'
#' Copyright (c) 2024
#' Email: rrbutler@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'   Future - revamp naming in all enrichment scripts/filenames to reference version and orthotype
#'   v0.0.6 - exporting functions to scripts for utility, now sources
#'   v0.0.6 - generating one2one orthologs by default
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
#'   - TreatAD: TREAT-AD pathway traced modules from Cary, et al. (2022)
#'   - Mostafavi: SpeakEasy gene modules from Mostafavi, et al. (2018)
#'   - Milind: Submodules of human/mouse co-expression from Milind, et al. (2020)
#'   - Wan: Submodules of human/mouse co-expression from Wan, et al. (2020)
#'   - TmtAD: TMT-LP proteomic modules from Johnson, et al. (2022)
#'
#' ---------------------------

library(data.table)
library(biomaRt)
library(readxl)
library(tidyverse)

wk_dir <- "/labs/flongo/reference/biodomains"
ensembl_ver <- "111"


# Functions -------------------------------------------------------------

# import functions 
source("/labs/flongo/scripts/convert_orthos.R")


# TREAT-AD ---------------------------------------------------------------------

# File from Synapse (syn25428992) and mouse (syn26592124) no longer need to convert
treatad <- readRDS(file.path(wk_dir, "annotated_biodomains_hsap.RDS")) %>%
  select(Biodomain, ensembl_id) %>% 
  unnest(ensembl_id) %>% 
  distinct()
names(treatad) <- c("module", "GENE")
fwrite(treatad, file = file.path(wk_dir, "TREAT-AD.Human.Genes.txt"), sep = "\t")

treatad <- readRDS(file.path(wk_dir, "annotated_biodomains_mmus_Feb24.rds")) %>%
  select(Biodomain, ensembl_id) %>% 
  unnest(ensembl_id) %>% 
  distinct()
names(treatad) <- c("module", "GENE")
fwrite(treatad, file = file.path(wk_dir, "TREAT-AD.Mouse.Genes.txt"), sep = "\t")


# Mostafavi, et al. ------------------------------------------------------------

# File from Table S3 (https://doi.org/10.1038/s41593-018-0154-9)
mostafavi <- data.table(
  read_xlsx(
    file.path(wk_dir, "NIHMS962042-supplement-Sup_Table_3.xlsx"), 
    skip = 4
  )
)

names(mostafavi) <- c("module", "symbol")

fwrite(
  convert_orthos(mostafavi, 
    in_col = "symbol",
    in_species = "hsapiens",
    in_type = "symbol", 
    out_col = "GENE",
    out_species = "hsapiens",
    out_type = "ensID", 
    ver = ensembl_ver
  ), 
  file = file.path(wk_dir, "Mostafavi_etal.Human.Genes.txt"), 
  sep = "\t"
)

fwrite(
  convert_orthos(mostafavi,
    in_col = "symbol",
    in_species = "hsapiens",
    in_type = "symbol", 
    out_col = "GENE",
    out_species = "mmusculus",
    out_type = "ensID", 
    ortho_type = "one2one",
    ver = ensembl_ver
  ), 
  file = file.path(wk_dir, "Mostafavi_etal.Mouse.Genes.txt"), 
  sep = "\t"
)


# Milind, et al. ---------------------------------------------------------------

# Files from Synapse (rosmap genes: syn23660902, mayo genes: syn23660898, mssm genes: syn23660900)
study <- c("rosmap", "mayo", "mssm")
names(study) <- c("ROSMAP", "Mayo", "MSSM")

# read in and join the tables
milind <- lapply(study, function(i) { 
  readRDS(file.path(wk_dir, paste0(i, "_submodule_sig_genes.RDS"))) 
})
milind <- rbindlist(milind, idcol = "study")

setnames(milind, c("Submodule", "Gene"), c("module", "GENE"))

milind_groups <- fread(file.path(wk_dir, "Milind_etal.group_annots.txt"))
milind <- milind[module %in% milind_groups$module]

fwrite(milind, file = file.path(wk_dir, "Milind_etal.Human.Genes.txt"), sep = "\t")
fwrite(
  convert_orthos(milind,
    in_col = "GENE",
    in_species = "hsapiens",
    in_type = "ensID", 
    out_col = "GENE",
    out_species = "mmusculus",
    out_type = "ensID", 
    ortho_type = "one2one",
    ver = ensembl_ver
  ), 
  file = file.path(wk_dir, "Milind_etal.Mouse.Genes.txt"), 
  sep = "\t"
)


# Wan, et al. ------------------------------------------------------------------

# Originally AD_consensus_clusers.tsv, found in /labs/flongo/ScAtt/combined_module_enrichment
# But current directory only contains AD_modules_full.tsv, which is same as AD_consensus_clusers.tsv but w/ 2 additional columns
ad_modules <- fread(file.path(wk_dir, "AD_modules_full.tsv")) %>% 
  select(-Consensus_Cluster, -Mmus_GeneID)

# From Figure S1 legend in Wan et al (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7428328/)
string <- list(
  "A" = c("TCXblue", "PHGyellow", "IFGyellow"),
  "B" = c("FPturquoise", "TCXturquoise", "STGblue", "PHGturquoise", "IFGturquoise", "DLPFCblue", "CBEturquoise"),
  "C" = c("CBEyellow", "PHGbrown", "TCXgreen", "DLPFCyellow", "STGbrown", "FPyellow", "IFGbrown"),
  "D" = c("CBEbrown", "DLPFCbrown", "PHGgreen", "STGyellow", "FPblue", "TCXyellow", "IFGblue"),
  "E" = c("FPbrown", "TCXbrown", "STGturquoise", "DLPFCturquoise", "PHGblue", "CBEblue")
)

# Map consensus cluster to data table
ccs <- reshape2::melt(string) 
setnames(ccs, c("Module", "Consensus_Cluster"))

wan <- merge(ad_modules, ccs, sort = FALSE)
wan[, c("method", "ModuleName", "ModuleNameFull") := NULL]

setnames(wan, c("Module", "GeneID", "external_gene_name", "brainRegion"), 
         c("module", "GENE", "symbol", "Brain_Region"))

fwrite(wan, file = file.path(wk_dir, "Wan_etal.Human.Genes.txt"), sep = "\t")
fwrite(
  convert_orthos(wan,
    in_col = "GENE",
    in_species = "hsapiens",
    in_type = "ensID", 
    out_col = "GENE",
    out_species = "mmusculus",
    out_type = "ensID", 
    ortho_type = "one2one",
    ver = ensembl_ver
  ), 
  file = file.path(wk_dir, "Wan_etal.Mouse.Genes.txt"), 
  sep = "\t"
)


# TMT-AD -----------------------------------------------------------------------

# File from Table S4 (https://doi.org/10.1038%2Fs41593-021-00999-y)
tmt <- data.table(read_xlsx(
  file.path(wk_dir, "41593_2021_999_MOESM3_ESM.xlsx"), 
  sheet = "4.Consensus.kMEs", 
  skip = 2
))

tmt <- tmt[, c("symbol", "uniprotID") := tstrsplit(UniqueID, "\\|")]
tmt <- tmt[`net$colors` != "grey"]
tmt <- tmt[, c("module", "color", "kME") := tstrsplit(kMEtableSortVector, "\\|")]
tmt <- tmt[, c("symbol", "uniprotID", "module", "color", "kME")]
# A gene can be part of multiple proteins / belong to multiple modules
# A gene symbol can also have multiple ensembl gene ID... so many-to-many join
# And here: https://github.com/Longo-Lab/scripts/blob/2e9fb01c9a7d4aef18542bbebbbdc1c55402a8c5/biodomain_TMT-AD.R#L117

# # if you want to drop the duplicates here:
# tmt <- unique(tmt, by=c("symbol", "module"))


fwrite(
  convert_orthos(tmt, 
    in_col = "symbol",
    in_species = "hsapiens",
    in_type = "symbol", 
    out_col = "GENE",
    out_species = "hsapiens",
    out_type = "ensID", 
    ver = ensembl_ver
  ), 
  file = file.path(wk_dir, "TMT-AD.Human.Genes.txt"), 
  sep = "\t"
)
fwrite(
  convert_orthos(tmt,
    in_col = "symbol",
    in_species = "hsapiens",
    in_type = "symbol", 
    out_col = "GENE",
    out_species = "mmusculus",
    out_type = "ensID", 
    ortho_type = "one2one",
    ver = ensembl_ver
  ), 
  file = file.path(wk_dir, "TMT-AD.Mouse.Genes.txt"), 
  sep = "\t"
)


# TMT-Resilience -----------------------------------------------------------------
# File from Table S3/4 (DOI:https://doi.org/10.1016/j.mcpro.2023.100542)
# tranferred commands from TMT/save_TMT_data.rbindlist

# Load module names
m_names <- fread(file.path(wk_dir, 'TMT/module_names.csv')) %>% 
  pull(name, module)

# BA6
ba6 <- data.table(read_excel(
  file.path(wk_dir, 'TMT/mmc11.xlsx'), 
  sheet = '3.BA6 module assignments', 
  skip = 2
))

# Check format is consistent
assertthat::are_equal(str_count(ba6$ProteinIDs, '\\|') %>% unique(), 1)
assertthat::are_equal(str_count(ba6$kMEtableSortVector1, '\\|') %>% unique(), 1)

# filter grey
ba6 <- ba6[`net$colors` != "grey"]

# split and add module names
ba6 <- ba6[, c("symbol", "uniprotID") := tstrsplit(ProteinIDs, "\\|")]
ba6 <- ba6[symbol != "0"]
ba6 <- ba6[, c("module", "color", "kME") := tstrsplit(kMEtableSortVector1, "(\\|| )")]
ba6[, module := paste(module, m_names[[module]]), by=.I]

ba6 <- ba6[, c("symbol", "uniprotID", "module", "color", "kME")]
# All 7787 unique proteins, but 7345 unique gene symbols (i.e., 442 dupe genes)
# ba6 %>% group_by(Symbol) %>% summarise(n = n()) %>% arrange(-n) %>% View()
# A gene can be part of multiple proteins / belong to multiple modules
# A gene symbol can also have multiple ensembl gene ID... so many-to-many join
# And here: https://github.com/Longo-Lab/scripts/blob/2e9fb01c9a7d4aef18542bbebbbdc1c55402a8c5/biodomain_TMT-Resilience.R#L138

# # if you want to drop the duplicates here:
# ba6 <- unique(ba6, by=c("symbol", "module"))

# convert and write to file
fwrite(
  convert_orthos(ba6, 
    in_col = "symbol",
    in_species = "hsapiens",
    in_type = "symbol", 
    out_col = "GENE",
    out_species = "hsapiens",
    out_type = "ensID", 
    ver = ensembl_ver
  ), 
  file = file.path(wk_dir, "TMT-Resilience.Human.Genes.txt"), 
  sep = "\t"
)
fwrite(
  convert_orthos(ba6,
    in_col = "symbol",
    in_species = "hsapiens",
    in_type = "symbol", 
    out_col = "GENE",
    out_species = "mmusculus",
    out_type = "ensID", 
    ortho_type = "one2one",
    ver = ensembl_ver
  ), 
  file = file.path(wk_dir, "TMT-Resilience.Mouse.Genes.txt"), 
  sep = "\t"
)








