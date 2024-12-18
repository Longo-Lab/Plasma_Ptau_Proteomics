#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 01-all_gene_lists.R
#'
#' Version: 0.0.10
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
#'   v0.0.10 - fixing the symbol lookups by creating a symbol column where possible
#'   v0.0.9 - adding TMTSoma-CSF
#'   v0.0.8 - add background gene lists for proteomic data
#'   v0.0.8 - revamp naming in all enrichment scripts/filenames to reference version and orthotype
#'   v0.0.7 - add PTX trial dataset
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
#'   - TMT-AD: TMT-LP proteomic modules from Johnson, et al. (2022)
#'   - TMTSoma-CSF: TMT-LP & Soma proteomic modules from Dammer, et al. (2024)
#'
#' ---------------------------

library(data.table)
library(biomaRt)
library(readxl)
library(tidyverse)

wk_dir <- "/labs/flongo/reference/biodomains"
ens_ver <- "113"
orth <- "one2one"


# Functions -------------------------------------------------------------

# import functions 
source("/labs/flongo/scripts/convert_orthos.R")


# TREAT-AD ---------------------------------------------------------------------

# File from Synapse (syn25428992) and mouse (syn26592124) no longer need to convert
treatad <- readRDS(file.path(wk_dir, "annotated_biodomains_hsap.RDS")) %>%
  select(Biodomain, symbol) %>% 
  unnest(symbol) %>% 
  distinct()
names(treatad) <- c("module", "symbol")
fwrite(
  convert_orthos(treatad, 
    in_col = "symbol",
    in_species = "hsapiens",
    in_type = "symbol", 
    out_col = "GENE",
    out_species = "hsapiens",
    out_type = "ensID", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("TREAT-AD", ens_ver, orth, "Human.Genes.txt", sep = ".")
  ), 
  sep = "\t"
)

treatad <- readRDS(file.path(wk_dir, "annotated_biodomains_mmus_Feb24.rds")) %>%
  select(Biodomain, symbol) %>% 
  unnest(symbol) %>% 
  distinct()
names(treatad) <- c("module", "symbol")
fwrite(
  convert_orthos(treatad, 
    in_col = "symbol",
    in_species = "mmusculus",
    in_type = "symbol", 
    out_col = "GENE",
    out_species = "mmusculus",
    out_type = "ensID", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("TREAT-AD", ens_ver, orth, "Mouse.Genes.txt", sep = ".")
  ), 
  sep = "\t"
)


# Mostafavi, et al. ------------------------------------------------------------

# File from Table S3 (https://doi.org/10.1038/s41593-018-0154-9)
mostafavi <- data.table(
  read_xlsx(
    file.path(wk_dir, "NIHMS962042-supplement-Sup_Table_3.xlsx"), 
    skip = 4
  )
)

names(mostafavi) <- c("module", "symbol")

# human genes lookup ensID
fwrite(
  convert_orthos(mostafavi, 
    in_col = "symbol",
    in_species = "hsapiens",
    in_type = "symbol", 
    out_col = "GENE",
    out_species = "hsapiens",
    out_type = "ensID", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("Mostafavi_etal", ens_ver, orth, "Human.Genes.txt", sep = ".")
  ), 
  sep = "\t"
)

# convert to mouse ensID
mostafavi <- convert_orthos(mostafavi,
  in_col = "symbol",
  in_species = "hsapiens",
  in_type = "symbol", 
  out_col = "GENE",
  out_species = "mmusculus",
  out_type = "ensID", 
  ortho_type = orth,
  ver = ens_ver
) 

# lookup mouse symbols
fwrite(
  convert_orthos(mostafavi,
    in_col = "GENE",
    in_species = "mmusculus",
    in_type = "ensID", 
    out_col = "symbol",
    out_species = "mmusculus",
    out_type = "symbol", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("Mostafavi_etal", ens_ver, orth, "Mouse.Genes.txt", sep = ".")
  ), 
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

# human genes lookup symbol
fwrite(
  convert_orthos(milind, 
    in_col = "GENE",
    in_species = "hsapiens",
    in_type = "ensID", 
    out_col = "symbol",
    out_species = "hsapiens",
    out_type = "symbol", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("Milind_etal", ens_ver, orth, "Human.Genes.txt", sep = ".")
  ), 
  sep = "\t"
)

# convert to mouse ensID
milind <- convert_orthos(milind,
  in_col = "GENE",
  in_species = "hsapiens",
  in_type = "ensID", 
  out_col = "GENE",
  out_species = "mmusculus",
  out_type = "ensID", 
  ortho_type = orth,
  ver = ens_ver
) 

# lookup mouse symbols
fwrite(
  convert_orthos(milind,
    in_col = "GENE",
    in_species = "mmusculus",
    in_type = "ensID", 
    out_col = "symbol",
    out_species = "mmusculus",
    out_type = "symbol", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("Milind_etal", ens_ver, orth, "Mouse.Genes.txt", sep = ".")
  ), 
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

# human already have symbol, write
fwrite(wan, 
  file = file.path(wk_dir, 
    paste("Wan_etal", ens_ver, orth, "Human.Genes.txt", sep = ".")
  ), 
  sep = "\t"
)

# convert to mouse ensID
wan <- convert_orthos(wan,
  in_col = "GENE",
  in_species = "hsapiens",
  in_type = "ensID", 
  out_col = "GENE",
  out_species = "mmusculus",
  out_type = "ensID", 
  ortho_type = orth,
  ver = ens_ver
) 
wan[, symbol := NULL]

# lookup mouse symbols
fwrite(
  convert_orthos(wan,
    in_col = "GENE",
    in_species = "mmusculus",
    in_type = "ensID", 
    out_col = "symbol",
    out_species = "mmusculus",
    out_type = "symbol", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("Wan_etal", ens_ver, orth, "Mouse.Genes.txt", sep = ".")
  ), 
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

# background genes
fwrite(
  list(unique(
    convert_orthos(tmt, 
      in_col = "symbol",
      in_species = "hsapiens",
      in_type = "symbol", 
      out_col = "GENE",
      out_species = "hsapiens",
      out_type = "ensID", 
      ortho_type = orth,
      ver = ens_ver
    )$GENE
  )),
  file = file.path(wk_dir, 
    paste("TMT-AD", ens_ver, orth, "Human.Background.txt", sep = ".")
  )
)
fwrite(
  list(unique(
    convert_orthos(tmt, 
      in_col = "symbol",
      in_species = "hsapiens",
      in_type = "symbol", 
      out_col = "GENE",
      out_species = "mmusculus",
      out_type = "ensID", 
      ortho_type = orth,
      ver = ens_ver
    )$GENE
  )),
  file = file.path(wk_dir, 
    paste("TMT-AD", ens_ver, orth, "Mouse.Background.txt", sep = ".")
  )
)

# filter for module members
tmt <- tmt[`net$colors` != "grey"]
tmt <- tmt[symbol != "0"]
tmt <- tmt[, c("module", "color", "kME") := tstrsplit(kMEtableSortVector, "\\|")]
tmt <- tmt[, c("symbol", "uniprotID", "module", "color", "kME")]
# A gene can be part of multiple proteins / belong to multiple modules
# A gene symbol can also have multiple ensembl gene ID... so many-to-many join
# And here: https://github.com/Longo-Lab/scripts/blob/2e9fb01c9a7d4aef18542bbebbbdc1c55402a8c5/biodomain_TMT-AD.R#L117

# # if you want to drop the duplicates here:
# tmt <- unique(tmt, by=c("symbol", "module"))

# human genes lookup ensID
fwrite(
  convert_orthos(tmt, 
    in_col = "symbol",
    in_species = "hsapiens",
    in_type = "symbol", 
    out_col = "GENE",
    out_species = "hsapiens",
    out_type = "ensID", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("TMT-AD", ens_ver, orth, "Human.Genes.txt", sep = ".")
  ),
  sep = "\t"
)

# convert to mouse ensID
tmt <- convert_orthos(tmt,
  in_col = "symbol",
  in_species = "hsapiens",
  in_type = "symbol", 
  out_col = "GENE",
  out_species = "mmusculus",
  out_type = "ensID", 
  ortho_type = orth,
  ver = ens_ver
)

# lookup mouse symbols
fwrite(
  convert_orthos(tmt,
    in_col = "GENE",
    in_species = "mmusculus",
    in_type = "ensID", 
    out_col = "symbol",
    out_species = "mmusculus",
    out_type = "symbol", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("TMT-AD", ens_ver, orth, "Mouse.Genes.txt", sep = ".")
  ),
  sep = "\t"
)


# TMT-Resilience -----------------------------------------------------------------
# File from Table S3/4 (DOI:https://doi.org/10.1016/j.mcpro.2023.100542)
# tranferred commands from TMT-Resilience/save_TMT_data.rbindlist

# Load module names
m_names <- fread(file.path(wk_dir, 'TMT-Resilience/module_names.csv')) %>% 
  pull(name, module)

# BA6
ba6 <- data.table(read_excel(
  file.path(wk_dir, 'TMT-Resilience/mmc11.xlsx'), 
  sheet = '3.BA6 module assignments', 
  skip = 2
))

# Check format is consistent
assertthat::are_equal(str_count(ba6$ProteinIDs, '\\|') %>% unique(), 1)
assertthat::are_equal(str_count(ba6$kMEtableSortVector1, '\\|') %>% unique(), 1)

ba6 <- ba6[, c("symbol", "uniprotID") := tstrsplit(ProteinIDs, "\\|")]

# background genes
fwrite(
  list(unique(
    convert_orthos(ba6, 
      in_col = "symbol",
      in_species = "hsapiens",
      in_type = "symbol", 
      out_col = "GENE",
      out_species = "hsapiens",
      out_type = "ensID", 
      ortho_type = orth,
      ver = ens_ver
    )$GENE
  )),
  file = file.path(wk_dir, 
    paste("TMT-Resilience", ens_ver, orth, "Human.Background.txt", sep = ".")
  )
)
fwrite(
  list(unique(
    convert_orthos(ba6, 
      in_col = "symbol",
      in_species = "hsapiens",
      in_type = "symbol", 
      out_col = "GENE",
      out_species = "mmusculus",
      out_type = "ensID", 
      ortho_type = orth,
      ver = ens_ver
    )$GENE
  )),
  file = file.path(wk_dir, 
    paste("TMT-Resilience", ens_ver, orth, "Mouse.Background.txt", sep = ".")
  )
)

# filter grey
ba6 <- ba6[`net$colors` != "grey"]

# split and add module names
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

# human genes lookup ensID
fwrite(
  convert_orthos(ba6, 
    in_col = "symbol",
    in_species = "hsapiens",
    in_type = "symbol", 
    out_col = "GENE",
    out_species = "hsapiens",
    out_type = "ensID", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("TMT-Resilience", ens_ver, orth, "Human.Genes.txt", sep = ".")
  ),
  sep = "\t"
)

# convert to mouse ensID
ba6 <- convert_orthos(ba6,
  in_col = "symbol",
  in_species = "hsapiens",
  in_type = "symbol", 
  out_col = "GENE",
  out_species = "mmusculus",
  out_type = "ensID", 
  ortho_type = orth,
  ver = ens_ver
) 

# lookup mouse symbols
fwrite(
  convert_orthos(ba6,
    in_col = "GENE",
    in_species = "mmusculus",
    in_type = "ensID", 
    out_col = "symbol",
    out_species = "mmusculus",
    out_type = "symbol", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("TMT-Resilience", ens_ver, orth, "Mouse.Genes.txt", sep = ".")
  ),
  sep = "\t"
)

# TMT-PTXtrial ----------------------------------------------------------------

# Capture background 
tmt <- fread("PTX_ModuleAssignments_ALL_ME.csv")
tmt <- tmt[, c("symbol", "uniprotID") := tstrsplit(UniqueID, "\\|")]
tmt <- tmt[, c("module", "color", "kME") := tstrsplit(kMEtableSortVector, " |\\|")]

# split complexes
tmt <- tmt[, 
  .(symbol = unlist(tstrsplit(symbol, "; "))), 
  by = c(names(tmt)[names(tmt) != "symbol"])
]
# drop IG@ suffixes after "-"
tmt[, symbol := gsub("(IG\\w+)-.*", "\\1", symbol)]
# remove duplicates per module
tmt <- unique(tmt, by = c("symbol", "module"))


# background genes
fwrite(
  list(unique(
    convert_orthos(tmt, 
      in_col = "symbol",
      in_species = "hsapiens",
      in_type = "symbol", 
      out_col = "GENE",
      out_species = "hsapiens",
      out_type = "ensID", 
      ortho_type = orth,
      ver = ens_ver
    )$GENE
  )),
  file = file.path(wk_dir, 
    paste("TMT-PTXtrial", ens_ver, orth, "Human.Background.txt", sep = ".")
  )
)
fwrite(
  list(unique(
    convert_orthos(tmt, 
      in_col = "symbol",
      in_species = "hsapiens",
      in_type = "symbol", 
      out_col = "GENE",
      out_species = "mmusculus",
      out_type = "ensID", 
      ortho_type = orth,
      ver = ens_ver
    )$GENE
  )),
  file = file.path(wk_dir, 
    paste("TMT-PTXtrial", ens_ver, orth, "Mouse.Background.txt", sep = ".")
  )
)

# Filter outside modules
tmt <- tmt[moduleColors != "grey"]
tmt <- tmt[symbol != "0"]
tmt <- tmt[, c("symbol", "uniprotID", "module", "color", "kME")]
# A gene can be part of multiple proteins / belong to multiple modules
# A gene symbol can also have multiple ensembl gene ID... so many-to-many join
# And here: https://github.com/Longo-Lab/scripts/blob/2e9fb01c9a7d4aef18542bbebbbdc1c55402a8c5/biodomain_TMT-AD.R#L117

# # if you want to drop the duplicates here:
# tmt <- unique(tmt, by=c("symbol", "module"))

# human genes lookup ensID
fwrite(
  convert_orthos(tmt, 
    in_col = "symbol",
    in_species = "hsapiens",
    in_type = "symbol", 
    out_col = "GENE",
    out_species = "hsapiens",
    out_type = "ensID", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("TMT-PTXtrial", ens_ver, orth, "Human.Genes.txt", sep = ".")
  ),
  sep = "\t"
)

# convert to mouse ensID
tmt <- convert_orthos(tmt,
  in_col = "symbol",
  in_species = "hsapiens",
  in_type = "symbol", 
  out_col = "GENE",
  out_species = "mmusculus",
  out_type = "ensID", 
  ortho_type = orth,
  ver = ens_ver
) 

# lookup mouse symbols
fwrite(
  convert_orthos(tmt,
    in_col = "GENE",
    in_species = "mmusculus",
    in_type = "ensID", 
    out_col = "symbol",
    out_species = "mmusculus",
    out_type = "symbol", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("TMT-PTXtrial", ens_ver, orth, "Mouse.Genes.txt", sep = ".")
  ),
  sep = "\t"
)

# TMTSoma-CSF ----------------------------------------------------------------
# File from Table S5/7 (DOI:https://doi.org/10.1126/scitranslmed.adn3504)
# tranferred commands from TMTSoma-CSF/*

# Load module names
m_names <- fread(file.path(wk_dir, 'TMTSoma-CSF/TMTSoma-CSF.names.csv')) 
m_names[, name := paste(module, name)]

# Capture background 
tmt <- fread("TMTSoma-CSF/adn3504_Data_file_S1_T7.csv")
tmt <- tmt[, c("symbol", "uniprotID", "platform", "source") := tstrsplit(UniqueID, "\\|")]
tmt <- tmt[, c("module", "color", "kME") := tstrsplit(`Module|kME`, "(\\|| )")]

# split complexes
tmt <- tmt[, 
  .(symbol = unlist(tstrsplit(symbol, "; "))), 
  by = c(names(tmt)[names(tmt) != "symbol"])
]
# remove duplicates per module
tmt <- unique(tmt, by = c("symbol", "module"))

# background genes
fwrite(
  list(unique(
    convert_orthos(tmt, 
      in_col = "symbol",
      in_species = "hsapiens",
      in_type = "symbol", 
      out_col = "GENE",
      out_species = "hsapiens",
      out_type = "ensID", 
      ortho_type = orth,
      ver = ens_ver
    )$GENE
  )),
  file = file.path(wk_dir, 
    paste("TMTSoma-CSF", ens_ver, orth, "Human.Background.txt", sep = ".")
  )
)
fwrite(
  list(unique(
    convert_orthos(tmt, 
      in_col = "symbol",
      in_species = "hsapiens",
      in_type = "symbol", 
      out_col = "GENE",
      out_species = "mmusculus",
      out_type = "ensID", 
      ortho_type = orth,
      ver = ens_ver
    )$GENE
  )),
  file = file.path(wk_dir, 
    paste("TMTSoma-CSF", ens_ver, orth, "Mouse.Background.txt", sep = ".")
  )
)

# filter grey
tmt <- tmt[`Member of Module` != "grey"]

# split and add module names
tmt <- tmt[symbol != "0"]
tmt$module <- plyr::mapvalues(
  x = tmt$module, 
  from = m_names$module, 
  to = m_names$name
)

tmt <- tmt[, c("symbol", "uniprotID", "module", "color", "kME")]
# A gene can be part of multiple proteins / belong to multiple modules
# A gene symbol can also have multiple ensembl gene ID... so many-to-many join
# And here: https://github.com/Longo-Lab/scripts/blob/2e9fb01c9a7d4aef18542bbebbbdc1c55402a8c5/biodomain_TMT-AD.R#L117

# # if you want to drop the duplicates here:
# tmt <- unique(tmt, by=c("symbol", "module"))

# human genes lookup ensID
fwrite(
  convert_orthos(tmt, 
    in_col = "symbol",
    in_species = "hsapiens",
    in_type = "symbol", 
    out_col = "GENE",
    out_species = "hsapiens",
    out_type = "ensID", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("TMTSoma-CSF", ens_ver, orth, "Human.Genes.txt", sep = ".")
  ),
  sep = "\t"
)

# convert to mouse ensID
tmt <- convert_orthos(tmt,
  in_col = "symbol",
  in_species = "hsapiens",
  in_type = "symbol", 
  out_col = "GENE",
  out_species = "mmusculus",
  out_type = "ensID", 
  ortho_type = orth,
  ver = ens_ver
) 

# lookup mouse symbols
fwrite(
  convert_orthos(tmt,
    in_col = "GENE",
    in_species = "mmusculus",
    in_type = "ensID", 
    out_col = "symbol",
    out_species = "mmusculus",
    out_type = "symbol", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("TMTSoma-CSF", ens_ver, orth, "Mouse.Genes.txt", sep = ".")
  ),
  sep = "\t"
)

# Soma-plasma ----------------------------------------------------------------
# File from Table S7 (DOI:https://doi.org/10.1126/science.aaq1327)
# tranferred to Soma-plasma/*

# Load module names
m_names <- fread(file.path(wk_dir, 'Soma-plasma/Soma-plasma.names.csv')) 
m_names[, name := paste(module, name)]

# Capture background 
tmt <- fread("Soma-plasma/aaq1327_excel_tables_T7.csv")
setnames(tmt, 
  old = c("Module Annotation", "Module Colour", "Gene Symbol"),
  new = c("module", "color", "symbol")
)

# split complexes harder with this one... some are " " or "  " one is ","
tmt <- tmt[, 
  .(symbol = unlist(tstrsplit(symbol, " |  |,"))), 
  by = c(names(tmt)[names(tmt) != "symbol"])
]
# remove duplicates per module
tmt <- unique(tmt, by = c("symbol", "module"))

# background genes
fwrite(
  list(unique(
    convert_orthos(tmt, 
      in_col = "symbol",
      in_species = "hsapiens",
      in_type = "symbol", 
      out_col = "GENE",
      out_species = "hsapiens",
      out_type = "ensID", 
      ortho_type = orth,
      ver = ens_ver
    )$GENE
  )),
  file = file.path(wk_dir, 
    paste("Soma-plasma", ens_ver, orth, "Human.Background.txt", sep = ".")
  )
)
fwrite(
  list(unique(
    convert_orthos(tmt, 
      in_col = "symbol",
      in_species = "hsapiens",
      in_type = "symbol", 
      out_col = "GENE",
      out_species = "mmusculus",
      out_type = "ensID", 
      ortho_type = orth,
      ver = ens_ver
    )$GENE
  )),
  file = file.path(wk_dir, 
    paste("Soma-plasma", ens_ver, orth, "Mouse.Background.txt", sep = ".")
  )
)

# filter grey
tmt <- tmt[color != "gray"]

# split and add module names
tmt[, module := gsub("PM", "M", module)]
tmt$module <- plyr::mapvalues(
  x = tmt$module, 
  from = m_names$module, 
  to = m_names$name
)

tmt <- tmt[, c("symbol", "UniProt", "module", "color", "kWithin")]
# A gene can be part of multiple proteins / belong to multiple modules
# A gene symbol can also have multiple ensembl gene ID... so many-to-many join
# And here: https://github.com/Longo-Lab/scripts/blob/2e9fb01c9a7d4aef18542bbebbbdc1c55402a8c5/biodomain_TMT-AD.R#L117

# # if you want to drop the duplicates here:
# tmt <- unique(tmt, by=c("symbol", "module"))

# human genes lookup ensID
fwrite(
  convert_orthos(tmt, 
    in_col = "symbol",
    in_species = "hsapiens",
    in_type = "symbol", 
    out_col = "GENE",
    out_species = "hsapiens",
    out_type = "ensID", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("Soma-plasma", ens_ver, orth, "Human.Genes.txt", sep = ".")
  ),
  sep = "\t"
)

# convert to mouse ensID
tmt <- convert_orthos(tmt,
  in_col = "symbol",
  in_species = "hsapiens",
  in_type = "symbol", 
  out_col = "GENE",
  out_species = "mmusculus",
  out_type = "ensID", 
  ortho_type = orth,
  ver = ens_ver
) 

# lookup mouse symbols
fwrite(
  convert_orthos(tmt,
    in_col = "GENE",
    in_species = "mmusculus",
    in_type = "ensID", 
    out_col = "symbol",
    out_species = "mmusculus",
    out_type = "symbol", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("Soma-plasma", ens_ver, orth, "Mouse.Genes.txt", sep = ".")
  ),
  sep = "\t"
)




# Expression correlation stat prep -------------------------------------

# TREAT-AD --------------------------
# table downloaded from synapse
# synapse get -q "SELECT * FROM syn22758536"
expr = fread("SYNAPSE_TABLE_QUERY_96970654.csv")

expr[, RNA_stat := -log10(RNA_fdr_CorPVal) * RNA_TE]
expr[, Pro_stat := -log10(Pro_fdr_CorPVal) * Pro_TE]
setnames(expr, c("ENSG", "GName"), c("GENE", "symbol"))
expr <- expr[, .(GENE, symbol, RNA_TE, RNA_fdr_CorPVal, RNA_stat, Pro_TE, Pro_fdr_CorPVal, Pro_stat)]

# human genes already have ids
fwrite(expr,
  file = file.path(wk_dir, 
    paste("TREAT-AD", ens_ver, orth, "Human.Corr.txt", sep = ".")
  ),
  sep = "\t"
)

# convert to mouse ensID
expr <- convert_orthos(expr,
  in_col = "GENE",
  in_species = "hsapiens",
  in_type = "ensID", 
  out_col = "GENE",
  out_species = "mmusculus",
  out_type = "ensID", 
  ortho_type = orth,
  ver = ens_ver
) 
expr[, symbol := NULL]

# lookup mouse symbols
fwrite(
  convert_orthos(expr, 
    in_col = "GENE",
    in_species = "mmusculus",
    in_type = "ensID", 
    out_col = "symbol",
    out_species = "mmusculus",
    out_type = "symbol", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("TREAT-AD", ens_ver, orth, "Mouse.Corr.txt", sep = ".")
  ),
  sep = "\t"
)


# TMTSoma-CSF ---------------------------------------------
expr <- fread("TMTSoma-CSF/adn3504_Data_file_S1_T5.csv")

# prep and filter
expr[, Pro_stat := -log10(`FDR (BH)`) * `log2FC\nAD-Control`]
expr <- expr[order(-`F-Value`)][!is.na(Pro_stat)][, head(.SD, 1), by = Symbol]
setnames(expr, 
  old = c("Symbol", "log2FC\nAD-Control", "FDR (BH)"), 
  new = c("symbol", "log2FC", "FDR")
)
expr <- expr[, .(symbol, log2FC, `F-Value`, FDR, Pro_stat)]

# human genes lookup ensID
fwrite(
  convert_orthos(expr, 
    in_col = "symbol",
    in_species = "hsapiens",
    in_type = "symbol", 
    out_col = "GENE",
    out_species = "hsapiens",
    out_type = "ensID", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("TMTSoma-CSF", ens_ver, orth, "Human.Corr.txt", sep = ".")
  ),
  sep = "\t"
)

# convert to mouse ensID
expr <- convert_orthos(expr,
  in_col = "symbol",
  in_species = "hsapiens",
  in_type = "symbol", 
  out_col = "GENE",
  out_species = "mmusculus",
  out_type = "ensID", 
  ortho_type = orth,
  ver = ens_ver
) 

# lookup mouse symbols
fwrite(
  convert_orthos(expr, 
    in_col = "GENE",
    in_species = "mmusculus",
    in_type = "ensID", 
    out_col = "symbol",
    out_species = "mmusculus",
    out_type = "symbol", 
    ortho_type = orth,
    ver = ens_ver,
    keep_original = TRUE
  ), 
  file = file.path(wk_dir, 
    paste("TMTSoma-CSF", ens_ver, orth, "Mouse.Corr.txt", sep = ".")
  ),
  sep = "\t"
)

# Soma-plasma ---------------------------------------------
# No AD status, nothing to corr



