#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 10-TREAT-AD_GSEA.R
##
## Version: 0.0.1
##
## Purpose of script: GSEA based Enrichment of TREAT-AD biodomains
## in Ptau geneset
##
## Author: Grey Cary, Robert B Butler III
##
## Date Created: 2023-03-01
##
## Copyright (c) 2023
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##   Currently hardcodes the input files (good target for options funtion)
##   looks at bioRxiv paper: https://www.biorxiv.org/content/10.1101/2021.04.05.438450v2.full
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
# library(ggplot2)
# library(ggpubr)
library(data.table)
library(synapser)
library(clusterProfiler)
library(tidyverse)
# library(biomaRt)

# directory locations
rootdir = '/mnt/box/plasma_ptau_proteomics'
wk_dir = paste0(rootdir, '/baseline')
setwd(wk_dir)
refdir = paste0(rootdir, '/biodomains')
prefix = "dataProt_SS-205063.pTau.baseline"
nameset = "TREAT-AD"

# Login to synapse
synLogin(authToken="eyJ0eXAiOiJKV1QiLCJraWQiOiJXN05OOldMSlQ6SjVSSzpMN1RMOlQ3TDc6M1ZYNjpKRU9VOjY0NFI6VTNJWDo1S1oyOjdaQ0s6RlBUSCIsImFsZyI6IlJTMjU2In0.eyJhY2Nlc3MiOnsic2NvcGUiOlsidmlldyIsImRvd25sb2FkIl0sIm9pZGNfY2xhaW1zIjp7fX0sInRva2VuX3R5cGUiOiJQRVJTT05BTF9BQ0NFU1NfVE9LRU4iLCJpc3MiOiJodHRwczovL3JlcG8tcHJvZC5wcm9kLnNhZ2ViYXNlLm9yZy9hdXRoL3YxIiwiYXVkIjoiMCIsIm5iZiI6MTcyMDU1Mjc2OCwiaWF0IjoxNzIwNTUyNzY4LCJqdGkiOiI5ODcwIiwic3ViIjoiMzUwODUyNyJ9.pocXji9g6HdmEVGwmmocFx3EnTwrW79Dw7UTjmfBm2hkzUlPoa-BRanhnT_XZ5v4zN4a-uiEn_l91vzvyl-jPXE8jU9b9En-eD1wgL2zhsnBDNTiV0DiGf5OgUzjkk3ifxoc3Qu9FS3uLFZUA0l2ia4tbLRi4q5e6taMlwsPFbSFJuWHsM8EDb6B9_jdcrSkOof_DPELq4_W0FYTLQoeW6QjzvX5eGT_nhVREqz1vvJVJDO_BtUJ3a8YVJjI-yAP_wZoXE-DFrVsssuabDogbkfgsv7T0speO3ZZBUkDVVYesLkYAdJWMDL84ZJ_kCAb8HXXy88YDqCVvGDATy8rqA")

# Functions ---------------------------
# function to count enriched terms by biodomain
bd.tally <- function( enrTerm, biodomDefTbl = biodom ){
  bdt <- bind_cols(
    biodomDefTbl %>% select(Biodomain, abbr, label, color) %>% distinct(),
    n_term = map_dbl( 
      unique(biodomDefTbl$Biodomain),
      ~ biodomDefTbl %>% filter(Biodomain == .x) %>% pull(GOterm_Name) %>% length()),
    n_sig_term = map_dbl( 
      unique(biodomDefTbl$Biodomain),
      ~ biodomDefTbl %>% filter( GOterm_Name %in% enrTerm, Biodomain == .x) %>% pull(GO_ID) %>% unique() %>% length())) %>% 
    mutate(n_sig_term = ifelse(
      label == 'none', 
      length( setdiff( enrTerm, biodom$GOterm_Name) ), 
      n_sig_term )) %>% 
    mutate(Biodomain = fct_reorder(Biodomain, n_sig_term, .desc = F)) %>% 
    arrange(Biodomain)
  
  bdt$prop <- bdt$n_sig_term / bdt$n_term
  
  return(bdt)
}

# Main ------------------------------------------
# DE file filtering and ordering
ptau = fread(paste(prefix, "csv", sep="."))
pfilt = ptau[!is.na(GENE) , .(beta, GENE)][order(-beta)]
pfilt = pull(pfilt, beta, GENE)

# load biodom reference table
biodom <- full_join(
  # biodomains
  readRDS(synGet('syn25428992')$path),
  # domain labels
  read_csv(synGet('syn26856828')$path, show_col_types = F),
  by = c('Biodomain'='domain')
) %>%
  mutate(Biodomain = case_when(Biodomain == 'none' ~ NA_character_, T ~ Biodomain))

# # Potentially filter biodom by only GOterms with proteins in Somalogic chip ----
# biofilt = data.table(biodom)
# biofilt[, isInMySet := lapply(ensembl_id, function(x) x %in% ptau$GENE)]
# biofilt[lapply(biofilt$isInMySet, sum) > 0]

# perform enrichment with clusterProfiler
enr <- gseGO(geneList = pfilt, ont = 'all', keyType = "ENSEMBL",
             OrgDb = org.Hs.eg.db::org.Hs.eg.db, # probably need to specify the mouse DB here: org.Mm.eg.db::org.Mm.eg.db
             eps = 0, pvalueCutoff = 1, scoreType = 'std')

# count the number of significant terms per biodomain
bdt <- bd.tally(enr@result %>% filter(p.adjust < 0.05) %>% pull(Description), biofilt)

# annotate results with biodomain information for plotting
enr@result <- enr@result %>% 
  left_join(., biodom %>% select(ID = GO_ID, Biodomain), by = 'ID') %>%
  full_join(., biodom %>% select(Biodomain, abbr, label, color) %>% distinct() , by = 'Biodomain' ) %>% 
  mutate(term = Description) %>% 
  mutate(label = fct_relevel(
    label,  
    bdt$label %>% str_subset(.,'none', negate = T) %>% c('none',.) )) %>% 
  arrange(label)

saveRDS(enr, file=paste(prefix, nameset, "GSEA_Biodomain.rds", sep="."))

# plot
plot <- ggplot(enr@result, aes( label, NES )) + theme_minimal()+
  scale_color_identity(guide = 'none')+ scale_fill_identity(guide = 'none')+
  scale_x_discrete(drop = F)+
  scale_size(
    limits = -log10(enr@result$p.adjust) %>% range(na.rm = T),
    range = c(0.05,5))+
  geom_jitter( data = subset(enr@result, p.adjust > 0.05 ),
               color = 'grey85',  size = .3 )+
  geom_violin(data = subset(enr@result, p.adjust <= 0.05 & NES < 0 ),
              scale='width', aes(col = color), lwd = .5 )+
  geom_violin(data = subset(enr@result, p.adjust <= 0.05 & NES > 0 ),
              scale='width', aes(col = color), lwd = .5 )+
  geom_jitter( data = subset(enr@result, p.adjust <= 0.05 ),
               aes(fill = color, size = -log10(p.adjust), text = term),
               color = 'grey20', alpha = .4, shape = 21 ) +
  labs(x='') + coord_flip() +
  theme(legend.position = 'none')

# use plotly to interactively explore plots
p = plotly::ggplotly(plot, tooltip = 'text')
htmlwidgets::saveWidget(p, paste(prefix, nameset, "GSEA_Biodomain.html", sep="."),
                        selfcontained=F)





