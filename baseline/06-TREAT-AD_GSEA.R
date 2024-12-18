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
## Author: Grey Cary, Robert B Butler III, Michael P Brown
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
library(optparse)
# library(biomaRt)

# Nameset helper function
nameset_default <- function(input){
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

# Login to synapse
synLogin()

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
# DE file filtering and ordering (dataProt_SS-205063.pTau.baseline.csv)
# `bm`~ biomarker
bm = fread(input)
pfilt = bm[!is.na(GENE) , .(beta, GENE)][order(-beta)]
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
bdt <- bd.tally(enr@result %>% filter(p.adjust < 0.05) %>% pull(Description))     # Removed the last field `, biofilt)` to test functionality, `biofilt` is commented atm

# annotate results with biodomain information for plotting
enr@result <- enr@result %>% 
  left_join(., biodom %>% select(ID = GO_ID, Biodomain), by = 'ID') %>%
  full_join(., biodom %>% select(Biodomain, abbr, label, color) %>% distinct() , by = 'Biodomain' ) %>% 
  mutate(term = Description) %>% 
  mutate(label = fct_relevel(
    label,  
    bdt$label %>% str_subset(.,'none', negate = T) %>% c('none',.) )) %>% 
  arrange(label)

saveRDS(enr, file=file.path(output, paste(nameset, "GSEA_Biodomain.rds", sep=".")))

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
htmlwidgets::saveWidget(p, file.path(output, paste(nameset, "GSEA_Biodomain.html", sep=".")),
                        selfcontained=F)





