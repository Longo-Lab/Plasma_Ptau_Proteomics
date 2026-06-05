
# setup -------------------------------------------------------------------

syn.client <- reticulate::import('synapseclient')
syn.utils <- reticulate::import('synapseutils')
syn <- syn.client$Synapse()
library(igraph)
library(tidyverse)

v_info <- function(net){
  v.attr <- tibble( na = igraph::vertex.attributes(net) ) %>% t() %>% 
    as_tibble(rownames = NA, .name_repair = ~igraph::vertex_attr_names(net) ) %>% 
    unnest(everything()) 
}


# download biodomain info -------------------------------------------------

# biological domain annotations
biodom <- syn$get('syn25428992') %>% .$path %>% readRDS() 
dom.cols <- syn$get('syn26856828') %>% .$path %>% read_csv( col_types=cols() )

biodom <- full_join(biodom, dom.cols, by = c('Biodomain'='domain') ) %>%
  mutate(
    Biodomain = if_else(Biodomain == 'none', NA_character_, Biodomain),
    Subdomain = if_else(is.na(Subdomain), 'none', Subdomain),
    sd = str_c(abbr,'_',Subdomain)
  )

# biodomain term annotation graph
bd.graph <- syn$get('syn51497240') %>% .$path %>% read_graph(format = 'graphml')
bd.graph.f <- delete_edges(bd.graph, edges = E(bd.graph)[E(bd.graph)$kappa < 0.5])
rm(bd.graph)


# summarize biodomains per gene -------------------------------------------

bd.genes <- syn$get('syn25428992') %>% .$path %>% readRDS(.) %>%
  select(Biodomain:Subdomain, hgnc_symbol) %>%
  unnest_longer(hgnc_symbol) %>%
  distinct() %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != '') %>%
  # NUMBER OF BIODOM/SUBDOM per GENE ##
  group_by(hgnc_symbol) %>%
  summarise(n_bd = length(unique(Biodomain)),
            n_sd = length(unique(paste0(Biodomain,' ',Subdomain)))) %>%
  ungroup() %>%
  
  # NUMBER OF BIODOM TERMS per GENE ##
  left_join(
    .,
    syn$get('syn25428992') %>% .$path %>% readRDS(.) %>%
      select(Biodomain:Subdomain, GO_ID, hgnc_symbol) %>%
      unnest_longer(hgnc_symbol) %>%
      filter(!is.na(hgnc_symbol), hgnc_symbol != '') %>%
      group_by(hgnc_symbol, Biodomain) %>%
      summarise(n_bd_term = length(unique(GO_ID))) %>%
      ungroup()
    ,
    by = 'hgnc_symbol'
  ) %>%
  
  # NUMBER OF SUBDOM TERMS per GENE ##
  left_join(
    .,
    syn$get('syn25428992') %>% .$path %>% readRDS(.) %>%
      select(Biodomain:Subdomain, GO_ID, hgnc_symbol) %>%
      unnest_longer(hgnc_symbol) %>%
      filter(!is.na(hgnc_symbol), hgnc_symbol != '') %>%
      group_by(hgnc_symbol, Biodomain, Subdomain) %>%
      summarise(n_sd_term = length(unique(GO_ID))) %>%
      ungroup()
    ,
    by = c('hgnc_symbol','Biodomain')
  ) %>%  
  
  # TOTAL NUMBER of TERMS per BIODOMAIN ##
  left_join(
    .,
    syn$get('syn25428992') %>% .$path %>% readRDS(.) %>%
      select(Biodomain, GO_ID) %>%
      group_by(Biodomain) %>%
      summarise(bd_terms = length(unique(GO_ID))) %>%
      ungroup()
    ,
    by = 'Biodomain'
  ) %>%
  mutate(bd_pct = 100 * (n_bd_term / bd_terms)) %>%
  
  # TOTAL NUMBER of TERMS per SUBDOMAIN ##
  left_join(
    .,
    syn$get('syn25428992') %>% .$path %>% readRDS(.) %>%
      select(Biodomain, Subdomain, GO_ID) %>%
      group_by(Biodomain,Subdomain) %>%
      summarise(sd_terms = length(unique(GO_ID))) %>%
      ungroup()
    ,
    by = c('Biodomain','Subdomain')
  ) %>%
  mutate(sd_pct = 100 * (n_sd_term / sd_terms)) %>%
  
  # Join with other info ##
  left_join(
    .,
    syn$get('syn25428992') %>% .$path %>% readRDS(.) %>%
      select(Biodomain, Subdomain, hgnc_symbol, n_hgncSymbol, GO_ID, GOterm_Name) %>%
      unnest_longer(hgnc_symbol),
    by = c('hgnc_symbol', 'Biodomain', 'Subdomain')
  )


# Compute top BD per gene -------------------------------------------------

gene.top.bd <- tibble(gene = unique(bd.genes$hgnc_symbol)) %>% 
  
  # Collect the table of terms per gene, the term annotation graph per gene, and extract term nodes per graph
  mutate(
    terms = furrr::future_map(gene, ~ bd.genes %>% filter(hgnc_symbol == .x)),
    graph = furrr::future_map(terms, ~ igraph::delete_vertices(bd.graph.f, igraph::V(bd.graph.f)[
      (igraph::V(bd.graph.f)$name %in% .x$GO_ID)==F])),
    nodes = furrr::future_map(graph, ~v_info(.x) %>% mutate(degree = degree(.x)) %>% 
                                relocate(degree) %>% arrange(desc(degree)))
  ) %>% 
  select(gene, nodes) %>% 
  unnest(nodes) %>% 
  select(hgnc_symbol = gene, degree, GO_ID = name) %>% 
  left_join(bd.genes, ., by = c('hgnc_symbol','GO_ID')) %>% 
  group_by(hgnc_symbol) %>% 
  
  # Calculate metrics to rank BIODOMAIN terms per GENE
  mutate(
    degree_rank = rank(degree, ties.method = 'min'),
    n_sym_rank = rank(1/n_hgncSymbol),
    bd_rank = bd_pct, #rank(n_bd_term),
    sd_rank = sd_pct,
    sum = degree_rank + bd_rank + 2*n_sym_rank
  ) %>% arrange(desc(sum)) %>% 
  
  # Nest all terms for each GENE
  nest(terms = Biodomain:sum) %>% 
  
  # Pull top term for each GENE
  mutate(top = map(terms, ~.x %>% slice(1))) %>% 
  unnest(top)

