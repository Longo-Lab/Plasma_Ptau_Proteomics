#!/usr/bin/env bash

# run on SCG and results copied over ----------------------------------------

# # grabbing up to date TREAT-AD gene sets
# synapse -p auth_token get syn26592124
# synapse -p auth_token get syn25428992
# synapse get -q "SELECT * FROM syn22758536"



# new version of gene lookups using 01-all_gene_lists.r ---------------------
#  uses the top hit for one2many/many2many and ensembl version 111 directly
ml R/4.2.2
Rscript 01-all_gene_lists.R 

# updated genes_info to v111
sbatch -J genes_info --mem=5G -c 2 -t 01:00:00 -p interactive \
  -o %A_genes_info_v111.log \
  --wrap "ml R/4.2.2; Rscript 02-genes_info.R 111"


