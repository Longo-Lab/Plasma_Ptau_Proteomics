#!/usr/bin/env bash

# filter genesets
Rscript 01-gene_lists.R \
  -i "Plasma_GFAP,Plasma_NFL,PlasmaPTau181,Plasma_pTau217,PlasmaAB140P,PlasmaAB142P,PlasmaABRatio,Plasma_pTau217_AB42Ratio"

# run padj gene set analysis
bash 02-magma_command_padj.sh
bash 03-magma_command_top10.sh

# run continuous property heritability enrichment
bash 04-magma_command_cont_prop.sh

Rscript 05-magma_plots.R
Rscript 06-magma_plots_continuous.R

