#!/usr/bin/env bash



Rscript 01-Top_enrichments.R  \
  -i "Plasma_GFAP,Plasma_NFL,PlasmaPTau181,Plasma_pTau217,PlasmaAB140P,PlasmaAB142P,PlasmaABRatio,Plasma_pTau217_AB42Ratio" \
  -n "dataProt_SS-205063.PC.related.all_bios" 
