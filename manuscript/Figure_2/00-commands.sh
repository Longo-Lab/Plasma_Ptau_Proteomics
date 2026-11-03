#!/usr/bin/env bash



Rscript 01-Top_enrichments.R  \
  -i "PlasmaGFAP,PlasmaNFL,PlasmaPTau181,PlasmapTau217,PlasmaAB140P,PlasmaAB142P,PlasmaABRatio,PlasmapTau217_AB42Ratio" \
  -n "dataProt_SS-205063.PC.related.all_bios" 

Rscript 02-combine_STRING_modules.R

