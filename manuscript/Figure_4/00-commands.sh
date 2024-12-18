#!/usr/bin/env bash


apoe=("" ".APOE")

parallel -j 2 Rscript 01-Charma_etal_correlation.R  \
  -i {1} \
  -r "../../cell_type_prot"
  -n "dataProt_SS-205063.PC.related.all_bios{2}" \
  ::: "Plasma_GFAP,Plasma_NFL,PlasmaPTau181,Plasma_pTau217,PlasmaAB140P,PlasmaAB142P,PlasmaABRatio,Plasma_pTau217_AB42Ratio" \
      "Plasma_GFAP.APOE,Plasma_NFL.APOE,PlasmaPTau181.APOE,Plasma_pTau217.APOE,PlasmaAB140P.APOE,PlasmaAB142P.APOE,PlasmaABRatio.APOE,Plasma_pTau217_AB42Ratio.APOE" \
  :::+ "${apoe[@]}"

