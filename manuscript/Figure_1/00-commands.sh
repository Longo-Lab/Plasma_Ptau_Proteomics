#!/usr/bin/env bash


apoe=("" ".APOE")

parallel -j 2 Rscript 01-Top_genes.R  \
  -i {1} \
  -n "dataProt_SS-205063.PC.related.all_bios{2}" \
  ::: "PlasmaGFAP,PlasmaNfL,PlasmaPTau181,PlasmapTau217,PlasmaAB140P,PlasmaAB142P,PlasmaABRatio,PlasmapTau217_AB42Ratio" \
      "PlasmaGFAP.APOE,PlasmaNfL.APOE,PlasmaPTau181.APOE,PlasmapTau217.APOE,PlasmaAB140P.APOE,PlasmaAB142P.APOE,PlasmaABRatio.APOE,PlasmapTau217_AB42Ratio.APOE" \
  :::+ "${apoe[@]}"

parallel -j 2 Rscript 02-Overlap_genes.R  \
  -i {1} \
  -n "dataProt_SS-205063.PC.related.all_bios{2}" \
  ::: "PlasmaGFAP,PlasmaNfL,PlasmaPTau181,PlasmapTau217,PlasmaAB140P,PlasmaAB142P,PlasmaABRatio,PlasmapTau217_AB42Ratio" \
      "PlasmaGFAP.APOE,PlasmaNfL.APOE,PlasmaPTau181.APOE,PlasmapTau217.APOE,PlasmaAB140P.APOE,PlasmaAB142P.APOE,PlasmaABRatio.APOE,PlasmapTau217_AB42Ratio.APOE" \
  :::+ "${apoe[@]}"

parallel -j 2 Rscript 03-Long_plots.R  \
  -i {1} \
  -n "dataProt_SS-205063.PC.related.all_bios{2}" \
  ::: "PlasmaGFAP,PlasmaNfL,PlasmaPTau181,PlasmapTau217,PlasmaAB140P,PlasmaAB142P,PlasmaABRatio,PlasmapTau217_AB42Ratio" \
      "PlasmaGFAP.APOE,PlasmaNfL.APOE,PlasmaPTau181.APOE,PlasmapTau217.APOE,PlasmaAB140P.APOE,PlasmaAB142P.APOE,PlasmaABRatio.APOE,PlasmapTau217_AB42Ratio.APOE" \
  :::+ "${apoe[@]}"
