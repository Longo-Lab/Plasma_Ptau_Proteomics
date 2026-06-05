#!/usr/bin/env bash


apoe=("" ".APOE")
modules=("TMT-AD" "Wan_etal" "TMTSoma-CSF" "Soma-plasma")

parallel -j 8 Rscript 01-Combined-Enrich.R  \
  -i {1} \
  -n "dataProt_SS-205063.PC.related.all_bios{2}" \
  -m "{3}" \
  -r "../../biodomains" \
  ::: "PlasmaGFAP,PlasmaNfL,PlasmaPTau181,PlasmapTau217,PlasmaAB140P,PlasmaAB142P,PlasmaABRatio,PlasmapTau217_AB42Ratio" \
      "PlasmaGFAP.APOE,PlasmaNfL.APOE,PlasmaPTau181.APOE,PlasmapTau217.APOE,PlasmaAB140P.APOE,PlasmaAB142P.APOE,PlasmaABRatio.APOE,PlasmapTau217_AB42Ratio.APOE" \
  :::+ "${apoe[@]}" \
  ::: "${modules[@]}"

Rscript 02-Long_plots.R \
  -i "ARL2.12587.65.3,S100A13.7223.60.3,TBCA.12501.10.3" \
  -n "dataProt_SS-205063.PC.related.ARL2-TBCA-S100A13"

# Biomarker heatmap -----------------------
Rscript 03-Combined-TREAT-AD_GSEA.R  \
  -i "PlasmaGFAP,PlasmaNfL,PlasmaPTau181,PlasmapTau217,PlasmaAB140P,PlasmaAB142P,PlasmaABRatio,PlasmapTau217_AB42Ratio" \
  --inpath "../../baseline_PCs_relatedness" \
  -r "../../biomarker_names.txt" \
  -s "all_bios"


