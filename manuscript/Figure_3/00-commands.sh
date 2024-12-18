#!/usr/bin/env bash


apoe=("" ".APOE")
modules=("TMT-AD" "TMT-Resilience" "TREAT-AD" "Wan_etal")

parallel -j 8 Rscript 01-Combined-Enrich.R  \
  -i {1} \
  -n "dataProt_SS-205063.PC.related.all_bios{2}" \
  -m "{3}" \
  -r "../../biodomains" \
  ::: "Plasma_GFAP,Plasma_NFL,PlasmaPTau181,Plasma_pTau217,PlasmaAB140P,PlasmaAB142P,PlasmaABRatio,Plasma_pTau217_AB42Ratio" \
      "Plasma_GFAP.APOE,Plasma_NFL.APOE,PlasmaPTau181.APOE,Plasma_pTau217.APOE,PlasmaAB140P.APOE,PlasmaAB142P.APOE,PlasmaABRatio.APOE,Plasma_pTau217_AB42Ratio.APOE" \
  :::+ "${apoe[@]}" \
  ::: "${modules[@]}"


# another set with csf/plasma proteomics
modules=("TMTSoma-CSF" "Soma-plasma" "TMT-PTXtrial")
parallel -j 8 Rscript 01-Combined-Enrich.R  \
  -i {1} \
  -n "dataProt_SS-205063.PC.related.all_bios{2}" \
  -m "{3}" \
  -r "../../biodomains" \
  ::: "Plasma_GFAP,Plasma_NFL,PlasmaPTau181,Plasma_pTau217,PlasmaAB140P,PlasmaAB142P,PlasmaABRatio,Plasma_pTau217_AB42Ratio" \
      "Plasma_GFAP.APOE,Plasma_NFL.APOE,PlasmaPTau181.APOE,Plasma_pTau217.APOE,PlasmaAB140P.APOE,PlasmaAB142P.APOE,PlasmaABRatio.APOE,Plasma_pTau217_AB42Ratio.APOE" \
  :::+ "${apoe[@]}" \
  ::: "${modules[@]}"
