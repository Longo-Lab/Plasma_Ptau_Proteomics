#!/usr/bin/env bash

  # "Merged.M13"
  # "Merged.M8"
  # "Merged.M27"
  # "Merged.M5"
  # "Merged.M2"
  # "Merged.M29"
  # "Merged.M1"
  # "Merged.M32"
  # "Merged.M15"
  # "Merged.M18"
  # "Merged.M21"
  # "Merged.M14"

# plotting dists
Rscript 01-Prediction_power.R  \
  -i "Merged.M13,Merged.M8,Merged.M27,Merged.M5,Merged.M2,Merged.M29,Merged.M1,Merged.M32,Merged.M15,Merged.M18,Merged.M21,Merged.M14"

# biodomain enrichments
Rscript 02-Combined-TREAT-AD_GSEA.R  \
  -i "Merged.M13,Merged.M8,Merged.M27,Merged.M5,Merged.M2,Merged.M29,Merged.M1,Merged.M32,Merged.M15,Merged.M18,Merged.M21,Merged.M14" \
  --inpath "../../WSS_predict" \
  -n "dataProt_SS-205063.follow-up" \
  -s "top10.pred"

