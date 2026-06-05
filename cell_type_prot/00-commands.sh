#!/usr/bin/env bash

# 41593_2015_BFnn4160_MOESM43_ESM.xlsx
# Supplementary Table 6
# Differentially expressed proteins in cultured CNS cell types. (XLSX 2671 kb)

# 41593_2015_BFnn4160_MOESM50_ESM.xlsx
# Supplementary Table 13
# Differentially expressed proteins in isolated CNS cell types. (XLSX 3391 kb)

# Rscript 01-proteomic_celltypes.R

# # correlation-based cell type enrichment
# Rscript 01-Sharma_etal_correlation.R  \
  # -i "PlasmaGFAP,PlasmaNfL,PlasmaPTau181,PlasmapTau217,PlasmaAB140P,PlasmaAB142P,PlasmaABRatio,PlasmapTau217_AB42Ratio" \
  # -r "../../cell_type_prot" \
  # -n "dataProt_SS-205063.PC.related.all_bios" 

# Rscript 02-correlation_topgenes.R \
  # -i "../../baseline_PCs_relatedness" \
  # -n "dataProt_SS-205063.PC.related" \
  # -r "../../cell_type_prot" \
  # -b PlasmaGFAP
  
# Rscript 02-correlation_topgenes.R \
  # -i "../../baseline_PCs_relatedness" \
  # -n "dataProt_SS-205063.PC.related" \
  # -r "../../cell_type_prot" \
  # --rtype "cultured" \
  # -b PlasmaGFAP

# Rscript 02-correlation_topgenes.R \
  # -i "../../baseline_PCs_relatedness" \
  # -n "dataProt_SS-205063.PC.related" \
  # -r "../../cell_type_prot" \
  # --rtype "specificity" \
  # -b PlasmaGFAP
  
# Rscript 02-correlation_topgenes.R \
  # -i "../../baseline_PCs_relatedness" \
  # -n "dataProt_SS-205063.PC.related" \
  # -r "../../cell_type_prot" \
  # -b PlasmaNfL

# Rscript 02-correlation_topgenes.R \
  # -i "../../baseline_PCs_relatedness" \
  # -n "dataProt_SS-205063.PC.related" \
  # -r "../../cell_type_prot" \
  # -b PlasmaPTau181

# Rscript 02-correlation_topgenes.R \
  # -i "../../baseline_PCs_relatedness" \
  # -n "dataProt_SS-205063.PC.related" \
  # -r "../../cell_type_prot" \
  # --rtype "cultured" \
  # -b PlasmapTau217

# EWCE ----------------------------------------------------------------
# in lieu of the normalized and scaled expression, instead using cell type specificity 
# to generate scores for samples. These were generated on SCG and copied over


# Specificity based cell type enrichment
Rscript 03-EWCE_enrichment.R  \
  -i "PlasmaGFAP,PlasmaNfL,PlasmaPTau181,PlasmapTau217,PlasmaAB140P,PlasmaAB142P,PlasmaABRatio,PlasmapTau217_AB42Ratio" \
  -r "../../cell_type_prot" \
  -t 16

# SEA-AD Neurons
parallel -j1 Rscript 04-EWCE_topgenes.R \
  -n "dataProt_SS-205063.PC.related.all_bios" \
  -c SEA-AD.Class \
  -g stat_padj \
  -t {1} \
  -l 2 \
  -r "ewce_SEA-AD.Class.1kDS.rda" \
  ::: "PlasmaGFAP" "PlasmaPTau181" "PlasmaAB140P" "PlasmaAB142P" "PlasmapTau217"


# SEA-AD Neurons
parallel -j1 Rscript 04-EWCE_topgenes.R \
  -n "dataProt_SS-205063.PC.related.all_bios" \
  -c SEA-AD.Neighborhood \
  -g stat_padj \
  -t {1} \
  -l 2 \
  -r "ewce_SEA-AD.Neighborhood.1kDS.rda" \
  ::: "PlasmaGFAP" "PlasmaPTau181" "PlasmaAB140P" "PlasmaAB142P" "PlasmaNfL"

# Run again for AB40, no labels showed
parallel -j1 Rscript 04-EWCE_topgenes.R \
  -n "dataProt_SS-205063.PC.related.all_bios" \
  -c SEA-AD.Neighborhood \
  -g stat_padj \
  -t {1} \
  -l 2 \
  -m 32 \
  -r "ewce_SEA-AD.Neighborhood.1kDS.rda" \
  ::: "PlasmaAB140P"



# Sharma et al types
parallel -j1 Rscript 04-EWCE_topgenes.R \
  -n "dataProt_SS-205063.PC.related.all_bios" \
  -c Sharma_etal.cultured \
  -g stat_padj \
  -t {1} \
  -r "Sharma_etal.EWCE.human.rda" \
  ::: "PlasmaGFAP" "PlasmaPTau181" "PlasmapTau217" "PlasmaAB140P"

# Plot heatmaps
Rscript 05-EWCE_heatmap.R  \
  -i "PlasmaGFAP,PlasmaNfL,PlasmaPTau181,PlasmapTau217,PlasmaAB140P,PlasmaAB142P,PlasmaABRatio,PlasmapTau217_AB42Ratio" \
  -n "dataProt_SS-205063.PC.related.all_bios" \
  -w 3.75 \
  -c "Sharma_etal.cultured" 

Rscript 05-EWCE_heatmap.R  \
  -i "PlasmaGFAP,PlasmaNfL,PlasmaPTau181,PlasmapTau217,PlasmaAB140P,PlasmaAB142P,PlasmaABRatio,PlasmapTau217_AB42Ratio" \
  -n "dataProt_SS-205063.PC.related.all_bios" \
  -l 2 \
  -w 5.5 \
  -c "SEA-AD.Neighborhood" 



# WSS_proteins ---------------------------------------------------------------
Rscript 06-EWCE_enrichment_WSS.R  \
  -i "Merged.M1,Merged.M11,Merged.M7,Merged.M9,Merged.M2,PlasmaPTau181.M8,Merged.M6,PlasmaAB140P.M17,Merged.M12,Merged.M4,Merged.M5,PlasmaAB140P.M20,Merged.M8,PlasmaPTau181.M35,Merged.M3,Merged.M26,Merged.M13,Merged.M10,PlasmaPTau181.M18,Merged.M14,Merged.M19,Merged.M16,PlasmaPTau181.M26,Merged.M29,Merged.M24,PlasmaPTau181.M13,PlasmaPTau181.M47,Merged.M25,PlasmaAB140P.M35,PlasmaPTau181.M46,Merged.M20,Merged.M32,Merged.M22,Merged.M15,Merged.M21,Merged.M18,PlasmaPTau181.M5,PlasmaAB140P.M15,Merged.M27,PlasmaAB140P.M40,Merged.M17,PlasmaNfL.M3,PlasmaAB140P.M18,PlasmaPTau181.M30,PlasmaGFAP.M12" \
  -n "dataProt_SS-205063.baseline_WSS" \
  -m 10 \
  -s "all_WSS"


# Plot heatmaps
Rscript 05-EWCE_heatmap.R  \
  -i "Merged.M1,Merged.M11,Merged.M7,Merged.M9,Merged.M2,PlasmaPTau181.M8,Merged.M6,PlasmaAB140P.M17,Merged.M12,Merged.M4,Merged.M5,PlasmaAB140P.M20,Merged.M8,PlasmaPTau181.M35,Merged.M3,Merged.M26,Merged.M13,Merged.M10,PlasmaPTau181.M18,Merged.M14,Merged.M19,Merged.M16,PlasmaPTau181.M26,Merged.M29,Merged.M24,PlasmaPTau181.M13,PlasmaPTau181.M47,Merged.M25,PlasmaAB140P.M35,PlasmaPTau181.M46,Merged.M20,Merged.M32,Merged.M22,Merged.M15,Merged.M21,Merged.M18,PlasmaPTau181.M5,PlasmaAB140P.M15,Merged.M27,PlasmaAB140P.M40,Merged.M17,PlasmaNfL.M3,PlasmaAB140P.M18,PlasmaPTau181.M30,PlasmaGFAP.M12" \
  -n "dataProt_SS-205063.baseline_WSS.all_WSS" \
  -g "WSS_mod" \
  -r "../module_scores/STRING.INT.WSS.terms.csv" \
  -w 5 \
  --height 5 \
  -c "Sharma_etal.cultured" 

Rscript 05-EWCE_heatmap.R  \
  -i "Merged.M1,Merged.M11,Merged.M7,Merged.M9,Merged.M2,PlasmaPTau181.M8,Merged.M6,PlasmaAB140P.M17,Merged.M12,Merged.M4,Merged.M5,PlasmaAB140P.M20,Merged.M8,PlasmaPTau181.M35,Merged.M3,Merged.M26,Merged.M13,Merged.M10,PlasmaPTau181.M18,Merged.M14,Merged.M19,Merged.M16,PlasmaPTau181.M26,Merged.M29,Merged.M24,PlasmaPTau181.M13,PlasmaPTau181.M47,Merged.M25,PlasmaAB140P.M35,PlasmaPTau181.M46,Merged.M20,Merged.M32,Merged.M22,Merged.M15,Merged.M21,Merged.M18,PlasmaPTau181.M5,PlasmaAB140P.M15,Merged.M27,PlasmaAB140P.M40,Merged.M17,PlasmaNfL.M3,PlasmaAB140P.M18,PlasmaPTau181.M30,PlasmaGFAP.M12" \
  -n "dataProt_SS-205063.baseline_WSS.all_WSS" \
  -g "WSS_mod" \
  -l 2 \
  -r "../module_scores/STRING.INT.WSS.terms.csv" \
  -w 6.5 \
  --height 7 \
  -c "SEA-AD.Neighborhood" 
