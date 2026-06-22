#!/usr/bin/env bash


# # plotting dists

# Rscript 01-Long_plots.R  \
  # -i "Adaptive immunity|Immunoglobulins|Cytokines,Amino acid metabolism,Cadherin|Cell adhesion,Carbon metabolism,Cell cycle,Chaperone|Protein folding,Collagen,Complement cascade,Dementia|HDL,ER-Golgi transport,Ephrin signaling,Extracellular space,Focal adhesion|Integrins,Glutathione metabolism|Oxidative stress,Glycosylation,HDL,NFkB|APC/C degradation,Neurexins|Neuroligins,Neuroactive ligand-receptor interaction,PI3K-Akt|JAK-STAT|MAPK signaling,Protein Ubiquitination,RHO GTPase|Neurotrophin signaling,RNA-binding,Synaptic vesicle|SNARE complex,Translation,Wnt signaling" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "STRING_alpha_scaled" \
  # -w 20 \
  # -b 6

Rscript 01-Long_plots.R  \
  -i "Merged.M1,Merged.M11,Merged.M7,Merged.M9,Merged.M2,PlasmaPTau181.M8,Merged.M6,PlasmaAB140P.M17,Merged.M12,Merged.M4,Merged.M5,PlasmaAB140P.M20,Merged.M8,PlasmaPTau181.M35,Merged.M3,Merged.M26,Merged.M13,Merged.M10,PlasmaPTau181.M18,Merged.M14,Merged.M19,Merged.M16,PlasmaPTau181.M26,Merged.M29,Merged.M24,PlasmaPTau181.M13,PlasmaPTau181.M47,Merged.M25,PlasmaAB140P.M35,PlasmaPTau181.M46,Merged.M20,Merged.M32,Merged.M22,Merged.M15,Merged.M21,Merged.M18,PlasmaPTau181.M5,PlasmaAB140P.M15,Merged.M27,PlasmaAB140P.M40,Merged.M17,PlasmaNfL.M3,PlasmaAB140P.M18,PlasmaPTau181.M30,PlasmaGFAP.M12" \
  -n "dataProt_SS-205063.baseline.WSS" \
  -s "STRING.INT.WSS" \
  -w 28 \
  -b 12 




# # Trying with ranksum --------------------------------------
# scores=('ranksum_scaled' 'ranksum_spec')
# for s in ${scores[*]}; do
  # Rscript 01-Long_plots.R  \
    # -i "APP Metabolism,Apoptosis,Autophagy,Cell Cycle,DNA Repair,Endolysosome,Epigenetic,Immune Response,Lipid Metabolism,Metal Binding and Homeostasis,Mitochondrial Metabolism,Myelination,Oxidative Stress,Proteostasis,Structural Stabilization,Synapse,Tau Homeostasis,Vasculature" \
    # -n "dataProt_SS-205063.baseline.all_bios" \
    # -s "treatad_${s}"
# done




# WSS heatmap ------------------------------
Rscript 02-Combined-TREAT-AD_GSEA.R  \
  -i "Merged.M1,Merged.M11,Merged.M7,Merged.M9,Merged.M2,PlasmaPTau181.M8,Merged.M6,PlasmaAB140P.M17,Merged.M12,Merged.M4,Merged.M5,PlasmaAB140P.M20,Merged.M8,PlasmaPTau181.M35,Merged.M3,Merged.M26,Merged.M13,Merged.M10,PlasmaPTau181.M18,Merged.M14,Merged.M19,Merged.M16,PlasmaPTau181.M26,Merged.M29,Merged.M24,PlasmaPTau181.M13,PlasmaPTau181.M47,Merged.M25,PlasmaAB140P.M35,PlasmaPTau181.M46,Merged.M20,Merged.M32,Merged.M22,Merged.M15,Merged.M21,Merged.M18,PlasmaPTau181.M5,PlasmaAB140P.M15,Merged.M27,PlasmaAB140P.M40,Merged.M17,PlasmaNfL.M3,PlasmaAB140P.M18,PlasmaPTau181.M30,PlasmaGFAP.M12" \
  -n "dataProt_SS-205063.baseline_WSS" \
  -s "all_WSS"

# WSS in GNPC data...
Rscript 02-Combined-TREAT-AD_GSEA.R  \
  -i "Merged.M1,Merged.M11,Merged.M7,Merged.M9,Merged.M2,PlasmaPTau181.M8,Merged.M6,PlasmaAB140P.M17,Merged.M12,Merged.M4,Merged.M5,PlasmaAB140P.M20,Merged.M8,PlasmaPTau181.M35,Merged.M3,Merged.M26,Merged.M13,Merged.M10,PlasmaPTau181.M18,Merged.M14,Merged.M19,Merged.M16,PlasmaPTau181.M26,Merged.M29,Merged.M24,PlasmaPTau181.M13,PlasmaPTau181.M47,Merged.M25,PlasmaAB140P.M35,PlasmaPTau181.M46,Merged.M20,Merged.M32,Merged.M22,Merged.M15,Merged.M21,Merged.M18,PlasmaPTau181.M5,PlasmaAB140P.M15,Merged.M27,PlasmaAB140P.M40,Merged.M17,PlasmaNfL.M3,PlasmaAB140P.M18,PlasmaPTau181.M30,PlasmaGFAP.M12" \
  --inpath ../../gnpc/export/WSS_prot \
  -n "GNPC.baseline_prot" \
  -s "all_WSS"


Rscript 03-Focused-TREAT-AD_GSEA.R  \
  -m "Merged.M8" \
  -n "dataProt_SS-205063.baseline_WSS" \
  -s "all_WSS"

