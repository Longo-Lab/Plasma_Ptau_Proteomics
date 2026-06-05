#!/usr/bin/env bash


# # plotting dists
# Rscript 01-Long_plots.R  \
  # -i "Astrocyte,Microglia,Neuron,Oligodendrocyte" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "Sharma.lvl1"

# Rscript 01-Long_plots.R  \
  # -i "Astrocyte,CGE,L2/3 IT,L4/5/6 IT-ET Car3,MGE,Microglia-PVM,NP/CT/L6b,Oligo,Vascular" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "SEA-AD.Neighborhood.1kDS.lvl2"

# Rscript 01-Long_plots.R  \
  # -i "APP Metabolism,Apoptosis,Autophagy,Cell Cycle,Endolysosome,Immune Response,Lipid Metabolism,Metal Binding and Homeostasis,Mitochondrial Metabolism,Myelination,Oxidative Stress,Proteostasis,Structural Stabilization,Synapse,Tau Homeostasis,Vasculature" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "treatad_BioDomain"

# scores=('harmonic_norm' 'harmonic_scaled' 'harmonic_spec' 'closeness_norm' 'alpha_scaled' 'alpha_spec' 'ranksum_scaled' 'ranksum_spec')
# for s in ${scores[*]}; do
  # Rscript 01-Long_plots.R  \
    # -i "APP Metabolism,Apoptosis,Autophagy,Cell Cycle,Endolysosome,Immune Response,Lipid Metabolism,Metal Binding and Homeostasis,Mitochondrial Metabolism,Myelination,Oxidative Stress,Proteostasis,Structural Stabilization,Synapse,Tau Homeostasis,Vasculature" \
    # -n "dataProt_SS-205063.baseline.all_bios" \
    # -s "treatad_${s}"
# done

# Rscript 01-Long_plots.R  \
  # -i "Adaptive immunity|Immunoglobulins|Cytokines,Amino acid metabolism,Cadherin|Cell adhesion,Carbon metabolism,Cell cycle,Chaperone|Protein folding,Collagen,Complement cascade,Dementia|HDL,ER-Golgi transport,Ephrin signaling,Extracellular space,Focal adhesion|Integrins,Glutathione metabolism|Oxidative stress,Glycosylation,HDL,NFkB|APC/C degradation,Neurexins|Neuroligins,Neuroactive ligand-receptor interaction,PI3K-Akt|JAK-STAT|MAPK signaling,Protein Ubiquitination,RHO GTPase|Neurotrophin signaling,RNA-binding,Synaptic vesicle|SNARE complex,Translation,Wnt signaling" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "STRING_alpha_scaled" \
  # -w 20 \
  # -b 6

# Rscript 01-Long_plots.R  \
  # -i "MEMerged-M1,MEMerged-M10,MEMerged-M11,MEMerged-M12,MEMerged-M13,MEMerged-M14,MEMerged-M15,MEMerged-M16,MEMerged-M17,MEMerged-M18,MEMerged-M19,MEMerged-M2,MEMerged-M20,MEMerged-M21,MEMerged-M22,MEMerged-M24,MEMerged-M25,MEMerged-M26,MEMerged-M27,MEMerged-M28,MEMerged-M29,MEMerged-M3,MEMerged-M30,MEMerged-M32,MEMerged-M33,MEMerged-M35,MEMerged-M36,MEMerged-M4,MEMerged-M5,MEMerged-M6,MEMerged-M7,MEMerged-M8,MEMerged-M9,MEPlasmaAB140P-M15,MEPlasmaAB140P-M17,MEPlasmaAB140P-M18,MEPlasmaAB140P-M22,MEPlasmaAB140P-M24,MEPlasmaAB140P-M28,MEPlasmaAB140P-M38,MEPlasmaAB140P-M40,MEPlasmaAB140P-M41,MEPlasmaAB140P-M44,MEPlasmaAB140P-M45,MEPlasmaAB140P-M48,MEPlasmaAB140P-M50,MEPlasmaAB142P-M25,MEPlasmaAB142P-M26,MEPlasmaAB142P-M28,MEPlasmaAB142P-M31,MEPlasmaAB142P-M34,MEPlasmaAB142P-M36,MEPlasmaGFAP-M12,MEPlasmaGFAP-M4,MEPlasmaGFAP-M9,MEPlasmaNfL-M21,MEPlasmaNfL-M25,MEPlasmaNfL-M26,MEPlasmaNfL-M3,MEPlasmaNfL-M32,MEPlasmaNfL-M33,MEPlasmaPTau181-M13,MEPlasmaPTau181-M18,MEPlasmaPTau181-M20,MEPlasmaPTau181-M23,MEPlasmaPTau181-M24,MEPlasmaPTau181-M26,MEPlasmaPTau181-M27,MEPlasmaPTau181-M30,MEPlasmaPTau181-M35,MEPlasmaPTau181-M37,MEPlasmaPTau181-M40,MEPlasmaPTau181-M43,MEPlasmaPTau181-M46,MEPlasmaPTau181-M47,MEPlasmaPTau181-M5,MEPlasmaPTau181-M8,MEPlasmapTau217-M3" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "STRING.eigengenes" \
  # -w 28 \
  # -b 13

# Rscript 01-Long_plots.R  \
  # -i "MEMerged-M1,MEMerged-M10,MEMerged-M11,MEMerged-M12,MEMerged-M13,MEMerged-M14,MEMerged-M15,MEMerged-M16,MEMerged-M17,MEMerged-M18,MEMerged-M19,MEMerged-M2,MEMerged-M20,MEMerged-M21,MEMerged-M22,MEMerged-M24,MEMerged-M25,MEMerged-M26,MEMerged-M27,MEMerged-M28,MEMerged-M29,MEMerged-M3,MEMerged-M30,MEMerged-M32,MEMerged-M33,MEMerged-M35,MEMerged-M36,MEMerged-M4,MEMerged-M5,MEMerged-M6,MEMerged-M7,MEMerged-M8,MEMerged-M9,MEPlasmaAB140P-M15,MEPlasmaAB140P-M17,MEPlasmaAB140P-M18,MEPlasmaAB140P-M22,MEPlasmaAB140P-M24,MEPlasmaAB140P-M28,MEPlasmaAB140P-M38,MEPlasmaAB140P-M40,MEPlasmaAB140P-M41,MEPlasmaAB140P-M44,MEPlasmaAB140P-M45,MEPlasmaAB140P-M48,MEPlasmaAB140P-M50,MEPlasmaAB142P-M25,MEPlasmaAB142P-M26,MEPlasmaAB142P-M28,MEPlasmaAB142P-M31,MEPlasmaAB142P-M34,MEPlasmaAB142P-M36,MEPlasmaGFAP-M12,MEPlasmaGFAP-M4,MEPlasmaGFAP-M9,MEPlasmaNfL-M21,MEPlasmaNfL-M25,MEPlasmaNfL-M26,MEPlasmaNfL-M3,MEPlasmaNfL-M32,MEPlasmaNfL-M33,MEPlasmaPTau181-M13,MEPlasmaPTau181-M18,MEPlasmaPTau181-M20,MEPlasmaPTau181-M23,MEPlasmaPTau181-M24,MEPlasmaPTau181-M26,MEPlasmaPTau181-M27,MEPlasmaPTau181-M30,MEPlasmaPTau181-M35,MEPlasmaPTau181-M37,MEPlasmaPTau181-M40,MEPlasmaPTau181-M43,MEPlasmaPTau181-M46,MEPlasmaPTau181-M47,MEPlasmaPTau181-M5,MEPlasmaPTau181-M8,MEPlasmapTau217-M3" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "STRING.INT.eigengenes" \
  # -w 28 \
  # -b 13 

# Rscript 01-Long_plots.R  \
  # -i "MEcelltype-Astrocyte,MEcelltype-CGE,MEcelltype-L2/3 IT,MEcelltype-L4/5/6 IT-ET Car3,MEcelltype-MGE,MEcelltype-Microglia-PVM,MEcelltype-NP/CT/L6b,MEcelltype-Oligo,MEcelltype-Vascular" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "SEA-AD.Neighborhood_0.9.INT.eigengenes" 

# Rscript 01-Long_plots.R  \
  # -i "MEcelltype-Astrocyte,MEcelltype-Microglia,MEcelltype-Neuron,MEcelltype-Oligodendrocyte" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "Sharma_etal_0.9.INT.eigengenes" 

Rscript 01-Long_plots.R  \
  -i "Merged.M1,Merged.M11,Merged.M7,Merged.M9,Merged.M2,PlasmaPTau181.M8,Merged.M6,PlasmaAB140P.M17,Merged.M12,Merged.M4,Merged.M5,PlasmaAB140P.M20,Merged.M8,PlasmaPTau181.M35,Merged.M3,Merged.M26,Merged.M13,Merged.M10,PlasmaPTau181.M18,Merged.M14,Merged.M19,Merged.M16,PlasmaPTau181.M26,Merged.M29,Merged.M24,PlasmaPTau181.M13,PlasmaPTau181.M47,Merged.M25,PlasmaAB140P.M35,PlasmaPTau181.M46,Merged.M20,Merged.M32,Merged.M22,Merged.M15,Merged.M21,Merged.M18,PlasmaPTau181.M5,PlasmaAB140P.M15,Merged.M27,PlasmaAB140P.M40,Merged.M17,PlasmaNfL.M3,PlasmaAB140P.M18,PlasmaPTau181.M30,PlasmaGFAP.M12" \
  -n "dataProt_SS-205063.baseline.WSS" \
  -s "STRING.INT.WSS" \
  -w 28 \
  -b 12 


# heatmaps
# Rscript 02-Clustering.R  \
  # -i "Astrocyte,Microglia,Neuron,Oligodendrocyte" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "Sharma.lvl1" \
  # -k 0

# Rscript 02-Clustering.R  \
  # -i "Astrocyte,CGE,L2/3 IT,L4/5/6 IT-ET Car3,MGE,Microglia-PVM,NP/CT/L6b,Oligo,Vascular" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "SEA-AD.Neighborhood.1kDS.lvl2" \
  # -k 10 # second best 4 proposed

# Rscript 02-Clustering.R  \
  # -i "APP Metabolism,Apoptosis,Autophagy,Cell Cycle,Endolysosome,Immune Response,Lipid Metabolism,Metal Binding and Homeostasis,Mitochondrial Metabolism,Myelination,Oxidative Stress,Proteostasis,Structural Stabilization,Synapse,Tau Homeostasis,Vasculature" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "treatad_BioDomain" \
  # -k 0

# scores=('harmonic_norm' 'harmonic_scaled' 'harmonic_spec' 'closeness_norm' 'alpha_scaled' 'alpha_spec')
# for s in ${scores[*]}; do
  # Rscript 02-Clustering.R  \
    # -i "APP Metabolism,Apoptosis,Autophagy,Cell Cycle,Endolysosome,Immune Response,Lipid Metabolism,Metal Binding and Homeostasis,Mitochondrial Metabolism,Myelination,Oxidative Stress,Proteostasis,Structural Stabilization,Synapse,Tau Homeostasis,Vasculature" \
    # -n "dataProt_SS-205063.baseline.all_bios" \
    # -s "treatad_${s}" \
    # -k 0
# done

# Rscript 02-Clustering.R  \
  # -i "Adaptive immunity|Immunoglobulins|Cytokines,Amino acid metabolism,Cadherin|Cell adhesion,Carbon metabolism,Cell cycle,Chaperone|Protein folding,Collagen,Complement cascade,Dementia|HDL,ER-Golgi transport,Ephrin signaling,Extracellular space,Focal adhesion|Integrins,Glutathione metabolism|Oxidative stress,Glycosylation,HDL,NFkB|APC/C degradation,Neurexins|Neuroligins,Neuroactive ligand-receptor interaction,PI3K-Akt|JAK-STAT|MAPK signaling,Protein Ubiquitination,RHO GTPase|Neurotrophin signaling,RNA-binding,Synaptic vesicle|SNARE complex,Translation,Wnt signaling" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "STRING_alpha_scaled" \
  # -k 0

# Rscript 02-Clustering.R  \
  # -i "Adaptive immunity|Immunoglobulins|Cytokines,Amino acid metabolism,Cadherin|Cell adhesion,Carbon metabolism,Cell cycle,Chaperone|Protein folding,Collagen,Complement cascade,Dementia|HDL,ER-Golgi transport,Ephrin signaling,Extracellular space,Focal adhesion|Integrins,Glutathione metabolism|Oxidative stress,Glycosylation,HDL,NFkB|APC/C degradation,Neurexins|Neuroligins,Neuroactive ligand-receptor interaction,PI3K-Akt|JAK-STAT|MAPK signaling,Protein Ubiquitination,RHO GTPase|Neurotrophin signaling,RNA-binding,Synaptic vesicle|SNARE complex,Translation,Wnt signaling" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "STRING_alpha_scaled" \
  # -k 15 # tied for third with 5

# Rscript 02-Clustering.R  \
  # -i "MEMerged-M1,MEMerged-M10,MEMerged-M11,MEMerged-M12,MEMerged-M13,MEMerged-M14,MEMerged-M15,MEMerged-M16,MEMerged-M17,MEMerged-M18,MEMerged-M19,MEMerged-M2,MEMerged-M20,MEMerged-M21,MEMerged-M22,MEMerged-M24,MEMerged-M25,MEMerged-M26,MEMerged-M27,MEMerged-M28,MEMerged-M29,MEMerged-M3,MEMerged-M30,MEMerged-M32,MEMerged-M33,MEMerged-M35,MEMerged-M36,MEMerged-M4,MEMerged-M5,MEMerged-M6,MEMerged-M7,MEMerged-M8,MEMerged-M9,MEPlasmaAB140P-M15,MEPlasmaAB140P-M17,MEPlasmaAB140P-M18,MEPlasmaAB140P-M22,MEPlasmaAB140P-M24,MEPlasmaAB140P-M28,MEPlasmaAB140P-M38,MEPlasmaAB140P-M40,MEPlasmaAB140P-M41,MEPlasmaAB140P-M44,MEPlasmaAB140P-M45,MEPlasmaAB140P-M48,MEPlasmaAB140P-M50,MEPlasmaAB142P-M25,MEPlasmaAB142P-M26,MEPlasmaAB142P-M28,MEPlasmaAB142P-M31,MEPlasmaAB142P-M34,MEPlasmaAB142P-M36,MEPlasmaGFAP-M12,MEPlasmaGFAP-M4,MEPlasmaGFAP-M9,MEPlasmaNfL-M21,MEPlasmaNfL-M25,MEPlasmaNfL-M26,MEPlasmaNfL-M3,MEPlasmaNfL-M32,MEPlasmaNfL-M33,MEPlasmaPTau181-M13,MEPlasmaPTau181-M18,MEPlasmaPTau181-M20,MEPlasmaPTau181-M23,MEPlasmaPTau181-M24,MEPlasmaPTau181-M26,MEPlasmaPTau181-M27,MEPlasmaPTau181-M30,MEPlasmaPTau181-M35,MEPlasmaPTau181-M37,MEPlasmaPTau181-M40,MEPlasmaPTau181-M43,MEPlasmaPTau181-M46,MEPlasmaPTau181-M47,MEPlasmaPTau181-M5,MEPlasmaPTau181-M8,MEPlasmapTau217-M3" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "STRING.INT.eigengenes" \
  # -k 0

# Rscript 02-Clustering.R  \
  # -i "MEMerged-M1,MEMerged-M10,MEMerged-M11,MEMerged-M12,MEMerged-M13,MEMerged-M14,MEMerged-M15,MEMerged-M16,MEMerged-M17,MEMerged-M18,MEMerged-M19,MEMerged-M2,MEMerged-M20,MEMerged-M21,MEMerged-M22,MEMerged-M24,MEMerged-M25,MEMerged-M26,MEMerged-M27,MEMerged-M28,MEMerged-M29,MEMerged-M3,MEMerged-M30,MEMerged-M32,MEMerged-M33,MEMerged-M35,MEMerged-M36,MEMerged-M4,MEMerged-M5,MEMerged-M6,MEMerged-M7,MEMerged-M8,MEMerged-M9,MEPlasmaAB140P-M15,MEPlasmaAB140P-M17,MEPlasmaAB140P-M18,MEPlasmaAB140P-M22,MEPlasmaAB140P-M24,MEPlasmaAB140P-M28,MEPlasmaAB140P-M38,MEPlasmaAB140P-M40,MEPlasmaAB140P-M41,MEPlasmaAB140P-M44,MEPlasmaAB140P-M45,MEPlasmaAB140P-M48,MEPlasmaAB140P-M50,MEPlasmaAB142P-M25,MEPlasmaAB142P-M26,MEPlasmaAB142P-M28,MEPlasmaAB142P-M31,MEPlasmaAB142P-M34,MEPlasmaAB142P-M36,MEPlasmaGFAP-M12,MEPlasmaGFAP-M4,MEPlasmaGFAP-M9,MEPlasmaNfL-M21,MEPlasmaNfL-M25,MEPlasmaNfL-M26,MEPlasmaNfL-M3,MEPlasmaNfL-M32,MEPlasmaNfL-M33,MEPlasmaPTau181-M13,MEPlasmaPTau181-M18,MEPlasmaPTau181-M20,MEPlasmaPTau181-M23,MEPlasmaPTau181-M24,MEPlasmaPTau181-M26,MEPlasmaPTau181-M27,MEPlasmaPTau181-M30,MEPlasmaPTau181-M35,MEPlasmaPTau181-M37,MEPlasmaPTau181-M40,MEPlasmaPTau181-M43,MEPlasmaPTau181-M46,MEPlasmaPTau181-M47,MEPlasmaPTau181-M5,MEPlasmaPTau181-M8,MEPlasmapTau217-M3" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "STRING.INT.eigengenes" \
  # -k 6 
  # # --final="6,1,2,4,3,5"


# Rscript 02-Clustering.R  \
  # -i "APP Metabolism,Apoptosis,Autophagy,Cell Cycle,Endolysosome,Immune Response,Lipid Metabolism,Metal Binding and Homeostasis,Mitochondrial Metabolism,Myelination,Oxidative Stress,Proteostasis,Structural Stabilization,Synapse,Tau Homeostasis,Vasculature" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "treatad_alpha_scaled" \
  # -k 7 # tie for best between 2 & 7 

# Rscript 02-Clustering.R  \
  # -i "APP Metabolism,Apoptosis,Autophagy,Cell Cycle,Endolysosome,Immune Response,Lipid Metabolism,Metal Binding and Homeostasis,Mitochondrial Metabolism,Myelination,Oxidative Stress,Proteostasis,Structural Stabilization,Synapse,Tau Homeostasis,Vasculature" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "treatad_harmonic_spec" \
  # -k 7 # tie for best between 2 & 7 

# Rscript 02-Clustering.R  \
  # -i "APP Metabolism,Apoptosis,Autophagy,Cell Cycle,Endolysosome,Immune Response,Lipid Metabolism,Metal Binding and Homeostasis,Mitochondrial Metabolism,Myelination,Oxidative Stress,Proteostasis,Structural Stabilization,Synapse,Tau Homeostasis,Vasculature" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "treatad_alpha_spec" \
  # -k 12 # second best after 2

# Rscript 02-Clustering.R  \
  # -i "Astrocyte,Microglia,Neuron,Oligodendrocyte|Astrocyte,CGE,L2/3 IT,L4/5/6 IT-ET Car3,MGE,Microglia-PVM,NP/CT/L6b,Oligo,Vascular|APP Metabolism,Apoptosis,Autophagy,Cell Cycle,Endolysosome,Immune Response,Lipid Metabolism,Metal Binding and Homeostasis,Mitochondrial Metabolism,Myelination,Oxidative Stress,Proteostasis,Structural Stabilization,Synapse,Tau Homeostasis,Vasculature" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "Sharma.lvl1,SEA-AD.Neighborhood.1kDS.lvl2,treatad_alpha_scaled" \
  # -k 9 # second best


# # Trying with ranksum --------------------------------------
# scores=('ranksum_scaled' 'ranksum_spec')
# for s in ${scores[*]}; do
  # Rscript 01-Long_plots.R  \
    # -i "APP Metabolism,Apoptosis,Autophagy,Cell Cycle,DNA Repair,Endolysosome,Epigenetic,Immune Response,Lipid Metabolism,Metal Binding and Homeostasis,Mitochondrial Metabolism,Myelination,Oxidative Stress,Proteostasis,Structural Stabilization,Synapse,Tau Homeostasis,Vasculature" \
    # -n "dataProt_SS-205063.baseline.all_bios" \
    # -s "treatad_${s}"
# done

# for s in ${scores[*]}; do
  # Rscript 02-Clustering.R  \
    # -i "APP Metabolism,Apoptosis,Autophagy,Cell Cycle,DNA Repair,Endolysosome,Epigenetic,Immune Response,Lipid Metabolism,Metal Binding and Homeostasis,Mitochondrial Metabolism,Myelination,Oxidative Stress,Proteostasis,Structural Stabilization,Synapse,Tau Homeostasis,Vasculature" \
    # -n "dataProt_SS-205063.baseline.all_bios" \
    # -s "treatad_${s}" \
    # -k 0
# done

# Rscript 02-Clustering.R  \
  # -i "APP Metabolism,Apoptosis,Autophagy,Cell Cycle,DNA Repair,Endolysosome,Epigenetic,Immune Response,Lipid Metabolism,Metal Binding and Homeostasis,Mitochondrial Metabolism,Myelination,Oxidative Stress,Proteostasis,Structural Stabilization,Synapse,Tau Homeostasis,Vasculature" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "treatad_ranksum_scaled" \
  # -k 5 # tied for first with 2

# Rscript 02-Clustering.R  \
  # -i "APP Metabolism,Apoptosis,Autophagy,Cell Cycle,DNA Repair,Endolysosome,Epigenetic,Immune Response,Lipid Metabolism,Metal Binding and Homeostasis,Mitochondrial Metabolism,Myelination,Oxidative Stress,Proteostasis,Structural Stabilization,Synapse,Tau Homeostasis,Vasculature" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "treatad_ranksum_spec" \
  # -k 10 # third best after 2 & 4

# Rscript 02-Clustering.R  \
  # -i "Astrocyte,Microglia,Neuron,Oligodendrocyte|Astrocyte,CGE,L2/3 IT,L4/5/6 IT-ET Car3,MGE,Microglia-PVM,NP/CT/L6b,Oligo,Vascular|APP Metabolism,Apoptosis,Autophagy,Cell Cycle,DNA Repair,Endolysosome,Epigenetic,Immune Response,Lipid Metabolism,Metal Binding and Homeostasis,Mitochondrial Metabolism,Myelination,Oxidative Stress,Proteostasis,Structural Stabilization,Synapse,Tau Homeostasis,Vasculature" \
  # -n "dataProt_SS-205063.baseline.all_bios" \
  # -s "Sharma.lvl1,SEA-AD.Neighborhood.1kDS.lvl2,treatad_ranksum_scaled" \
  # -k 6 # tied for third with 14 after 2 & 3

# have not done clustering with WSS, skip


# WSS heatmap ------------------------------
Rscript 03-Combined-TREAT-AD_GSEA.R  \
  -i "Merged.M1,Merged.M11,Merged.M7,Merged.M9,Merged.M2,PlasmaPTau181.M8,Merged.M6,PlasmaAB140P.M17,Merged.M12,Merged.M4,Merged.M5,PlasmaAB140P.M20,Merged.M8,PlasmaPTau181.M35,Merged.M3,Merged.M26,Merged.M13,Merged.M10,PlasmaPTau181.M18,Merged.M14,Merged.M19,Merged.M16,PlasmaPTau181.M26,Merged.M29,Merged.M24,PlasmaPTau181.M13,PlasmaPTau181.M47,Merged.M25,PlasmaAB140P.M35,PlasmaPTau181.M46,Merged.M20,Merged.M32,Merged.M22,Merged.M15,Merged.M21,Merged.M18,PlasmaPTau181.M5,PlasmaAB140P.M15,Merged.M27,PlasmaAB140P.M40,Merged.M17,PlasmaNfL.M3,PlasmaAB140P.M18,PlasmaPTau181.M30,PlasmaGFAP.M12" \
  -n "dataProt_SS-205063.baseline_WSS" \
  -s "all_WSS"

