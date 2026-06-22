#!/usr/bin/env bash

###################################################################
#Script Name  : 00-commands.sh 
#Description  : Backup GNPC post analysis wrapper script
#Usage        : bash 00-commands.sh
#Author       : Robert R Butler III
#Date Created : 2026-06-10
#Email        : rrbutler@stanford.edu
#Copyright (c) 2026
###################################################################

# Notes       :
# Setup for the AD workbench shutdown for figures


###################### GNPC backup analysis ############################
# WSS_x_prot ------------------------------
# should just run with the tweaked association scripts
names=(
  "Merged.M1" 
  "Merged.M11" 
  "Merged.M7" 
  "Merged.M9" 
  "Merged.M2" 
  "PlasmaPTau181.M8" 
  "Merged.M6" 
  "PlasmaAB140P.M17" 
  "Merged.M12" 
  "Merged.M4" 
  "Merged.M5" 
  "PlasmaAB140P.M20" 
  "Merged.M8" 
  "PlasmaPTau181.M35" 
  "Merged.M3" 
  "Merged.M26" 
  "Merged.M13" 
  "Merged.M10" 
  "PlasmaPTau181.M18" 
  "Merged.M14" 
  "Merged.M19" 
  "Merged.M16" 
  "PlasmaPTau181.M26" 
  "Merged.M29" 
  "Merged.M24" 
  "PlasmaPTau181.M13" 
  "PlasmaPTau181.M47" 
  "Merged.M25" 
  "PlasmaAB140P.M35" 
  "PlasmaPTau181.M46" 
  "Merged.M20" 
  "Merged.M32" 
  "Merged.M22" 
  "Merged.M15" 
  "Merged.M21" 
  "Merged.M18" 
  "PlasmaPTau181.M5" 
  "PlasmaAB140P.M15" 
  "Merged.M27" 
  "PlasmaAB140P.M40" 
  "Merged.M17" 
  "PlasmaNfL.M3" 
  "PlasmaAB140P.M18" 
  "PlasmaPTau181.M30" 
  "PlasmaGFAP.M12"
)



# 03-Protein_QC.R ------------------------------------------
# Rscript 03-Protein_QC.R \
  # -i "WSS_prot/Merged.M3/dataProt_SS-205063.baseline_WSS.Merged.M3.lme4qtl.csv.gz"  \
  # -o "WSS_prot/Merged.M3" \
  # -r "2024-08-30_Somalogic_reference_Key.txt" \
  # -n "dataProt_SS-205063.baseline_WSS.Merged.M3"
  # -p "cohens_d"

parallel -j 18 Rscript 03-Protein_QC.R \
  -i "WSS_prot/{}/GNPC.baseline_prot.{}.linreg_glm.csv.gz" \
  -o "WSS_prot/{}" \
  -n "GNPC.baseline_prot.{}" \
  ::: "${names[@]}" 

# 06-TREAT-AD_GSEA.R ---------------------------------------
# be sure to have a functioning synapse client
# Rscript 06-TREAT-AD_GSEA.R \
  # -i "WSS_prot/Merged.M3/dataProt_SS-205063.baseline_WSS.Merged.M3.csv"  \
  # -o "WSS_prot/Merged.M3" \
  # -n "dataProt_SS-205063.baseline_WSS.Merged.M3" \
  # -p "cohens_d"

# conda activate synpase
parallel -j 5 Rscript 06-TREAT-AD_GSEA.R   \
  -i "WSS_prot/{}/GNPC.baseline_prot.{}.csv"  \
  -o "WSS_prot/{}" \
  -n "GNPC.baseline_prot.{}" \
  -k "SYMBOL" \
  ::: "${names[@]}" 


mapfile -t missing < <(for dir in WSS_prot/*/; do
  if ! ls "$dir"*.GSEA_Biodomain.rds 1>/dev/null 2>&1; then
    basename "$dir"
  fi
done)

parallel -j 5 Rscript 06-TREAT-AD_GSEA.R   \
  -i "WSS_prot/{}/GNPC.baseline_prot.{}.csv"  \
  -o "WSS_prot/{}" \
  -n "GNPC.baseline_prot.{}" \
  -k "SYMBOL" \
  ::: "${missing[@]}" 



# Binary variables -----------------------
names=(
  "HYPERCHO"
  "HYPERTEN"
  "DEP2YRS"
  "ANXIETY"
  "DIABETES"
  "Diagnosis_bin"
  "last_diag_bin"
)


# Continuous variables -----------------------
names=(
  "MoCA"
  "CDRSUM"
  "last_MoCA"
  "last_CDRSUM"
  "MoCA_velocity"
  "MoCA_v_floored"
)






