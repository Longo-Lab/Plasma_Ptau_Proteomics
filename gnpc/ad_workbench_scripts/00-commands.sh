#!/usr/bin/env bash

###################################################################
#Script Name  : 00-commands.sh 
#Description  : Main association analysis wrapper script
#Usage        : bash 00-commands.sh
#Author       : Robert R Butler III
#Date Created : 2025-10-28
#Email        : rrbutler@stanford.edu
#Copyright (c) 2025
###################################################################

# Notes       :
# must have run scripts in the main folder & associated specificity files
# cell_type_prot - ewce cell-type specificity matrices for Sharma et al & SEA-AD
# treat-ad - biodomain specificity matrix
# string - STRING module specificity matrix
 

# 01-Dataset_Exploration.Rmd
# run interactively
grps=("A"   "B"   "C"   "D"   "E"   "F"   "I"   "L"   "N"   "P"   "Q"   "Mix")

# 02-GNPC_imputation.R -----------------------------
# GNPC_K.csv all -1s are for IRF6, which appears to be failed probe. Dropped. 
# In the future could add to Mix and use that to impute IRF6.
for k in ${grps[@]}; do
  Rscript 02-GNPC_imputation.R \
    -i "GNPC_${k}.csv" \
    -n "GNPC_${k}.imputed"
done

# merge to GNPC_table.imputed.csv
# Extract GENE column from first file, then append all other columns
paste GNPC_A.imputed.csv \
      <(cut -d',' -f2- GNPC_B.imputed.csv) \
      <(cut -d',' -f2- GNPC_C.imputed.csv) \
      <(cut -d',' -f2- GNPC_D.imputed.csv) \
      <(cut -d',' -f2- GNPC_E.imputed.csv) \
      <(cut -d',' -f2- GNPC_F.imputed.csv) \
      <(cut -d',' -f2- GNPC_I.imputed.csv) \
      <(cut -d',' -f2- GNPC_L.imputed.csv) \
      <(cut -d',' -f2- GNPC_N.imputed.csv) \
      <(cut -d',' -f2- GNPC_P.imputed.csv) \
      <(cut -d',' -f2- GNPC_Q.imputed.csv) \
      <(cut -d',' -f2- GNPC_Mix.imputed.csv) \
      -d',' > GNPC.imputed.csv


# 03-WSS_construction.R ------------------------------
for k in ${grps[@]}; do
  Rscript 03-WSS_construction.R \
    -i "cell_type_prot/STRING_merged_scaled_scores.csv" \
    -n "GNPC_${k}.STRING" \
    -r "GNPC_${k}.imputed.csv" \
    --INT
done

# merge to GNPC_WSS.csv
awk 'FNR==1 && NR!=1{next;}{print}' celltype_scores/GNPC_*.STRING.INT.WSS.csv \
  > celltype_scores/GNPC.STRING.INT.WSS.csv


# WSS with Log2
for k in ${grps[@]}; do
  Rscript 03-WSS_construction.R \
    -i "cell_type_prot/STRING_merged_scaled_scores.csv" \
    -n "GNPC_${k}.STRING" \
    -r "GNPC_${k}.imputed.csv"
done

# merge to GNPC_WSS.csv
awk 'FNR==1 && NR!=1{next;}{print}' celltype_scores/GNPC_*.STRING.WSS.csv \
  > celltype_scores/GNPC.STRING.WSS.csv

# 04-WSS_validation.R -------------------------------
for k in ${grps[@]}; do
  Rscript 04-WSS_validation.R \
    -n "ADRCvsGNPC.$k" \
    -t 2 \
    -v "GNPC_$k"
done

# 05-Validation_stats.Rmd ---------------------------
# Run interactively


# 06-WSS_x_clin -------------------------------------
# Not much in the way of covariate selection at the GNPC level. 
# We can continue to do svs, and look at their & race's  respective importance 

# binary variables ------------
names=(
  "tia"
  "tbi"
  "cancer"
  "diabetes"
  "chf"
  "copd"
  "mi"
  "afib"
  "angina"
  "hyperlipidaemia"
  "hypertension"
  "depression"
  "anxiety"
  "diagnosis_bin"
  "last_diag_bin"
)

# QC and surrogates
parallel -j 1 Rscript 06-sample_QC_WSS_clin.R \
  -o "WSS_clin/{}"   \
  -v "age_at_visit,sex,contributor_code,race" \
  -c {} \
  -n "GNPC.baseline_clin_binary.{}" \
  --binomial \
  ::: "${names[@]}" 

names=(
  # "tia"
  # "tbi"
  # "cancer"
  "diabetes"
  "chf"
  "copd"
  # "mi"
  "afib"
  # "angina"
  # "hyperlipidaemia"
  # "hypertension"
  "depression"
  # "anxiety"
  "diagnosis_bin"
  # "last_diag_bin"
)

# QC and surrogates (rerun fails with no sv calculation)
parallel -j 1 Rscript 06-sample_QC_WSS_clin.R \
  -o "WSS_clin/{}"   \
  -v "age_at_visit,sex,contributor_code,race" \
  -c {} \
  -n "GNPC.baseline_clin_binary.{}" \
  --binomial \
  --no-svs \
  ::: "${names[@]}" 

names=(
  "tia"
  "tbi"
  "cancer"
  "diabetes"
  "chf"
  "copd"
  "mi"
  "afib"
  "angina"
  "hyperlipidaemia"
  "hypertension"
  "depression"
  "anxiety"
  "diagnosis_bin"
  "last_diag_bin"
)

# 08-linreg_lme4qtl.R
parallel -j 1 Rscript 08-linreg_glm.R \
  -i "WSS_clin/{}/GNPC.baseline_clin_binary.{}.WSSnorm.csv.gz" \
  -o "WSS_clin/{}"   \
  -p "WSS_names.txt" \
  -v "age_at_visit,sex,contributor_code,race" \
  -s "age_at_visit" \
  -c "{}" \
  -n "GNPC.baseline_clin_binary.{}" \
  --binomial \
  -f "0,1" \
  ::: "${names[@]}"


# diagnosis_bin as +/-1
names=(
  "diagnosis_bin"
)

# 08-linreg_lme4qtl.R
parallel -j 1 Rscript 08-linreg_glm.R \
  -i "WSS_clin/{}/GNPC.baseline_clin_binary.{}.WSSnorm.csv.gz" \
  -o "WSS_clin/{}"   \
  -p "WSS_names.txt" \
  -v "age_at_visit,sex,contributor_code,race" \
  -s "age_at_visit" \
  -c "{}" \
  -n "GNPC.baseline_clin_continuous.{}" \
  ::: "${names[@]}"



# continuous variables --------
names=(
  "MoCA"
  "last_MoCA"
  "MoCA_velocity"
  "MoCA_v_floored"
)

# QC and surrogates
parallel -j 1 Rscript 06-sample_QC_WSS_clin.R \
  -o "WSS_clin/{}"   \
  -v "age_at_visit,sex,contributor_code,race" \
  -c {} \
  -n "GNPC.baseline_clin_continuous.{}" \
  --no-svs \
  ::: "${names[@]}" 

# 08-linreg_glm.R
parallel -j 1 Rscript 08-linreg_glm.R \
  -i "WSS_clin/{}/GNPC.baseline_clin_continuous.{}.WSSnorm.csv.gz" \
  -o "WSS_clin/{}"   \
  -p "WSS_names.txt" \
  -v "age_at_visit,sex,contributor_code,race" \
  -s "age_at_visit" \
  -c "{}" \
  -n "GNPC.baseline_clin_continuous.{}" \
  ::: "${names[@]}" 


# 07-WSS_x_prot -------------------------------------
names=(
  # "Merged.M1" "Merged.M11" "Merged.M7" "Merged.M9" "Merged.M2" 
  "PlasmaPTau181.M8" "Merged.M6" "PlasmaAB140P.M17" "Merged.M12"
  "Merged.M4" "Merged.M5" "PlasmaAB140P.M20" "Merged.M8"
  "PlasmaPTau181.M35" "Merged.M3" "Merged.M26" "Merged.M13"
  "Merged.M10" "PlasmaPTau181.M18" "Merged.M14" "Merged.M19"
  "Merged.M16" "PlasmaPTau181.M26" "Merged.M29" "Merged.M24"
  "PlasmaPTau181.M13" "PlasmaPTau181.M47" "Merged.M25" "PlasmaAB140P.M35"
  "PlasmaPTau181.M46" "Merged.M20" "Merged.M32" "Merged.M22"
  "Merged.M15" "Merged.M21" "Merged.M18" "PlasmaPTau181.M5"
  "PlasmaAB140P.M15" "Merged.M27" "PlasmaAB140P.M40" "Merged.M17"
  "PlasmaNfL.M3" "PlasmaAB140P.M18" "PlasmaPTau181.M30" "PlasmaGFAP.M12"
)

# QC and surrogates (AD workbench is too much a shitshow to be able to compute svs)
parallel -j 1 Rscript 07-sample_QC_WSS_prot.R \
  -i "GNPC_baseline.csv" \
  -o "WSS_prot/{}"   \
  -v "age_at_visit,sex,contributor_code,race" \
  -c {} \
  --no-svs \
  -n "GNPC.baseline_prot.{}" \
  ::: "${names[@]}" 

# 08-linreg_glm.R
parallel -j 1 Rscript 08-linreg_glm.R \
  -i "WSS_prot/{}/GNPC.baseline_prot.{}.INTnorm.csv.gz" \
  -o "WSS_prot/{}"   \
  -p "prot_names.txt" \
  -v "age_at_visit,sex,contributor_code,race" \
  -s "age_at_visit" \
  -c "{}" \
  -n "GNPC.baseline_prot.{}" \
  ::: "${names[@]}" 


# AD workbench shutting down for two weeks(!?), need to copy linreg_results into a lifeboat
mkdir -p export
# zip -r linreg_results.zip WSS*/*/*linreg_glm.csv.gz
cp -r WSS*/*/*linreg_glm.csv.gz export/
# try with sphere
cp -r sphere/WSS_clin/*/synthetic_output* export/

