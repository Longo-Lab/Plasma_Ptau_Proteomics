#!/usr/bin/env bash
###################################################################
#Script Name  : 00-commands.sh
#Description  : Replication cohort wrapper for Figure 6 (WSS module
#               prediction) association analyses
#Usage        : bash 00-commands.sh
#Author       : Robert R Butler III
#Date Created : 2026-08-10
#Email        : rrbutler@stanford.edu
#Copyright (c) 2026
###################################################################
# Notes:
# Adapted from lines 610-880 of the main project's 00-commands.sh for an
# independent, cross-sectional external collaborator cohort. Before
# running this script:
#   1. Run 01-WSSmod_prep.Rmd -> produces dataProt_replication.csv.gz
#   2. Run 02-sample_QC_prediction.R, inspect the resulting
#      *.feature_importance.pdf (and beta/SE distributions if you used
#      --svs), and edit the covs=(...) placeholders below to
#      match what you decide is worth controlling for. See README.md
#      ("Choosing covariates" and "Surrogate variables") before you do.
#   3. Then run the remainder of the scripts (03 & 04).
# Differences from the internal "follow-up cohort" pipeline:
#   - Uses 03-linreg_glm.R (plain glm) instead of 02-linreg_lme4qtl.R,
#     since a cross-sectional cohort has no repeated measures to model
#   - No internal batch covariates (Storage_days/Study/Subarray/
#     Fed.Fasted/HybControlNormScale) or genomic PCs/SVA baked in --
#     those are chosen per-cohort via 02-sample_QC_prediction.R instead
#   - Skips 06-TREAT-AD_GSEA.R (BioDomain enrichment) entirely
# Requires GNU parallel (or run each Rscript call individually).
# run with R v4.4.1
###################################################################

INPUT="dataProt_replication.csv.gz"
REF_ANNOTATION="2024-08-30_Somalogic_reference_Key.txt"

###################### Replication cohort prot predictions #######################
# module names to be tested (all R2 > 0.25 in the original model fit)
names=(
  "Merged.M13"
  "Merged.M8"
  "Merged.M27"
  "Merged.M5"
  "Merged.M2"
  "Merged.M29"
  "Merged.M1"
  "Merged.M32"
  "Merged.M15"
  "Merged.M18"
  "Merged.M21"
  "Merged.M14"
)

# actual or predicted (.pred) module score as outcome
pred=("" ".pred")

# All modules use the same covariates here -- 
# set this for your dataset's known covariates to run feature importance estimation
COVS="Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,Storage,PC1,PC2..."

# 02-sample_QC_prediction.R (stop and inspect covariates before progressing) ---
# Rscript 02-sample_QC_prediction.R \
  # -i "dataProt_replication.csv.gz" \
  # -p "prot_names.txt" \
  # -v "Age,Gender" \
  # -c "Merged.M13" \
  # --svs 6 \
  # -n "dataProt_replication.Merged.M13" \
  # -o "WSS_predict/Merged.M13"

parallel -j 4 Rscript 02-sample_QC_prediction.R \
  -o "WSS_predict/{1}{2}"   \
  -v "$COVS" \
  -c {1}{2} \
  -n "dataProt_replication.{1}{2}" \
  ::: "${names[@]}" \
  ::: "${pred[@]}"


# 03-linreg_glm.R --------------------------------
# After above, you can include important covs & svs 
#   this can be an array of different cov strings for each names:pred combo (see below line 151)
#   or a fixed string to apply to all (as was the case in follow-up)
covs="Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,Storage"

# Rscript 03-linreg_glm.R \
  # -i "WSS_predict/Merged.M13/dataProt_replication.Merged.M13.INTnorm.csv.gz" \
  # -p "prot_names.txt" \
  # -v "Age,Gender" \
  # -c "Merged.M13" \
  # -n "dataProt_replication.Merged.M13" \
  # -o "WSS_predict/Merged.M13"

parallel -j 4 Rscript 03-linreg_glm.R \
  -i "WSS_predict/{1}{2}/dataProt_replication.{1}{2}.INTnorm.csv.gz" \
  -p "prot_names.txt" \
  -v "$covs" \
  -c "{1}{2}" \
  -n "dataProt_replication.{1}{2}" \
  -o "WSS_predict/{1}{2}" \
  ::: "${names[@]}" \
  ::: "${pred[@]}"


# 04-Protein_QC.R ------------------------------------------
# ix here is a SomaLogic probe ID, so -r maps probe -> gene symbol
# (needs internet access for the biomaRt Ensembl lookup)
# Rscript 04-Protein_QC.R \
  # -i "WSS_predict/Merged.M13/dataProt_replication.Merged.M13.linreg_glm.csv.gz" \
  # -o "WSS_predict/Merged.M13" \
  # -r "$REF_ANNOTATION" \
  # -n "dataProt_replication.Merged.M13"

parallel -j 4 Rscript 04-Protein_QC.R \
  -i "WSS_predict/{1}{2}/dataProt_replication.{1}{2}.linreg_glm.csv.gz" \
  -o "WSS_predict/{1}{2}" \
  -r "$REF_ANNOTATION" \
  -n "dataProt_replication.{1}{2}" \
  ::: "${names[@]}" \
  ::: "${pred[@]}"


###################### Replication cohort clin predictions #######################
# Continuous variable: Cognitive Score (MoCA or equivalent) -----------------------
cog="CognitiveScore"
# set this for your dataset's known covariates to run feature importance estimation
COVS="Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,Storage,PC1,PC2..."

# 02-sample_QC_prediction.R (stop and inspect covariates before progressing) ---
# Rscript 02-sample_QC_prediction.R \
  # -i "dataProt_replication.csv.gz" \
  # -p "WSS_names.txt" \
  # -v "Age,Gender" \
  # -c "CognitiveScore" \
  # --svs 6 \
  # --clin \
  # -n "dataProt_replication.CognitiveScore" \
  # -o "WSS_clin/CognitiveScore"

parallel -j 2 Rscript 02-sample_QC_prediction.R \
  -o "WSS_clin/${cog}{2}" \
  -p {1} \
  -v "$COVS" \
  -c "$cog" \
  --clin \
  -n "dataProt_replication.clin_continuous.${cog}{2}" \
  ::: "WSS_names.txt" "WSS_names_pred.txt" \
  :::+ "${pred[@]}"


# WSS_names.txt (actual) and WSS_names_pred.txt (predicted) may warrant
# different covariates/SVs -- covs is linked (:::+) to that pair, so set
# each from your own 02-sample_QC_prediction.R -p WSS_names[_pred].txt
# -c CognitiveScore ... runs
covs=(
  "Age,Gender,..."   # actual   (WSS_names.txt)
  "Age,Gender,..."   # predicted (WSS_names_pred.txt)
)

# 03-linreg_glm.R --------------------------------
parallel -j 2 Rscript 03-linreg_glm.R \
  -i "WSS_clin/${cog}{2}/dataProt_replication.clin_continuous.${cog}{2}.INTnorm.csv.gz" \
  -p "{1}" \
  -v "{3}" \
  -c "$cog" \
  -n "dataProt_replication.clin_continuous.${cog}{2}" \
  -o "WSS_clin/${cog}{2}" \
  ::: "WSS_names.txt" "WSS_names_pred.txt" \
  :::+ "" ".pred" \
  :::+ "${covs[@]}"

# Binary variable: AD Diagnosis -----------------------
diag="Diagnosis_bin" # assumes 0/1 encoding otherwise specify with -f
# set this for your dataset's known covariates to run feature importance estimation
COVS="Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,Storage,PC1,PC2..."

# 02-sample_QC_prediction.R (stop and inspect covariates before progressing) ---
# Rscript 02-sample_QC_prediction.R \
  # -i "dataProt_replication.csv.gz" \
  # -p "WSS_names.txt" \
  # -v "Age,Gender" \
  # -c "Diagnosis_bin" \
  # -b \
  # --svs 6 \
  # --clin \
  # -n "dataProt_replication.Diagnosis_bin" \
  # -o "WSS_clin/Diagnosis_bin"

parallel -j 2 Rscript 02-sample_QC_prediction.R \
  -o "WSS_clin/${diag}{2}" \
  -p {1} \
  -v "$COVS" \
  -c "$diag" \
  -b \
  --clin \
  -n "dataProt_replication.clin_binary.${diag}{2}" \
  ::: "WSS_names.txt" "WSS_names_pred.txt" \
  :::+ "${pred[@]}"

# WSS_names.txt (actual) and WSS_names_pred.txt (predicted) may warrant
# different covariates/SVs -- covs is linked (:::+) to that pair, so set
# each from your own 02-sample_QC_prediction.R -p WSS_names[_pred].txt
# -c Diagnosis_bin ... runs
covs=(
  "Age,Gender,..."   # actual   (WSS_names.txt)
  "Age,Gender,..."   # predicted (WSS_names_pred.txt)
)

# 03-linreg_glm.R --------------------------------
parallel -j 2 Rscript 03-linreg_glm.R \
  -i "WSS_clin/${diag}{2}/dataProt_replication.clin_binary.${diag}{2}.INTnorm.csv.gz" \
  -p "{1}" \
  -v "{3}" \
  -c "$diag" \
  -b \
  -f "0,1" \
  -n "dataProt_replication.clin_binary.${diag}{2}" \
  -o "WSS_clin/${diag}{2}" \
  ::: "WSS_names.txt" "WSS_names_pred.txt" \
  :::+ "" ".pred" \
  :::+ "${covs[@]}"


