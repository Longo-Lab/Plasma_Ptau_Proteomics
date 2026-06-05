#!/usr/bin/env bash

# # # # sbatch Rscript on the fly with --wrap
# # # sbatch -J ptau --mem=12G -c 16 -t 06:00:00 -p interactive \
  # # # -o %A_02-linreg_lme4.log \
  # # # --wrap "ml R/4.2.2; Rscript 02-linreg_lme4.R"
# # # run local faster

# # Rscript 06-Milind_etal.R
# # Rscript 07-Mostafavi_etal.R

names=(
  "PlasmaPTau181"
  "PlasmaAB140P" 
  "PlasmaAB142P" 
  "Plasma_NFL" 
  "Plasma_GFAP" 
  "Plasma_pTau217" 
  "PlasmaABRatio" 
  "Plasma_pTau217_AB42Ratio"
)
recodes=(
  'PlateId,P0028462,P0028452,P0028449,P0028428,P0028565,P0028463,P0028437'
  'PlateId,P0028451,P0028463,P0028462,P0028438,P0028430,P0028437'
  'PlateId,P0028438,P0028439,P0028452,P0028449,P0028565,P0028451,P0028447,P0028430,P0028437'
  'PlateId,P0028430,P0028439,P0028463,P0028438,P0028428,P0028447,P0028462,P0028448'
  'PlateId,P0028448,P0028438,P0028565,P0028437'
  'PlateId,P0028428,P0028439,P0028462,P0028437,P0028430,P0028448,P0028452,P0028447'
  'PlateId,P0028448,P0028451,P0028438,P0028565,P0028449,P0028452,P0028428,P0028447'
  'PlateId,P0028437,P0028438,P0028463,P0028430,P0028462,P0028439,P0028428,P0028447,P0028448'
)
covs=(
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
)
apoe=("" ".APOE")
modules=("TMT-AD" "TMT-Resilience" "TREAT-AD" "Wan_etal")

# # 01-sample_QC.R  ------------------------------------------
# for a in "${apoe[@]}"; do
#   for name in "${names[@]}"; do
#     if [[ "${name}${a}" == *"Ratio"* ]]; then
#       if [[ "${name}${a}" == *"APOE"* ]]; then
#         covs='Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage,APOE'
#       else  
#         covs='Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage'
#       fi
#     elif [[ "${name}${a}" == *"APOE"* ]]; then
#       covs='Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,APOE'
#     else
#       covs='Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage'
#     fi

#     Rscript ../baseline/01-sample_QC.R \
#     -i "../dataProt_SS-205063.merged.QC.csv" \
#     -o "${name}${a}" \
#     -p "../prot_names.txt" \
#     -v "$covs" \
#     -c "$name" \
#     -n "dataProt_SS-205063.PC.${name}${a}" \
#     -m "../genomic/mypcair.PCs.csv"
#   done
# done

# # 02-linreg_glm.R ------------------------------------------
# parallel -j 0 -k Rscript ../baseline/02-linreg_glm.R \
#   -i "{1}{4}/dataProt_SS-205063.PC.{1}{4}.INTnorm.csv.gz"   \
#   -o {1}{4}   \
#   -p "../prot_names.txt" \
#   -v {3}{=4 s/\./\,/ =}  \
#   -d {2}  \
#   -c "{1}"   \
#   -n "dataProt_SS-205063.PC.{1}{4}" \
#   ::: "${names[@]}" \
#   :::+ "${recodes[@]}"  \
#   :::+ "${covs[@]}"  \
#   ::: "${apoe[@]}"

# for a in "${apoe[@]}"; do
#   if [[ "$a" != "" ]]; then
#     covs='Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,APOE,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
#   else
#     covs='Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
#   fi
#   Rscript ../baseline/02-linreg_glm.R \
#     -i "PlasmaPTau181$a/dataProt_SS-205063.PC.PlasmaPTau181$a.INTnorm.csv.gz"   \
#     -o "PlasmaPTau181$a"   \
#     -p "../prot_names.txt" \
#     -v "$covs"   \
#     -d 'PlateId,P0028462,P0028452,P0028449,P0028428,P0028565,P0028463,P0028437'  \
#     -c "PlasmaPTau181"   \
#     -n "dataProt_SS-205063.PC.PlasmaPTau181$a"

#   Rscript ../baseline/02-linreg_glm.R \
#     -i "PlasmaAB140P$a/dataProt_SS-205063.PC.PlasmaAB140P$a.INTnorm.csv.gz"   \
#     -o "PlasmaAB140P$a"   \
#     -p "../prot_names.txt" \
#     -v "$covs"   \
#     -d 'PlateId,P0028451,P0028463,P0028462,P0028438,P0028430,P0028437'  \
#     -c "PlasmaAB140P"   \
#     -n "dataProt_SS-205063.PC.PlasmaAB140P$a"

#   Rscript ../baseline/02-linreg_glm.R \
#     -i "PlasmaAB142P$a/dataProt_SS-205063.PC.PlasmaAB142P$a.INTnorm.csv.gz"   \
#     -o "PlasmaAB142P$a"   \
#     -p "../prot_names.txt" \
#     -v "$covs"   \
#     -d 'PlateId,P0028438,P0028439,P0028452,P0028449,P0028565,P0028451,P0028447,P0028430,P0028437'  \
#     -c "PlasmaAB142P"   \
#     -n "dataProt_SS-205063.PC.PlasmaAB142P$a"

#   Rscript ../baseline/02-linreg_glm.R \
#     -i "Plasma_NFL$a/dataProt_SS-205063.PC.Plasma_NFL$a.INTnorm.csv.gz"   \
#     -o "Plasma_NFL$a"   \
#     -p "../prot_names.txt" \
#     -v "$covs"   \
#     -d 'PlateId,P0028430,P0028439,P0028463,P0028438,P0028428,P0028447,P0028462,P0028448'  \
#     -c "Plasma_NFL"   \
#     -n "dataProt_SS-205063.PC.Plasma_NFL$a"

#   Rscript ../baseline/02-linreg_glm.R \
#     -i "Plasma_GFAP$a/dataProt_SS-205063.PC.Plasma_GFAP$a.INTnorm.csv.gz"   \
#     -o "Plasma_GFAP$a"   \
#     -p "../prot_names.txt" \
#     -v "$covs"   \
#     -d 'PlateId,P0028448,P0028438,P0028565,P0028437'  \
#     -c "Plasma_GFAP"   \
#     -n "dataProt_SS-205063.PC.Plasma_GFAP$a"

#   Rscript ../baseline/02-linreg_glm.R \
#     -i "Plasma_pTau217$a/dataProt_SS-205063.PC.Plasma_pTau217$a.INTnorm.csv.gz"   \
#     -o "Plasma_pTau217$a"   \
#     -p "../prot_names.txt" \
#     -v "$covs"   \
#     -d 'PlateId,P0028428,P0028439,P0028462,P0028437,P0028430,P0028448,P0028452,P0028447'  \
#     -c "Plasma_pTau217"   \
#     -n "dataProt_SS-205063.PC.Plasma_pTau217$a"
# done

# # 02-linreg_glm.R for Ratios -----------------------------------
# for a in "${apoe[@]}"; do
#   if [[ "$a" != "" ]]; then
#     covs='Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage,APOE,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
#   else
#     covs='Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
#   fi
#   Rscript ../baseline/02-linreg_glm.R \
#     -i "PlasmaABRatio$a/dataProt_SS-205063.PC.PlasmaABRatio$a.INTnorm.csv.gz"   \
#     -o "PlasmaABRatio$a"   \
#     -p "../prot_names.txt" \
#     -v "$covs"   \
#     -d 'PlateId,P0028448,P0028451,P0028438,P0028565,P0028449,P0028452,P0028428,P0028447'  \
#     -c "PlasmaABRatio"   \
#     -n "dataProt_SS-205063.PC.PlasmaABRatio$a"

#   Rscript ../baseline/02-linreg_glm.R \
#     -i "Plasma_pTau217_AB42Ratio$a/dataProt_SS-205063.PC.Plasma_pTau217_AB42Ratio$a.INTnorm.csv.gz"   \
#     -o "Plasma_pTau217_AB42Ratio$a"   \
#     -p "../prot_names.txt" \
#     -v "$covs"   \
#     -d 'PlateId,P0028437,P0028438,P0028463,P0028430,P0028462,P0028439,P0028428,P0028447,P0028448'  \
#     -c "Plasma_pTau217_AB42Ratio"   \
#     -n "dataProt_SS-205063.PC.Plasma_pTau217_AB42Ratio$a"
# done

# # 03-Protein_QC.R ------------------------------------------
# for a in "${apoe[@]}"; do
#   for name in "${names[@]}"; do
#     Rscript ../baseline/03-Protein_QC.R \
#     -i "${name}${a}/dataProt_SS-205063.PC.${name}${a}.linreg_glm.csv.gz"  \
#     -o "${name}${a}" \
#     -r "../2024-08-30_Somalogic_reference_Key.txt" \
#     -n "dataProt_SS-205063.PC.${name}${a}"
#   done
# done

# # Rscript 03-Protein_QC.R \
# #   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.linreg_glm.csv.gz"  \
# #   -o "PlasmaAB140P" \
# #   -r "../Somalogic_reference_Key.txt" \
# #   -n "dataProt_SS-205063.PlasmaAB140P"

# # 04-TREAT-AD_correlation.R --------------------------------
# for a in "${apoe[@]}"; do
#   for name in "${names[@]}"; do
#     Rscript ../baseline/04-TREAT-AD_correlation.R \
#     -i "${name}${a}/dataProt_SS-205063.PC.${name}${a}.csv"  \
#     -o "${name}${a}" \
#     -n "dataProt_SS-205063.PC.${name}${a}" \
#     -c "$name"
#   done
# done

# # Rscript 04-TREAT-AD_correlation.R \
# #   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
# #   -o "PlasmaAB140P" \
# #   -n "dataProt_SS-205063.PlasmaAB140P" \
# #   -c "PlasmaAB140P"

# # 05-Combined-Enrich.R -------------------------------------
# for a in "${apoe[@]}"; do
#   for name in "${names[@]}"; do
#     for mod in "${modules[@]}"; do
#     Rscript ../baseline/05-Combined-Enrich.R \
#     -i "${name}${a}/dataProt_SS-205063.PC.${name}${a}.csv"  \
#     -o "${name}${a}" \
#     -n "dataProt_SS-205063.PC.${name}${a}" \
#     -m "$mod"
#     done
#   done
# done

# Rscript 05-Combined-Enrich.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
#   -o "PlasmaAB140P" \
#   -n "dataProt_SS-205063.PlasmaAB140P" \
#   -m "TMT-AD"

# 06-TREAT-AD_GSEA.R ---------------------------------------
for a in "${apoe[@]}"; do
  for name in "${names[@]}"; do
    Rscript ../baseline/06-TREAT-AD_GSEA.R \
    -i "${name}${a}/dataProt_SS-205063.PC.${name}${a}.csv"  \
    -o "${name}${a}" \
    -n "dataProt_SS-205063.PC.${name}${a}"
  done
done

# Rscript ../baseline/06-TREAT-AD_GSEA.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PC.PlasmaAB140P.csv"  \
#   -o "PlasmaAB140P" \
#   -n "dataProt_SS-205063.PC.PlasmaAB140P" 

# 07-Sharma_etal_correlation.R -----------------------------
for a in "${apoe[@]}"; do
  for name in "${names[@]}"; do
    Rscript ../baseline/07-Sharma_etal_correlation.R \
    -i "${name}${a}/dataProt_SS-205063.PC.${name}${a}.csv"  \
    -o "${name}${a}" \
    -n "dataProt_SS-205063.PC.${name}${a}" \
    -c "$name" 
  done
done

# Rscript 07-Sharma_etal_correlation.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
#   -o "PlasmaAB140P" \
#   -n "dataProt_SS-205063.PlasmaAB140P" \
#   -c "PlasmaAB140P" 
