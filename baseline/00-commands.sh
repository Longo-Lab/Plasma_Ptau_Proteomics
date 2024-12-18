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

#     Rscript 01-sample_QC.R \
#     -i "../dataProt_SS-205063.merged.QC.csv" \
#     -o "${name}${a}" \
#     -p "../prot_names.txt" \
#     -v "$covs" \
#     -c "$name" \
#     -n "dataProt_SS-205063.${name}${a}"
#   done
# done

# Rscript 01-sample_QC.R \
#   -i "../dataProt_SS-205063.merged.QC.csv"   \
#   -o "PlasmaAB140P"   \
#   -p "../prot_names.txt" \
#   -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage'   \
#   -c "PlasmaAB140P"   \
#   -n "dataProt_SS-205063.PlasmaAB140P"

# 02-linreg_glm.R ------------------------------------------
for a in "${apoe[@]}"; do
  if [[ "$a" != "" ]]; then
    covs='Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,APOE,sv1,sv2,sv3,sv4,sv5,sv6'
  else
    covs='Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6'
  fi
  Rscript 02-linreg_glm.R \
    -i "PlasmaPTau181$a/dataProt_SS-205063.PlasmaPTau181$a.INTnorm.csv.gz"   \
    -o "PlasmaPTau181$a"   \
    -p "../prot_names.txt" \
    -v "$covs"   \
    -d 'PlateId,P0028452,P0028439,P0028437,P0028449,P0028465,P0028565,P0028466'  \
    -c "PlasmaPTau181"   \
    -n "dataProt_SS-205063.PlasmaPTau181$a"

  Rscript 02-linreg_glm.R \
    -i "PlasmaAB140P$a/dataProt_SS-205063.PlasmaAB140P$a.INTnorm.csv.gz"   \
    -o "PlasmaAB140P$a"   \
    -p "../prot_names.txt" \
    -v "$covs"   \
    -d 'PlateId,P0028451,P0028447,P0028439,P0028463,P0028462,P0028438,P0028430,P0028437'  \
    -c "PlasmaAB140P"   \
    -n "dataProt_SS-205063.PlasmaAB140P$a"

  Rscript 02-linreg_glm.R \
    -i "PlasmaAB142P$a/dataProt_SS-205063.PlasmaAB142P$a.INTnorm.csv.gz"   \
    -o "PlasmaAB142P$a"   \
    -p "../prot_names.txt" \
    -v "$covs"   \
    -d 'PlateId,P0028462,P0028449,P0028439,P0028452,P0028565,P0028451,P0028430,P0028447,P0028437'  \
    -c "PlasmaAB142P"   \
    -n "dataProt_SS-205063.PlasmaAB142P$a"

  Rscript 02-linreg_glm.R \
    -i "Plasma_NFL$a/dataProt_SS-205063.Plasma_NFL$a.INTnorm.csv.gz"   \
    -o "Plasma_NFL$a"   \
    -p "../prot_names.txt" \
    -v "$covs"   \
    -d 'PlateId,P0028565,P0028463,P0028438,P0028447,P0028448,P0028462'  \
    -c "Plasma_NFL"   \
    -n "dataProt_SS-205063.Plasma_NFL$a"

  Rscript 02-linreg_glm.R \
    -i "Plasma_GFAP$a/dataProt_SS-205063.Plasma_GFAP$a.INTnorm.csv.gz"   \
    -o "Plasma_GFAP$a"   \
    -p "../prot_names.txt" \
    -v "$covs"   \
    -d 'PlateId,P0028463,P0028430,P0028448,P0028438,P0028565,P0028437,P0028447'  \
    -c "Plasma_GFAP"   \
    -n "dataProt_SS-205063.Plasma_GFAP$a"

  Rscript 02-linreg_glm.R \
    -i "Plasma_pTau217$a/dataProt_SS-205063.Plasma_pTau217$a.INTnorm.csv.gz"   \
    -o "Plasma_pTau217$a"   \
    -p "../prot_names.txt" \
    -v "$covs"   \
    -d 'PlateId,P0028438,P0028430,P0028448,P0028462,P0028439,P0028452,P0028447'  \
    -c "Plasma_pTau217"   \
    -n "dataProt_SS-205063.Plasma_pTau217$a"
done

# 02-linreg_glm.R for Ratios -----------------------------------
for a in "${apoe[@]}"; do
  if [[ "$a" != "" ]]; then
    covs='Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage,APOE,sv1,sv2,sv3,sv4,sv5,sv6'
  else
    covs='Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6'
  fi
  Rscript 02-linreg_glm.R \
    -i "PlasmaABRatio$a/dataProt_SS-205063.PlasmaABRatio$a.INTnorm.csv.gz"   \
    -o "PlasmaABRatio$a"   \
    -p "../prot_names.txt" \
    -v "$covs"   \
    -d 'PlateId,P0028448,P0028438,P0028437,P0028451,P0028565,P0028449,P0028447,P0028428,P0028452'  \
    -c "PlasmaABRatio"   \
    -n "dataProt_SS-205063.PlasmaABRatio$a"

  Rscript 02-linreg_glm.R \
    -i "Plasma_pTau217_AB42Ratio$a/dataProt_SS-205063.Plasma_pTau217_AB42Ratio$a.INTnorm.csv.gz"   \
    -o "Plasma_pTau217_AB42Ratio$a"   \
    -p "../prot_names.txt" \
    -v "$covs"   \
    -d 'PlateId,P0028452,P0028463,P0028428,P0028430,P0028462,P0028438,P0028447,P0028439,P0028448'  \
    -c "Plasma_pTau217_AB42Ratio"   \
    -n "dataProt_SS-205063.Plasma_pTau217_AB42Ratio$a"
done

# 03-Protein_QC.R ------------------------------------------
for a in "${apoe[@]}"; do
  for name in "${names[@]}"; do
    Rscript 03-Protein_QC.R \
    -i "${name}${a}/dataProt_SS-205063.${name}${a}.linreg_glm.csv.gz"  \
    -o "${name}${a}" \
    -r "../2024-08-30_Somalogic_reference_Key.txt" \
    -n "dataProt_SS-205063.${name}${a}"
  done
done

# Rscript 03-Protein_QC.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.linreg_glm.csv.gz"  \
#   -o "PlasmaAB140P" \
#   -r "../Somalogic_reference_Key.txt" \
#   -n "dataProt_SS-205063.PlasmaAB140P"

# 04-TREAT-AD_correlation.R --------------------------------
for a in "${apoe[@]}"; do
  for name in "${names[@]}"; do
    Rscript 04-TREAT-AD_correlation.R \
    -i "${name}${a}/dataProt_SS-205063.${name}${a}.csv"  \
    -o "${name}${a}" \
    -n "dataProt_SS-205063.${name}${a}" \
    -c "$name"
  done
done

# Rscript 04-TREAT-AD_correlation.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
#   -o "PlasmaAB140P" \
#   -n "dataProt_SS-205063.PlasmaAB140P" \
#   -c "PlasmaAB140P"

# 05-Combined-Enrich.R -------------------------------------
for a in "${apoe[@]}"; do
  for name in "${names[@]}"; do
    for mod in "${modules[@]}"; do
    Rscript 05-Combined-Enrich.R \
    -i "${name}${a}/dataProt_SS-205063.${name}${a}.csv"  \
    -o "${name}${a}" \
    -n "dataProt_SS-205063.${name}${a}" \
    -m "$mod"
    done
  done
done

# Rscript 05-Combined-Enrich.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
#   -o "PlasmaAB140P" \
#   -n "dataProt_SS-205063.PlasmaAB140P" \
#   -m "TMT-AD"

# 06-TREAT-AD_GSEA.R ---------------------------------------
for a in "${apoe[@]}"; do
  for name in "${names[@]}"; do
    Rscript 06-TREAT-AD_GSEA.R \
    -i "${name}${a}/dataProt_SS-205063.${name}${a}.csv"  \
    -o "${name}${a}" \
    -n "dataProt_SS-205063.${name}${a}"
  done
done

# Rscript 06-TREAT-AD_GSEA.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
#   -o "PlasmaAB140P" \
#   -n "dataProt_SS-205063.PlasmaAB140P" 

# 07-Sharma_etal_correlation.R -----------------------------
for a in "${apoe[@]}"; do
  for name in "${names[@]}"; do
    Rscript 07-Sharma_etal_correlation.R \
    -i "${name}${a}/dataProt_SS-205063.${name}${a}.csv"  \
    -o "${name}${a}" \
    -n "dataProt_SS-205063.${name}${a}" \
    -c "$name" 
  done
done

# Rscript 07-Sharma_etal_correlation.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
#   -o "PlasmaAB140P" \
#   -n "dataProt_SS-205063.PlasmaAB140P" \
#   -c "PlasmaAB140P" 
