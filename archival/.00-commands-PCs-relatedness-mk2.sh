#!/usr/bin/env bash

names=(
  # "PlasmaPTau181"
  "PlasmaAB140P" 
  # "PlasmaAB142P" 
  # "Plasma_NFL" 
  # "Plasma_GFAP" 
  # "Plasma_pTau217"
  # "PlasmaABRatio" 
  # "Plasma_pTau217_AB42Ratio" 
)
recodes=(
  # 'PlateId,P0028462,P0028452,P0028449,P0028428,P0028565,P0028463,P0028437'
  'PlateId,P0028451,P0028463,P0028462,P0028438,P0028430,P0028437'
  # 'PlateId,P0028438,P0028439,P0028452,P0028449,P0028565,P0028451,P0028447,P0028430,P0028437'
  # 'PlateId,P0028430,P0028439,P0028463,P0028438,P0028428,P0028447,P0028462,P0028448'
  # 'PlateId,P0028448,P0028438,P0028565,P0028437'
  # 'PlateId,P0028428,P0028439,P0028462,P0028437,P0028430,P0028448,P0028452,P0028447'
  # 'PlateId,P0028448,P0028451,P0028438,P0028565,P0028449,P0028452,P0028428,P0028447'
  # 'PlateId,P0028437,P0028438,P0028463,P0028430,P0028462,P0028439,P0028428,P0028447,P0028448'
)
covs=(
  # 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  # 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  # 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  # 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  # 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  # 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  # 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
)
apoe=("" ".APOE")
modules=("TMT-AD" "TMT-Resilience" "TREAT-AD" "Wan_etal")

# 02-linreg_lme4qtl.R ------------------------------------------
# Rscript 02-linreg_lme4qtl.R \
#   -i "../baseline_PCs/Plasma_pTau217.APOE/dataProt_SS-205063.PC.Plasma_pTau217.APOE.INTnorm.csv.gz"   \
#   -o "Plasma_pTau217.APOE" \
#   -p "../prot_names.txt" \
#   -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,APOE'   \
#   -d 'PlateId,P0028428,P0028439,P0028462,P0028437,P0028430,P0028448,P0028452,P0028447'  \
#   -c "Plasma_pTau217" \
#   -n "dataProt_SS-205063.PC.related.Plasma_pTau217.APOE" \
#   -r "../genomic/mypcrelate.Plasma_pTau217.rds" \

parallel -j 2 Rscript 02-linreg_lme4qtl.R \
  -i "../baseline_PCs/{1}{4}/dataProt_SS-205063.PC.{1}{4}.INTnorm.csv.gz"   \
  -o {1}{4}   \
  -p "../prot_names.txt" \
  -v {3}{=4 s/\./\,/ =}  \
  -d {2}  \
  -t 1 \
  -c "{1}" \
  -n "dataProt_SS-205063.PC.related.{1}{4}" \
  -r "../genomic/mypcrelate.{1}.rds" \
  ::: "${names[@]}" \
  :::+ "${recodes[@]}"  \
  :::+ "${covs[@]}"  \
  ::: "${apoe[@]}"

# 03-Protein_QC.R ------------------------------------------
# ----------
parallel -j 8 Rscript ../baseline/03-Protein_QC.R   \
  -i "{1}{2}/dataProt_SS-205063.PC.related.{1}{2}.lme4qtl.csv.gz"  \
  -o "{1}{2}" \
  -r "../2024-08-30_Somalogic_reference_Key.txt" \
  -n "dataProt_SS-205063.PC.related.{1}{2}" \
  ::: "${names[@]}" \
  ::: "${apoe[@]}"
# ----------

# for a in "${apoe[@]}"; do
#   for name in "${names[@]}"; do
#     Rscript ../baseline/03-Protein_QC.R \
#     -i "${name}${a}/dataProt_SS-205063.PC.related.${name}${a}.linreg_lme4qtl.csv.gz"  \
#     -o "${name}${a}" \
#     -r "../2024-08-30_Somalogic_reference_Key.txt" \
#     -n "dataProt_SS-205063.PC.related.${name}${a}"
#   done
# done

# Rscript 03-Protein_QC.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.linreg_glm.csv.gz"  \
#   -o "PlasmaAB140P" \
#   -r "../Somalogic_reference_Key.txt" \
#   -n "dataProt_SS-205063.PlasmaAB140P"

# 04-TREAT-AD_correlation.R --------------------------------
# ----------
parallel -j 6 Rscript ../baseline/04-TREAT-AD_correlation.R   \
  -i "{1}{2}/dataProt_SS-205063.PC.related.{1}{2}.csv"  \
  -o "{1}{2}" \
  -n "dataProt_SS-205063.PC.related.{1}{2}" \
  -c "{1}"  \
  ::: "${names[@]}" \
  ::: "${apoe[@]}"
# ----------

# for a in "${apoe[@]}"; do
#   for name in "${names[@]}"; do
#     Rscript ../baseline/04-TREAT-AD_correlation.R \
#     -i "${name}${a}/dataProt_SS-205063.PC.related.${name}${a}.csv"  \
#     -o "${name}${a}" \
#     -n "dataProt_SS-205063.PC.related.${name}${a}" \
#     -c "$name"
#   done
# done

# Rscript 04-TREAT-AD_correlation.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
#   -o "PlasmaAB140P" \
#   -n "dataProt_SS-205063.PlasmaAB140P" \
#   -c "PlasmaAB140P"

# 05-Combined-Enrich.R -------------------------------------
# ----------
parallel -j 6 Rscript ../baseline/05-Combined-Enrich.R  \
  -i "{1}{3}/dataProt_SS-205063.PC.related.{1}{3}.csv"  \
  -o "{1}{3}" \
  -n "dataProt_SS-205063.PC.related.{1}{3}" \
  -m "{2}" \
  ::: "${names[@]}" \
  ::: "${modules[@]}" \
  ::: "${apoe[@]}"
# ----------

# for a in "${apoe[@]}"; do
#   for name in "${names[@]}"; do
#     for mod in "${modules[@]}"; do
#     Rscript ../baseline/05-Combined-Enrich.R \
#     -i "${name}${a}/dataProt_SS-205063.PC.related.${name}${a}.csv"  \
#     -o "${name}${a}" \
#     -n "dataProt_SS-205063.PC.related.${name}${a}" \
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
# ----------
parallel -j 6 Rscript ../baseline/06-TREAT-AD_GSEA.R   \
  -i "{1}{2}/dataProt_SS-205063.PC.related.{1}{2}.csv"  \
  -o "{1}{2}" \
  -n "dataProt_SS-205063.PC.related.{1}{2}"  \
  ::: "${names[@]}" \
  ::: "${apoe[@]}"
# ----------

# for a in "${apoe[@]}"; do
#   for name in "${names[@]}"; do
#     Rscript ../baseline/06-TREAT-AD_GSEA.R \
#     -i "${name}${a}/dataProt_SS-205063.PC.related.${name}${a}.csv"  \
#     -o "${name}${a}" \
#     -n "dataProt_SS-205063.PC.related.${name}${a}"
#   done
# done

# Rscript ../baseline/06-TREAT-AD_GSEA.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PC.PlasmaAB140P.csv"  \
#   -o "PlasmaAB140P" \
#   -n "dataProt_SS-205063.PC.PlasmaAB140P" 

# 07-Sharma_etal_correlation.R -----------------------------
# ----------
parallel -j 6 Rscript ../baseline/07-Sharma_etal_correlation.R   \
  -i "{1}{2}/dataProt_SS-205063.PC.related.{1}{2}.csv"  \
  -o "{1}{2}" \
  -n "dataProt_SS-205063.PC.related.{1}{2}" \
  -c "{1}" \
  ::: "${names[@]}" \
  ::: "${apoe[@]}"
# ----------

# for a in "${apoe[@]}"; do
#   for name in "${names[@]}"; do
#     Rscript ../baseline/07-Sharma_etal_correlation.R \
#     -i "${name}${a}/dataProt_SS-205063.PC.related.${name}${a}.csv"  \
#     -o "${name}${a}" \
#     -n "dataProt_SS-205063.PC.related.${name}${a}" \
#     -c "$name" 
#   done
# done

# Rscript 07-Sharma_etal_correlation.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
#   -o "PlasmaAB140P" \
#   -n "dataProt_SS-205063.PlasmaAB140P" \
#   -c "PlasmaAB140P" 
