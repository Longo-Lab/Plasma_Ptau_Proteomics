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

modules=("TMT-AD" "TMT-Resilience" "TREAT-AD" "Wan_etal")

# # 01-sample_QC.R  ------------------------------------------
# for name in "${names[@]}"; do
#   if [["$name" == "PlasmaABRatio" || "$name" == "Plasma_pTau217_AB42Ratio"]]; then
#     covs='Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage'
#   else
#     covs='Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage'
#   fi

#   Rscript ../baseline/01-sample_QC.R \
#   -i "../dataProt_SS-205063.merged.QC.csv" \
#   -o "$name" \
#   -p "../prot_names.txt" \
#   -v "$covs" \
#   -c "$name" \
#   -n "dataProt_SS-205063.PC.$name" \
#   -m "../genomic/mypcair.PCs.csv"
# done

# # Rscript 01-sample_QC.R \
# #   -i "../dataProt_SS-205063.merged.QC.csv"   \
# #   -o "PlasmaAB140P"   \
# #   -p "../prot_names.txt" \
# #   -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage'   \
# #   -c "PlasmaAB140P"   \
# #   -n "dataProt_SS-205063.PlasmaAB140P"

# 02-linreg_lme4qtl.R ------------------------------------------
# test -------------------------------------
# parallel -j1 Rscript 02-linreg_lme4qtl.R \
  # -i "../baseline_PCs/PlasmaAB140P/dataProt_SS-205063.PC.PlasmaAB140P.INTnorm.csv.gz"   \
  # -o "PlasmaAB140P" \
  # -p "../prot_names_test.txt" \
  # -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'   \
  # -d 'PlateId,P0028451,P0028463,P0028462,P0028438,P0028430,P0028437'  \
  # -c "PlasmaAB140P" \
  # -n "dataProt_SS-205063.PC.related.PlasmaAB140P.t{}" \
  # -t {} \
  # -r "../genomic/mypcrelate.PlasmaAB140P.rds" \
  # ::: 1 2 4 8 16


parallel -j16 -k echo 02-linreg_lme4qtl.R \
  -i "../baseline_PCs/{1}{3}/dataProt_SS-205063.PC.{1}{3}.INTnorm.csv.gz"   \
  -o {1}{3}   \
  -p "../prot_names.txt" \
  -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'   \
  -d {2}  \
  -c "{1}"   \
  -n "dataProt_SS-205063.PC.related.{1}{3}" \
  -t 1 \
  -r "../genomic/mypcrelate.{1}.rds" \
  ::: "${names[@]}" \
  :::+ 'PlateId,P0028462,P0028452,P0028449,P0028428,P0028565,P0028463,P0028437'  \
       'PlateId,P0028451,P0028463,P0028462,P0028438,P0028430,P0028437' \
       'PlateId,P0028438,P0028439,P0028452,P0028449,P0028565,P0028451,P0028447,P0028430,P0028437' \
       'PlateId,P0028430,P0028439,P0028463,P0028438,P0028428,P0028447,P0028462,P0028448'  \
       'PlateId,P0028448,P0028438,P0028565,P0028437'  \
       'PlateId,P0028428,P0028439,P0028462,P0028437,P0028430,P0028448,P0028452,P0028447'  \
       'PlateId,P0028448,P0028451,P0028438,P0028565,P0028449,P0028452,P0028428,P0028447'  \
       'PlateId,P0028437,P0028438,P0028463,P0028430,P0028462,P0028439,P0028428,P0028447,P0028448'  \
  ::: '' '.APOE'

# ------------------------------------------


Rscript 02-linreg_lme4qtl.R \
  -i "../baseline_PCs/PlasmaPTau181/dataProt_SS-205063.PC.PlasmaPTau181.INTnorm.csv.gz"   \
  -o "PlasmaPTau181"   \
  -p "../prot_names.txt" \
  -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'   \
  -d 'PlateId,P0028462,P0028452,P0028449,P0028428,P0028565,P0028463,P0028437'  \
  -c "PlasmaPTau181"   \
  -n "dataProt_SS-205063.PC.related.PlasmaPTau181" \
  -r "../genomic/mypcrelate.PlasmaPTau181.rds"  

Rscript 02-linreg_lme4qtl.R \
  -i "../baseline_PCs/PlasmaAB140P/dataProt_SS-205063.PC.PlasmaAB140P.INTnorm.csv.gz"   \
  -o "PlasmaAB140P"   \
  -p "../prot_names.txt" \
  -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'   \
  -d 'PlateId,P0028451,P0028463,P0028462,P0028438,P0028430,P0028437'  \
  -c "PlasmaAB140P"   \
  -n "dataProt_SS-205063.PC.related.PlasmaAB140P" \
  -r "../genomic/mypcrelate.PlasmaAB140P.rds"

Rscript 02-linreg_lme4qtl.R \
  -i "../baseline_PCs/PlasmaAB142P/dataProt_SS-205063.PC.PlasmaAB142P.INTnorm.csv.gz"   \
  -o "PlasmaAB142P"   \
  -p "../prot_names.txt" \
  -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'   \
  -d 'PlateId,P0028438,P0028439,P0028452,P0028449,P0028565,P0028451,P0028447,P0028430,P0028437'  \
  -c "PlasmaAB142P"   \
  -n "dataProt_SS-205063.PC.related.PlasmaAB142P" \
  -r "../genomic/mypcrelate.PlasmaAB142P.rds"

Rscript 02-linreg_lme4qtl.R \
  -i "../baseline_PCs/Plasma_NFL/dataProt_SS-205063.PC.Plasma_NFL.INTnorm.csv.gz"   \
  -o "Plasma_NFL"   \
  -p "../prot_names.txt" \
  -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'   \
  -d 'PlateId,P0028430,P0028439,P0028463,P0028438,P0028428,P0028447,P0028462,P0028448'  \
  -c "Plasma_NFL"   \
  -n "dataProt_SS-205063.PC.related.Plasma_NFL" \
  -r "../genomic/mypcrelate.Plasma_NFL.rds"

Rscript 02-linreg_lme4qtl.R \
  -i "../baseline_PCs/Plasma_GFAP/dataProt_SS-205063.PC.Plasma_GFAP.INTnorm.csv.gz"   \
  -o "Plasma_GFAP"   \
  -p "../prot_names.txt" \
  -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'   \
  -d 'PlateId,P0028448,P0028438,P0028565,P0028437'  \
  -c "Plasma_GFAP"   \
  -n "dataProt_SS-205063.PC.related.Plasma_GFAP" \
  -r "../genomic/mypcrelate.Plasma_GFAP.rds"

Rscript 02-linreg_lme4qtl.R \
  -i "../baseline_PCs/Plasma_pTau217/dataProt_SS-205063.PC.Plasma_pTau217.INTnorm.csv.gz"   \
  -o "Plasma_pTau217"   \
  -p "../prot_names.txt" \
  -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'   \
  -d 'PlateId,P0028428,P0028439,P0028462,P0028437,P0028430,P0028448,P0028452,P0028447'  \
  -c "Plasma_pTau217"   \
  -n "dataProt_SS-205063.PC.related.Plasma_pTau217" \
  -r "../genomic/mypcrelate.Plasma_pTau217.rds"
 
 Rscript 02-linreg_lme4qtl.R \
  -i "../baseline_PCs/PlasmaABRatio/dataProt_SS-205063.PC.PlasmaABRatio.INTnorm.csv.gz"   \
  -o "PlasmaABRatio"   \
  -p "../prot_names.txt" \
  -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'   \
  -d 'PlateId,P0028448,P0028451,P0028438,P0028565,P0028449,P0028452,P0028428,P0028447'  \
  -c "PlasmaABRatio"   \
  -n "dataProt_SS-205063.PC.related.PlasmaABRatio" \
  -r "../genomic/mypcrelate.PlasmaABRatio.rds"

Rscript 02-linreg_lme4qtl.R \
  -i "../baseline_PCs/Plasma_pTau217_AB42Ratio/dataProt_SS-205063.PC.Plasma_pTau217_AB42Ratio.INTnorm.csv.gz"   \
  -o "Plasma_pTau217_AB42Ratio"   \
  -p "../prot_names.txt" \
  -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'   \
  -d 'PlateId,P0028437,P0028438,P0028463,P0028430,P0028462,P0028439,P0028428,P0028447,P0028448'  \
  -c "Plasma_pTau217_AB42Ratio"   \
  -n "dataProt_SS-205063.PC.related.Plasma_pTau217_AB42Ratio" \
  -r "../genomic/mypcrelate.Plasma_pTau217_AB42Ratio.rds"

# 03-Protein_QC.R ------------------------------------------
for name in "${names[@]}"; do
  Rscript ../baseline/03-Protein_QC.R \
  -i "$name/dataProt_SS-205063.PC.related.$name.lme4qtl.csv.gz"  \
  -o "$name" \
  -r "../Somalogic_reference_Key.txt" \
  -n "dataProt_SS-205063.PC.related.$name"
done

# Rscript 03-Protein_QC.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.linreg_glm.csv.gz"  \
#   -o "PlasmaAB140P" \
#   -r "../Somalogic_reference_Key.txt" \
#   -n "dataProt_SS-205063.PlasmaAB140P"

# 04-TREAT-AD_correlation.R --------------------------------
for name in "${names[@]}"; do
  Rscript ../baseline/04-TREAT-AD_correlation.R \
  -i "$name/dataProt_SS-205063.PC.related.$name.csv"  \
  -o "$name" \
  -n "dataProt_SS-205063.PC.related.$name" \
  -c "$name"
done

# Rscript 04-TREAT-AD_correlation.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
#   -o "PlasmaAB140P" \
#   -n "dataProt_SS-205063.PlasmaAB140P" \
#   -c "PlasmaAB140P"

# 05-Combined-Enrich.R -------------------------------------
for name in "${names[@]}"; do
  for mod in "${modules[@]}"; do
  Rscript ../baseline/05-Combined-Enrich.R \
  -i "$name/dataProt_SS-205063.PC.related.$name.csv"  \
  -o "$name" \
  -n "dataProt_SS-205063.PC.related.$name" \
  -m "$mod"
  done
done

# Rscript 05-Combined-Enrich.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
#   -o "PlasmaAB140P" \
#   -n "dataProt_SS-205063.PlasmaAB140P" \
#   -m "TMT-AD"

# 06-TREAT-AD_GSEA.R ---------------------------------------
for name in "${names[@]}"; do
  Rscript ../baseline/06-TREAT-AD_GSEA.R \
  -i "$name/dataProt_SS-205063.PC.related.$name.csv"  \
  -o "$name" \
  -n "dataProt_SS-205063.PC.related.$name"
done

# Rscript 06-TREAT-AD_GSEA.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
#   -o "PlasmaAB140P" \
#   -n "dataProt_SS-205063.PlasmaAB140P" 

# 07-Sharma_etal_correlation.R -----------------------------
for name in "${names[@]}"; do
  Rscript ../baseline/07-Sharma_etal_correlation.R \
  -i "$name/dataProt_SS-205063.PC.related.$name.csv"  \
  -o "$name" \
  -n "dataProt_SS-205063.PC.related.$name" \
  -c "$name" 
done

# Rscript 07-Sharma_etal_correlation.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
#   -o "PlasmaAB140P" \
#   -n "dataProt_SS-205063.PlasmaAB140P" \
#   -c "PlasmaAB140P" 
