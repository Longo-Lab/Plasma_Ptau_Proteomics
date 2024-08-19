#!/usr/bin/env bash

# # # # sbatch Rscript on the fly with --wrap
# # # sbatch -J ptau --mem=12G -c 16 -t 06:00:00 -p interactive \
  # # # -o %A_02-linreg_lme4.log \
  # # # --wrap "ml R/4.2.2; Rscript 02-linreg_lme4.R"
# # # run local faster

# # Rscript 06-Milind_etal.R
# # Rscript 07-Mostafavi_etal.R

# demoing new overlap enrichment
# mkdir -p test
# cd test
# Rscript ../12-Combined-Enrich.R \
#   -i ../dataProt_SS-205063.pTau.baseline.csv \
#   -m TREAT-AD \
#   -r "../../biodomains"
# Rscript ../12-Combined-Enrich.R \
#   -i ../dataProt_SS-205063.pTau.baseline.csv \
#   -m TMT-AD \
#   -r "../../biodomains"
# Rscript ../12-Combined-Enrich.R \
#   -i ../dataProt_SS-205063.pTau.baseline.csv \
#   -m Wan_etal \
#   -r "../../biodomains"
# Rscript ../12-Combined-Enrich.R \
#   -i ../dataProt_SS-205063.pTau.baseline.csv \
#   -m TMT-Resilience \
#   -r "../../biodomains"


names=(
  "PlasmaPTau181"
  "PlasmaAB140P" 
  "PlasmaAB142P" 
  "PlasmaABRatio" 
  "Plasma_NFL" 
  "Plasma_GFAP" 
  "Plasma_pTau217" 
  "Plasma_pTau217_AB42Ratio"
)

modules=("TMT-AD" "TMT-Resilience" "TREAT-AD" "Wan_etal")

# 01-sample_QC.R  ------------------------------------------
for name in "${names[@]}"; do
  Rscript 01-sample_QC.R \
  -i "../dataProt_SS-205063.merged.QC.csv" \
  -o "$name" \
  -p "../prot_names.txt" \
  -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage' \
  -c "$name" \
  -n "dataProt_SS-205063.$name"
done

# # Rscript 01-sample_QC.R \
# #   -i "../dataProt_SS-205063.merged.csv"   \
# #   -o "PlasmaAB140P"   \
# #   -p "../prot_names.txt" \
# #   -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage'   \
# #   -c "PlasmaAB140P"   \
# #   -n "dataProt_SS-205063.PlasmaAB140P"

# # 02-$name-linreg_glm.R ------------------------------------------
# ### This is for running biomarker specific linreg scripts
# for name in "${names[@]}"; do
#   Rscript $name/02-$name-linreg_glm.R \
#   -i "$name/dataProt_SS-205063.$name.INTnorm.csv.gz"   \
#   -o "$name"   \
#   -p "../prot_names.txt" \
#   -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage'   \
#   -c "$name"   \
#   -n "dataProt_SS-205063.$name"
# done

# # Rscript 02-linreg_glm.R \
# #   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.INTnorm.csv.gz"   \
# #   -o "PlasmaAB140P"   \
# #   -p "../prot_names.txt" \
# #   -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage'   \
# #   -c "PlasmaAB140P"   \
# #   -n "dataProt_SS-205063.PlasmaAB140P"

# # 03-Protein_QC.R ------------------------------------------
# for name in "${names[@]}"; do
#   Rscript 03-Protein_QC.R \
#   -i "$name/dataProt_SS-205063.$name.linreg_glm.csv.gz"  \
#   -o "$name" \
#   -r "../Somalogic_reference_Key.txt" \
#   -n "dataProt_SS-205063.$name"
# done

# # Rscript 03-Protein_QC.R \
# #   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.linreg_glm.csv.gz"  \
# #   -o "PlasmaAB140P" \
# #   -r "../Somalogic_reference_Key.txt" \
# #   -n "dataProt_SS-205063.PlasmaAB140P"

# # 04-TREAT-AD_correlation.R --------------------------------
# for name in "${names[@]}"; do
#   Rscript 04-TREAT-AD_correlation.R \
#   -i "$name/dataProt_SS-205063.$name.csv"  \
#   -o "$name" \
#   -n "dataProt_SS-205063.$name" \
#   -c "$name"
# done

# # Rscript 04-TREAT-AD_correlation.R \
# #   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
# #   -o "PlasmaAB140P" \
# #   -n "dataProt_SS-205063.PlasmaAB140P" \
# #   -c "PlasmaAB140P"

# # 05-Combined-Enrich.R -------------------------------------
# for name in "${names[@]}"; do
#   for mod in "${modules[@]}"; do
#   Rscript 05-Combined-Enrich.R \
#   -i "$name/dataProt_SS-205063.$name.csv"  \
#   -o "$name" \
#   -n "dataProt_SS-205063.$name" \
#   -m "$mod"
#   done
# done

# # Rscript 05-Combined-Enrich.R \
# #   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
# #   -o "PlasmaAB140P" \
# #   -n "dataProt_SS-205063.PlasmaAB140P" \
# #   -m "TMT-AD"

# # 06-TREAT-AD_GSEA.R ---------------------------------------
# for name in "${names[@]}"; do
#   Rscript 06-TREAT-AD_GSEA.R \
#   -i "$name/dataProt_SS-205063.$name.csv"  \
#   -o "$name" \
#   -n "dataProt_SS-205063.$name"
# done

# # Rscript 06-TREAT-AD_GSEA.R \
# #   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
# #   -o "PlasmaAB140P" \
# #   -n "dataProt_SS-205063.PlasmaAB140P" 

# # 07-Sharma_etal_correlation.R -----------------------------
# for name in "${names[@]}"; do
#   Rscript 07-Sharma_etal_correlation.R \
#   -i "$name/dataProt_SS-205063.$name.csv"  \
#   -o "$name" \
#   -n "dataProt_SS-205063.$name" \
#   -c "$name" 
# done

# # Rscript 07-Sharma_etal_correlation.R \
# #   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
# #   -o "PlasmaAB140P" \
# #   -n "dataProt_SS-205063.PlasmaAB140P" \
# #   -c "PlasmaAB140P" 
