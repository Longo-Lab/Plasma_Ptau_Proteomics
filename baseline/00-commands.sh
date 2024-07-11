#!/usr/bin/env bash

Rscript 01-sample_QC.R

# # sbatch Rscript on the fly with --wrap
# sbatch -J ptau --mem=12G -c 16 -t 06:00:00 -p interactive \
  # -o %A_02-linreg_lme4.log \
  # --wrap "ml R/4.2.2; Rscript 02-linreg_lme4.R"
# run local faster
Rscript 02-linreg_glm.R

Rscript 03-Protein_QC.R

Rscript 04-TREAT-AD_correlation.R
Rscript 05-TREAT-AD_enrich.R
Rscript 06-Milind_etal.R
Rscript 07-Mostafavi_etal.R
Rscript 08-TMT-AD.R
Rscript 09-Wan_etal.R
Rscript 10-TREAT-AD_GSEA.R
Rscript 11-Sharma_etal_correlation.R

# demoing new overlap enrichment
mkdir -p test
cd test
Rscript ../12-Combined-Enrich.R -i ../dataProt_SS-205063.pTau.baseline.csv -m TREAT-AD -r "../../biodomains"
Rscript ../12-Combined-Enrich.R -i ../dataProt_SS-205063.pTau.baseline.csv -m TMT-AD -r "../../biodomains"
Rscript ../12-Combined-Enrich.R -i ../dataProt_SS-205063.pTau.baseline.csv -m Wan_etal -r "../../biodomains"
Rscript ../12-Combined-Enrich.R -i ../dataProt_SS-205063.pTau.baseline.csv -m TMT-Resilience -r "../../biodomains"