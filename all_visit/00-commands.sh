#!/usr/bin/env bash

Rscript 01-sample_QC.R

# # sbatch Rscript on the fly with --wrap
# sbatch -J ptau --mem=12G -c 16 -t 06:00:00 -p interactive \
  # -o %A_02-linreg_lme4.log \
  # --wrap "ml R/4.2.2; Rscript 02-linreg_lme4.R"
# run local faster
Rscript 02-linreg_glm.R

Rscript 03-Protein_QC.R
