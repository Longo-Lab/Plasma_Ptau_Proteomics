#!/usr/bin/env bash

###################################################################
#Script Name  : 01-magma_command_cont_prop.sh 
#Description  : Runs hmagma on continuous properties using 0kb window plasma proteomics
#Usage        : sbatch 08-magma_command_cont_prop.sh
#Author       : Robert R Butler III
#Date Created : 2024-09-18
#Email        : rrbutler@stanford.edu
#Copyright (c) 2024
###################################################################

# Condition for hla only --------------------------------------
grep -P '^hla' padj.geneset > covars.geneset

annotname='AdultBrain'

# running geneset only magma and then plotting results (two-sided)
parallel -j 8 magma \
  --gene-results {1}_${annotname}/{1}_${annotname}.genes.raw \
  --gene-covar {2}.txt \
  --set-annot covars.geneset \
  --model condition-hide=hla.filtered \
  --out {1}_${annotname}_{2} \
  ::: 'ALS' 'ALZ' 'FTD' 'LBD' 'PD' 'rAOO' 'BMI' \
  ::: 'betaLOGP' 't_stat' 'signLOGP'



