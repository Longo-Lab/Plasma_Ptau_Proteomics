#!/usr/bin/env bash

###################################################################
#Script Name  : 03-magma_command_top20.sh 
#Description  : Runs hmagma on top20 proteins in plasma proteomics
#Usage        : sbatch 03-magma_command_top20.sh
#Author       : Robert R Butler III
#Date Created : 2024-09-18
#Email        : rrbutler@stanford.edu
#Copyright (c) 2024
###################################################################
set -e

parallel="parallel --delay 2 -j 1 --link"

# where and which gene lists to use for *.geneset file
geneset_dir="."
geneset_prefix='top20'
genesets=('PlasmaGFAP' 'PlasmaNfL' 'PlasmaPTau181' 'PlasmapTau217' 'PlasmaAB140P' 'PlasmaAB142P' 'PlasmaABRatio' 'PlasmapTau217_AB42Ratio')

# where and which covars to use (if not in this dir, will check in geneset dir)
covar_dir='ref'
covars=('hla.filtered')

# which batchname to give this run
batchname=('top20')

# hmagma based annot selection (requires GRCh37)
ensg='ref/ENSG_coord_vGRCh37_PC.GRCh37.filtered.txt'
annotname='AdultBrain'

# make geneset file
# modified to check if geneset exists before making
if [ ! -f "${batchname}.geneset" ]; then
  # grabbing the genesets (flexible suffix), if not there, will try covar_dir
  for k in ${genesets[*]}; do
    if [ -f "${geneset_dir}/${geneset_prefix}_${k}.geneset" ]; then
      cp ${geneset_dir}/${geneset_prefix}_${k}.geneset ${batchname}_${k}.geneset
    elif [ -f "${geneset_dir}/${geneset_prefix}_${k}.txt" ]; then
      cp ${geneset_dir}/${geneset_prefix}_${k}.txt ${batchname}_${k}.geneset
    else
      cp ${covar_dir}/cov_${k}.geneset ${batchname}_${k}.geneset
    fi
  done
  # grabbing covariates, if not in covar_dir, will try geneset_dir
  for k in ${covars[*]}; do
    if [ -f "${covar_dir}/cov_${k}.geneset" ]; then
      cp ${covar_dir}/cov_${k}.geneset ${batchname}_${k}.geneset
    else
      cp ${geneset_dir}/cov_${k}.geneset ${batchname}_${k}.geneset
    fi
  done
  # need to make geneset file from pieces
  # rm ${batchname}.geneset
  for k in ${batchname}_*.geneset; do
    genes=$(cat $k)
    out=$(echo $k | sed "s/${batchname}_//" | sed "s/.geneset//")
    echo $out $genes >> ${batchname}.geneset
  done
fi


# running geneset only magma 
parallel -j 8 magma \
  --gene-results {1}_${annotname}/{1}_${annotname}.genes.raw \
  --set-annot ${batchname}.geneset \
  --model condition-hide=hla.filtered \
  --out {1}_${annotname}_${batchname} \
  ::: 'ALS_g1000' 'ALZ_g1000' 'FTD_g1000' 'LBD_g1000' 'PD_g1000' 'rAOO_g1000' 'BMI_g1000'

# cleanup genesets
rm -f ${batchname}_*.geneset

# squeue -j $SLURM_JOBID --Format=TimeUsed

# QC
Rscript 02b-magma_posthoc_QC.R
# Rscript magma_ranges.R ${batchname}
# Rscript --vanilla new_h_magma_top_scores.R
# Rscript --vanilla h_magma_heatmap.R

# output runtime
# squeue -j $SLURM_JOBID --Format=TimeUsed


