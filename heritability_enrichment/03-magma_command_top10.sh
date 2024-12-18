#!/usr/bin/env bash

###################################################################
#Script Name  : 03-magma_command_top10.sh 
#Description  : Runs hmagma on top10 proteins in plasma proteomics
#Usage        : sbatch 03-magma_command_top10.sh
#Author       : Robert R Butler III
#Date Created : 2024-09-18
#Email        : rrbutler@stanford.edu
#Copyright (c) 2024
###################################################################

parallel="parallel --delay 2 -j 1 --link"

# where and which gene lists to use for *.geneset file
geneset_dir="."
geneset_prefix='top10'
genesets=('Plasma_GFAP' 'Plasma_NFL' 'PlasmaPTau181' 'Plasma_pTau217' 'PlasmaAB140P' 'PlasmaAB142P' 'PlasmaABRatio' 'Plasma_pTau217_AB42Ratio')

# where and which covars to use (if not in this dir, will check in geneset dir)
covar_dir='ref'
covars=('hla.filtered')

# which batchname to give this run
batchname=('top10')

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
  ::: 'ALS' 'ALZ' 'FTD' 'LBD' 'PD' 'rAOO' 'BMI'

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


