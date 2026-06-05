#!/usr/bin/env bash

# # run on scg
# #SBATCH --job-name=mv_hmagma
# #SBATCH --nodes=1
# #SBATCH --ntasks-per-node=1
# #SBATCH --cpus-per-task=9
# #SBATCH --partition=interactive
# #SBATCH --account=default
# #SBATCH --time=24:00:00
# #SBATCH --mem=64G
# #SBATCH -o %A_hmagma_stim.log

# # load modules
# ml parallel R/4.2.2 python/3.7.8
# set -ex
set -e

# # Instead of $SLURM_NTASKS, use $SLURM_NNODES to determine how
# # many jobs should be run simultaneously.
# srun="srun -N1 -n1 --cpus-per-task $SLURM_CPUS_PER_TASK"
parallel="parallel --delay 2 -j 1 --link"
# root_dir='/labs/flongo/PS19_C31_stim/Gene_level/DGE'

# where and which gene lists to use for *.geneset file
geneset_dir="."
geneset_prefix='padj'
genesets=('PlasmaGFAP' 'PlasmaNfL' 'PlasmaPTau181' 'PlasmapTau217' 'PlasmaAB140P' 'PlasmaAB142P' 'PlasmaABRatio' 'PlasmapTau217_AB42Ratio')

# where and which covars to use (if not in this dir, will check in geneset dir)
covar_dir='ref'
covars=('hla.filtered')

# which batchname to give this run
batchname=('padj')

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
fi

# making annots
awk 'NR==FNR{A[$1];next}($1 in A)' $ensg ref/${annotname}.genes.annot \
  > ${annotname}_temp_annot
parallel -j 7 cp ${annotname}_temp_annot \
  {}_${annotname}.genes.annot \
  ::: 'ALS_g1000' 'ALZ_g1000' 'FTD_g1000' 'LBD_g1000' 'PD_g1000' 'rAOO_g1000' 'BMI_g1000'
rm ${annotname}_temp_annot

# running magma and then plotting results
$parallel python3 02a-magma_geneset.py \
  --threads 0 \
  -u {1} -d {2} \
  --grch37 \
  -p \
  --gene_model 'multi=snp-wise' \
  -a ${annotname} \
  -b ${batchname} \
  -c ${covars[*]} \
  -o magma_geneset_${batchname}_MAGMA \
  --gwas-dir 'ref' \
  --gwas 'ALS_g1000' 'ALZ_g1000' 'FTD_g1000' 'LBD_g1000' 'PD_g1000' 'rAOO_g1000' 'BMI_g1000' \
  --gref 'ref/g1000_eur/g1000_eur' \
  ::: '0' \
  ::: '0'

# Collect gene level data in folder
for k in *.genes.raw
  do
    out=$(echo $k | sed 's/\.genes\.raw//')
    [ ! -d "${out}" ] && mkdir -p ${out}
    mv ${out}\.[gl]* ${out}/
done

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
