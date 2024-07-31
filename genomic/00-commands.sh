#!/usr/bin/env bash

# prepping data for plink
Rscript 01-genomic_cleaning.R

# plink reorganization ---------------------------------------------
wgs='/mnt/box/Stanford_ADRC/SRS_WGS'
lrs='/mnt/box/Stanford_ADRC/LRS_WGS'

# find the shared SNPs
if [ ! -f vars.txt ]; then
  awk 'NR==FNR { lines[$2]=1; next } $2 in lines {print $2}' \
    "${wgs}/srs_temp.bim" \
    "${lrs}/lrs_temp.bim" \
    > vars.txt
fi

# SRS sample filter
plink2 \
  --bfile $wgs/srs_temp \
  --make-bed \
  --remove $wgs/srs_remove.txt \
  --extract vars.txt \
  --snps-only 'just-acgt' \
  --out $wgs/srs_temp2

# LRS sample filter
plink2 \
  --bfile $lrs/lrs_temp \
  --make-bed \
  --remove $lrs/lrs_remove.txt \
  --extract vars.txt \
  --snps-only 'just-acgt' \
  --out $lrs/lrs_temp2

# output lrs Apoe types
printf "19:44908684:T:C\n19:44908822:C:T\n" > apoe_vars.txt
plink2 \
  --bfile $lrs/lrs_temp \
  --extract apoe_vars.txt \
  --recode A \
  --out $lrs/lrs_temp

# Join files
# printf "$wgs/srs_temp2\n$lrs/lrs_temp2\n" > bfiles.txt
plink \
  --bfile $wgs/srs_temp2 \
  --bmerge $lrs/lrs_temp2 \
  --out $wgs/all_temp

# filters by MAF > 0.02 (2%), Sample missingness (99%), sliding window LD
plink2 \
  --bfile $wgs/all_temp \
  --make-bed \
  --maf 0.02 \
  --mind 0.01 \
  --hwe midp 1e-05 \
  --indep-pairwise 200 100 0.2 \
  --out $wgs/srs-lrs_comb.filt

# only pruned
plink2 \
  --bfile $wgs/all_temp_filt \
  --make-bed \
  --extract $wgs/all_temp_filt.prune.in \
  --out $wgs/srs-lrs_comb.filt

# # apoe set # wgs seq dosn't have Apoe4 locus
# plink2 \
  # --bfile $wgs/all_temp \
  # --make-bed \
  # --extract apoe_vars.txt \
  # --out $wgs/srs-lrs_comb.apoe_type


# do some cleanup
ls $wgs/srs_temp*
ls $lrs/lrs_temp*
ls $wgs/all_temp*
rm $wgs/srs_temp*
rm $lrs/lrs_temp*
rm $wgs/all_temp*

Rscript 02-genomic_relatedness.R

