#!/usr/bin/env bash

###################################################################
#Script Name  : 00-commands.sh 
#Description  : Main association analysis wrapper script
#Usage        : bash 00-commands.sh
#Author       : Robert R Butler III
#Date Created : 2026-03-03
#Email        : rrbutler@stanford.edu
#Copyright (c) 2026
###################################################################

# Notes       :
# must have run scripts in the parent folder to get dataProt_SS-205063.merged.QC.csv
# 
# now kill yourself cause you have to deal with Matlab...  a god damned license, what losers

# # Clone this repository (may just be me but had write issues)
# sudo git clone https://github.com/WashUFluidBiomarkers/plasma_ptau217_time.git
# sudo git clone https://github.com/Betthauser-Neuro-Lab/SILA-AD-Biomarker.git
# cp plasma_ptau217_time/sila_code.m sila_code_pTau217-AB42.m

# SILA model had best pTau217/AB42 R^2 and correlation.
Rscript 01-sira_age_onset_estimation.R
Rscript 01-sira_age_onset_estimation.R -f "cooksd" -n "dataProt_SS-205063.SILA_cooksd"

