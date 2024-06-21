#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
library(ggridges)
library(patchwork)
library(RNOmni)

options(width=100)
set.seed(1234) # setting the seed so I pick the same 6 random samples for all

# load the samples unscaled to external reference set --------------------------
a = fread("anml.QC/dataProt_SS-205063.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck_ADRC_Feb2021.csv.csv")

# dropping non-samples
a = a[SampleType == "Sample"]
length(a)

# checking minimum value of all spots 
table(a[,4:length(a)] <= 1)
## If it is less than 1, log transform will produce negative values
## Turns out there will be some negatives values
## Using log transform adding a pseudo-count is advised

# transposing 6 samples for easy viewing, and dropping the "V1" and "SampleType" columns
rows = sample(nrow(a), 6)
sprintf("Grabbing random samples: %s for testing...", paste(a[rows, SampleId], collapse=", "))
raw = transpose(a[rows, c(2, 4:7599)], keep.names="prot", make.names="SampleId")

# load the samples scaled to external reference --------------------------------
a = fread("anml.QC.anml.SMP/dataProt_SS-205063.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP_ADRC_Feb2021.csv")

# dropping non-samples
a = a[SampleType == "Sample"]
length(a)

# transposing 6 samples for easy viewing, and dropping the "V1" and "SampleType" columns
sprintf("Grabbing random samples: %s for testing...", paste(a[rows, SampleId], collapse=", "))
scaled = transpose(a[rows, c(2, 4:7599)], keep.names="prot", make.names="SampleId")

# check the scaling of the two -------------------------------------------------
sapply(scaled[,2:6], sum)
sapply(raw[,2:6], sum)
sapply(scaled[,2:6], summary)
sapply(raw[,2:6], summary)
## medians for the scaled set are at approximately 1000
## sum RFU for all spots per each sample is also noticeably scaled to the 30-50 mil range
## definitely want to use the scaled 
## not sure if it is possible to get any better scaling since these are intensities not counts
## should have a dramatic impact on inflated correlations

# Normalized values ------------------------------------------------------------
## log scaling base is really more of a preference. 
## log10 and log2 the shape is identical, only the numbers on the axis change
## but in log2, an increase of 1 represents a doubling, log10, a ten-fold increase
## generally if observing biological shifts like protein, a two-fold is more likely
cols = names(scaled)[-1]

# log 2 transform (add 1 because some values less than 1)
scaled2 = copy(scaled)
scaled2[, (cols) := lapply(.SD, function(i) log2(i + 1)), .SDcols=cols]

# log 2 transform (add 1 because some values less than 1)
raw2 = copy(raw)
raw2[, (cols) := lapply(.SD, function(i) log2(i + 1)), .SDcols=cols]

# check the normality/convenience of the two vs ten log scales
sapply(scaled2[,2:6], summary)
sapply(raw2[,2:6], summary)

# Plotting the distrubutions ---------------------------------------------------
# convert the tables for ridge plots 
Scaled = melt(scaled2, id.vars="prot", variable.name="SampleId", value.name="Log2RFU")
Unscaled = melt(raw2, id.vars="prot", variable.name="SampleId", value.name="Log2RFU")

# plot 
p.list = lapply(c("Scaled", "Unscaled"), function(i) {
  ggplot(get(i), aes(x=Log2RFU, y=SampleId, fill=SampleId)) +
    geom_density_ridges() + 
    ggtitle(sprintf("%s to external reference SomaLogic RFUs (log transformed)", i)) +
    theme_classic() + theme(plot.title=element_text(hjust=0.5))
})
wrap_plots(p.list)
## Okay, so clear difference with the scaling to external, will help a lot!!!
## Note that even with the Log transform though, there is a positive skew
## * May be a good opportunity to try inverse rank normalization
## this will potentially still cause issues, but luckily we can row-scale each 
## protein to get a z-score measure 


# Final scaling of data for lmm regression -------------------------------------
# drop unnecessary rows
a = fread("anml.QC.anml.SMP/dataProt_SS-205063.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP_ADRC_Feb2021.csv")
a = a[SampleType == "Sample"]
a[, c("V1", "SampleType") := NULL]
cols = names(a)[-1]

# log 2 transform of external scaled (add 1 because some values less than 1)
a[, (cols) := lapply(.SD, function(i) log2(i + 1)), .SDcols=cols]

# scale each gene column to z-scale
a[, (cols) := lapply(.SD, scale), .SDcols=cols]
## not normally distributed via shapiro-wilk, nor plots, though some approaching 
lapply(a[, 2:5], shapiro.test)
protz = melt(a[, 1:5], id.vars="SampleId", variable.name="Prot", value.name="z-scaled Log2RFU")
p1 = ggplot(protz, aes(x=`z-scaled Log2RFU`, y=Prot, fill=Prot)) +
  geom_density_ridges() + xlim(-4, 4) + 
  ggtitle("Z-scaled proteins across samples") +
  theme_classic() + theme(plot.title=element_text(hjust=0.5))

## may have to utilize inverse rank normalization ------------------------------
a = fread("anml.QC.anml.SMP/dataProt_SS-205063.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP_ADRC_Feb2021.csv")
a = a[SampleType == "Sample"]
a[, c("V1", "SampleType") := NULL]
# a[, (cols) := lapply(.SD, function(i) log2(i + 1)), .SDcols=cols]

# Inverse rank normalization (INT) to z-scale
a[, (cols) := lapply(.SD, RankNorm), .SDcols=cols]
## are normally distributed via shapiro-wilk, and plots Bingo
lapply(a[, 2:5], shapiro.test)
protINT = melt(a[, 1:5], id.vars="SampleId", variable.name="Prot", value.name="INT Log2RFU")
p2 = ggplot(protINT, aes(x=`INT Log2RFU`, y=Prot, fill=Prot)) +
  geom_density_ridges() + xlim(-4, 4) + 
  ggtitle("Rank-based INT scaled proteins across samples") +
  theme_classic() + theme(plot.title=element_text(hjust=0.5))
p1 + p2

###########################IMPORTANT#####################################
## In fact, since it is rank-based, log transform is totally unnecessary, 
## values are exactly the same either way as long as it is the scaled dataset.
## Rank-based INT PROTEIN VALUES CANNOT BE COMPARED WITHIN SAMPLE, ONLY ACROSS SAMPLES
## similarly COMPARING DIFFERENT PROTEINS IN THE SAME SAMPLE MUST USE SCALED LOG2(RFU + 1)
#########################################################################
# So, the final INTnorm dataset:
a = fread("anml.QC.anml.SMP/dataProt_SS-205063.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP_ADRC_Feb2021.csv")
a = a[SampleType == "Sample"]
a[, c("V1", "SampleType") := NULL]

# Inverse rank normalization (INT) to z-scale
a[, (cols) := lapply(.SD, RankNorm), .SDcols=cols]

# fwrite(a, file="dataProt_SS-205063.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.INTnorm.across_samples_only_ADRC_Mar2023.csv")

# The final Log2RFU dataset:
a = fread("anml.QC.anml.SMP/dataProt_SS-205063.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP_ADRC_Feb2021.csv")
a = a[SampleType == "Sample"]
a[, c("V1", "SampleType") := NULL]

# log 2 transform of external scaled (add 1 because some values less than 1)
a[, (cols) := lapply(.SD, function(i) log2(i + 1)), .SDcols=cols]

# fwrite(a, file="dataProt_SS-205063.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.Log2RFU.with_sample_only_ADRC_Mar2023.csv")
