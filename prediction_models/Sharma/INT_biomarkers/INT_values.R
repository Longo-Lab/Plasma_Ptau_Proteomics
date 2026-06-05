install.packages("RNOmni")
library(RNOmni)

cols_to_norm <- c("PlasmaPTau181","PlasmaAB142P","PlasmaAB140P","PlasmaABRatio",          
                  "PlasmapTau217","PlasmapTau217_AB42Ratio","PlasmaGFAP","PlasmaNfL")

pheno2_norm <- test_data

pheno2_norm[paste0(cols_to_norm, "_norm")] <-
  lapply(pheno2_norm[cols_to_norm], function(x) RNOmni::RankNorm(as.numeric(x)))
pheno2_norm[cols_to_norm] <- NULL

#Did the transformation work correctly? Yes!!
nm <- cols_to_norm[1]
summary(pheno2[[nm]])
summary(pheno2_norm[[nm]])
summary(pheno2[[nm]])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.780   1.310   1.620   1.889   2.120   6.130 
summary(pheno2_norm[[nm]])                 # or [[paste0(nm,"_norm")]]
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.780   1.310   1.620   1.889   2.120   6.130 
cor(test_data[[nm]], pheno2_norm[[nm]], method="spearman", use="complete.obs")
#1

feature <- c(
  "PlasmaPTau181_norm", "PlasmaAB142P_norm","PlasmaAB140P_norm","PlasmaABRatio_norm", 
  "PlasmapTau217_norm","PlasmapTau217_AB42Ratio_norm","PlasmaGFAP_norm","PlasmaNfL_norm"
)

cols <- c("SampleID", outcome, covs, feature)

test_data <- pheno2_norm

output <- "INT_biomarkers"

# z-scale scores-------------------------------------------------------------------------
#scale on training data and use same scale for test data 

cols_to_z <- outcome

pheno2_z <- pheno

# Fit scaling on training data
scaled_block <- scale(pheno[, cols_to_z])

pheno2_z[, cols_to_z] <- scaled_block

# Save parameters for later
z_center <- attr(scaled_block, "scaled:center")
z_scale  <- attr(scaled_block, "scaled:scale")


#did it work?
apply(pheno2_z[, cols_to_z], 2, mean, na.rm = TRUE)
apply(pheno2_z[, cols_to_z], 2, sd,   na.rm = TRUE)


pheno2 <- pheno2_z


#need to apply this scaling to test_data

test_data_z <- test_data
test_data2 <- test_data

test_data_z[, cols_to_z] <-
  sweep(test_data[, cols_to_z], 2, z_center, "-")

test_data_z[, cols_to_z] <-
  sweep(test_data_z[, cols_to_z], 2, z_scale, "/")

test_data <- test_data_z
